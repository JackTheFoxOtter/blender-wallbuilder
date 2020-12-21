# ##### BEGIN GPL LICENSE BLOCK #####
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# ##### END GPL LICENSE BLOCK #####

bl_info = {
    "name": "WallBuilder",
    "description": "Utilities to support wall and room building.",
    "author": "JackTheFoxOtter",
    "version": (1, 2),
    "blender": (2, 90, 1),
    "location": "View3D > Object > WallBuilder",
    "warning": "",
    "wiki_url": "https://github.com/JackTheFoxOtter/blender-wallbuilder",
    "tracker_url": "https://github.com/JackTheFoxOtter/blender-wallbuilder/issues",
    "support": "COMMUNITY",
    "category": "Object",
}

"""
Utilities to support wall and room building.
"""

from bpy.utils import register_class, unregister_class
from bpy.props import FloatProperty, BoolProperty
from functools import cmp_to_key
from mathutils import Vector
import colorsys
import bmesh
import math
import bpy

def vec_to_str(vec):
    """
    Returns a formatted string for the given vector.
    """
    return f"<{' '.join([str(x) for x in vec])}>"


def get_signed_angle(v, vBase=Vector([0.0, 1.0])):
    """
    Returns angle of 2D-Vector v relative to 2D-Vector vBase.
    Defaults to base X=0.0, Y=1.0 if vBase isn't specified.
    Raises an exception if angle couldn't be determined
    """
    v = v.to_2d()
    vBase = vBase.to_2d()
    
    angle = vBase.angle_signed(v, None)
    
    if angle is None:
        raise Exception(f"Couldn't determine angle for vector {vec_to_str(v)} and base {vec_to_str(vBase)}!")
        
    return angle


def get_adjesant_edge(origin_edge, origin_vert, inverted=False):
    """
    Returns edge connected to origin_vert which's angle is closest to origin_edge. (Except origin_edge)
    Returns None if no further edges are connected to origin_vert.
    The parameter inverted specifies the direction of origin_edge, and therefore the role of origin_vert:
        False - origin_vert is the start vertex of origin_edge
        True - origin_vert is the end vertex of origin_edge
    The direction of the returned adjesant edge follows the direction of origin_edge.
    """
    adjesant_edges = [edge for edge in origin_vert.link_edges if edge != origin_edge]
    if len(adjesant_edges) == 0: return None # No adjesant edges
    if len(adjesant_edges) == 1: return adjesant_edges[0] # Only one edge, return that

    if inverted:
        origin_edge_direction = (origin_vert.co - origin_edge.other_vert(origin_vert).co).to_2d().normalized()
    else:
        origin_edge_direction = (origin_edge.other_vert(origin_vert).co - origin_vert.co).to_2d().normalized()

    # Determine directions of all connected edges
    edge_directions = []
    for edge in adjesant_edges:
        if inverted:
            direction = (edge.other_vert(origin_vert).co - origin_vert.co).to_2d().normalized()
        else:
            direction = (origin_vert.co - edge.other_vert(origin_vert).co).to_2d().normalized()
            
        edge_directions.append((edge, direction))

    # Sort connected edges by angle
    edge_directions.sort(key=lambda x: get_signed_angle(x[1], origin_edge_direction), reverse=True)
    
    return edge_directions[0][0] if inverted else edge_directions[len(edge_directions)-1][0]


def get_corner_position(origin, edge1, edge2, wall_thickness):
    """
    Returns the position of the corner of edge1 and edge2, connected by origin.
    Corner is offset of each edge as defined by wall_thickness.
    P = O ± v; v = ||a||*l - ||b||*l; l = 1/sin(α)
    """
    # Determine direction vectors of edges
    dir1 = (origin.co - edge1.other_vert(origin).co).normalized()
    dir2 = (edge2.other_vert(origin).co - origin.co).normalized()
    
    # Determine angle in between edges
    angle = get_signed_angle(dir1, dir2)
    if abs(angle) < 0.0001:
        # Edges are parallel
        orthogonal = Vector([dir1.y, -dir1.x, 0.0])
        return origin.co - orthogonal * wall_thickness
    
    # Calculate length of edge vector
    length = 1 / math.sin(angle)
    offset = (dir2 - dir1) * length * wall_thickness
    return origin.co + offset


def generate_wall_mesh_data(reference_obj, mesh, wall_thickness=0.125, fill_rims=True):
    """
    Generates a wall geometry based on the reference object and applies it to the specified target mesh.
    The reference object should contain a mesh with a 2D wireframe representation of the wall layout, without faces.
    """
    data = reference_obj.data
    bm = bmesh.new()
    bm.from_mesh(data)

    # Dynamic variables
    new_verts = []
    new_edges = []
    new_faces = []
    corner_vertices = {}

    # Generate mesh data for walls
    for edge in bm.edges:
        # Iterate through each edge in the mesh
        end_vertex_indices = []
        
        for start, end in [(edge.verts[0], edge.verts[1]), (edge.verts[1], edge.verts[0])]:
            # Runs twice per edge, for both possible directions
            start_edge = get_adjesant_edge(edge, start, False)
            end_edge = get_adjesant_edge(edge, end, True)
            
            if start_edge:
                # Edge exists
                corner_position = get_corner_position(start, edge, start_edge, wall_thickness).to_3d()
                corner_position.z = start.co.z
                new_verts.append(corner_position.to_tuple())
                # Add corner to corner_vertices list to create corner caps later
                corner_vertices.setdefault(start, []).append(corner_position)
            else:
                # Edge doesn't exist -> end segment
                direction = (edge.other_vert(start).co - start.co).normalized()
                orthogonal = Vector([direction.y, -direction.x, 0.0])
                corner_position = (start.co + orthogonal * wall_thickness).to_3d()
                corner_position.z = start.co.z
                new_verts.append(corner_position.to_tuple())
                # Store index of added end vertex so we can close it later
                end_vertex_indices.append(len(new_verts)-1)
            
            if end_edge:
                # Edge exists
                corner_position = get_corner_position(end, end_edge, edge, wall_thickness).to_3d()
                corner_position.z = end.co.z
                new_verts.append(corner_position.to_tuple())
            else:
                # Edge doesn't exist -> end segment
                direction = (end.co - edge.other_vert(end).co).normalized()
                orthogonal = Vector([direction.y, -direction.x, 0.0])
                corner_position = (end.co + orthogonal * wall_thickness).to_3d()
                corner_position.z = end.co.z
                new_verts.append(corner_position.to_tuple())
                # Store index of added end vertex so we can close it later
                end_vertex_indices.append(len(new_verts)-1)
        
        if fill_rims:
            # Fill the face between the created vertices
            new_faces.append(list(range(len(new_verts)-4, len(new_verts))))
        else:
            # Connect the created vertices through edges
            new_edges.append([len(new_verts)-4, len(new_verts)-3])
            new_edges.append([len(new_verts)-2, len(new_verts)-1])
            # Connect vertices of end segments
            if len(end_vertex_indices) == 2:
                # One end segment
                new_edges.append([end_vertex_indices[0], end_vertex_indices[1]])
            elif len(end_vertex_indices) == 4:
                # Two end segments
                new_edges.append([end_vertex_indices[0], end_vertex_indices[3]])
                new_edges.append([end_vertex_indices[1], end_vertex_indices[2]])

    if fill_rims:
        # Fill corner caps
        for corner, verts in corner_vertices.items():
            if len(verts) > 2:
                # Corner needs cap (face between vertices)
                # Sort vertices around corner_position
                # We do this to ensure the faces aren't messed up later
                vert_directions = []
                for vert in verts:
                    direction = (corner.co - vert).to_2d().normalized()
                    vert_directions.append((vert, direction))
                vert_directions.sort(key=lambda x: get_signed_angle(x[1]), reverse=True)
                
                # Append sorted vertices to new_verts for new_mesh
                for vert, direction in vert_directions:
                    new_verts.append(vert.to_3d().to_tuple())
                
                # Create corner cap faces
                new_faces.append(list(range(len(new_verts)-len(verts), len(new_verts))))
        
    # Fill new mesh with generated mesh data
    mesh.from_pydata(new_verts, new_edges, new_faces)
    # Free BMesh of reference object
    bm.free()


def select_loop_faces(loop):
    """
    Select all faces on the specified loop.
    Returns a list of the selected faces in order.
    Only works for quads!
    """
    faces = []
    
    while not loop.face.select:
        # If radial loop links back here, we're boundary, thus done
        if loop.link_loop_radial_next == loop:
            break
        
        # Remember and mark current face
        loop.face.select = True
        faces.append(loop.face)
        
        # Jump to adjacent face and walk two edges forward
        loop = loop.link_loop_radial_next.link_loop_next.link_loop_next

    return faces


def get_loop_direction(loop):
    """
    Returns the direction of a given edge loop relative to the face it's connected to.
    The direction is defined through the 2D angle of the v1 to v2, where v1 is the direction of the
    loop's current edge to the opposite one on the same face, and v2 is the normal direction of the loops face.
    Positive Angle -> 1.0, Negative Angle -> -1.0.
    Only works for quads!
    """
    face = loop.face
    edge1 = loop.edge
    edge2 = loop.link_loop_next.link_loop_next.edge
    
    direction = edge2.verts[0].co - edge1.verts[0].co
    angle = get_signed_angle(direction.to_2d(), face.normal.to_2d())
    
    return 1.0 if angle > 0 else -1.0

    
def get_horizontal_face_rings(bm):
    """
    Returns all face rings along walls (horizontal faces) for a given mesh object.
    (Requires edit mode, probably)
    """
    face_rings = []
    
    for edge in [e for e in bm.edges if abs((e.verts[0].co - e.verts[1].co).z) > 0.0001]:
        # Loop through all vertical edges
        if not edge.select:
            edge.select = True
            # Append all faces of edge's forward loop to list. Already processed edges / faces are selected.
            forward_loop = edge.link_loops[0] if get_loop_direction(edge.link_loops[0]) > 0 else edge.link_loops[1]
            loop_faces = select_loop_faces(forward_loop)
            face_rings.append(loop_faces)
    
    return face_rings


def uv_unwrap_walls(mesh_object, material_per_face_ring=False):
    """
    Adds UV information to all walls (horizontal faces) of the specified mesh object.
    Will group UVs together by continuous face loops (inner "rooms" or outer perimeter in case of walls).
    This way, the amount of UV seams is minimized (one per continuous edge loop).
    1.0 in UV space is mapped to 1.0 in 3D-Space.
    """
    me = mesh_object.data
    bm = bmesh.new()
    bm = bmesh.from_edit_mesh(me)
    
    uv_layer = bm.loops.layers.uv.active
    if not uv_layer: uv_layer = bm.loops.layers.uv.new()
    
    face_rings = get_horizontal_face_rings(bm)
    for i in range(len(face_rings)):
        face_ring = face_rings[i]
        
        # Determine which material to use for the faces of this face ring
        material_index = 0
        if material_per_face_ring:
            material_index = add_wall_material(i+1, me)
            
        width_offset = 0.0
        for face in face_ring:
            # Set material of face
            face.material_index = material_index
            # Determine length and height of face
            face_width = face.edges[0].calc_length()
            face_height = face.edges[1].calc_length()
            for i in range(0, len(face.loops)):
                # Calculate UV position for loop of each vertex on face
                if i == 0:
                    # Loop for bottom-left vertex
                    uv = (width_offset, 0)
                elif i == 1:
                    # Loop for bottom-right vertex
                    uv = (width_offset+face_width, 0)
                elif i == 2:
                    # Loop for top-right vertex
                    uv = (width_offset+face_width, face_height)
                elif i == 3:
                    # Loop for top-left vertex
                    uv = (width_offset, face_height)
                
                # Assign UV position to UV-map
                face.loops[i][uv_layer].uv = uv
            width_offset += face_width


def create_default_material(material_name="Default Material", hue=0.0, saturation=0.0):
    """
    Returns the default checkerboard material for walls.
    Hue and Saturation parameters can be used to specif a color tint.
    """
    # Create the new material
    material = bpy.data.materials.new(name=material_name)
    material.use_nodes = True
    
    # Create nodes
    nodes = material.node_tree.nodes
    node_texture_coordinates = nodes.new("ShaderNodeTexCoord")
    node_texture_coordinates.location = (-540, 0)
    node_checkter_texture_1 = nodes.new("ShaderNodeTexChecker")
    node_checkter_texture_1.location = (-360, 100)
    node_checkter_texture_1.inputs[1].default_value = colorsys.hsv_to_rgb(hue, saturation, 0.75) + (1.0,)
    node_checkter_texture_1.inputs[2].default_value = colorsys.hsv_to_rgb(hue, saturation, 0.5) + (1.0,)
    node_checkter_texture_1.inputs[3].default_value = 6.0
    node_checkter_texture_2 = nodes.new("ShaderNodeTexChecker")
    node_checkter_texture_2.location = (-360, -100)
    node_checkter_texture_2.inputs[1].default_value = (1, 1, 1, 1)
    node_checkter_texture_2.inputs[2].default_value = (0, 0, 0, 1)
    node_checkter_texture_2.inputs[3].default_value = 6.0
    node_bump = nodes.new("ShaderNodeBump")
    node_bump.location = (-180, -100)
    node_bump.inputs[0].default_value = 0.5
    node_principled = nodes.get("Principled BSDF")
    
    # Link nodes
    links = material.node_tree.links
    links.new(node_texture_coordinates.outputs[2], node_checkter_texture_1.inputs[0])
    links.new(node_texture_coordinates.outputs[2], node_checkter_texture_2.inputs[0])
    links.new(node_checkter_texture_2.outputs[1], node_bump.inputs[2])
    links.new(node_checkter_texture_1.outputs[0], node_principled.inputs[0])
    links.new(node_bump.outputs[0], node_principled.inputs[19])
    
    return material


def add_wall_material(index, mesh):
    """
    Adds the wall material with the specified index to the mesh and returns it's index.
    0 adds the base material with no tint, >= 1 adds indexed materials with hue tint.
    """
    if index == 0:
        material_name = "Wall Material"
        material_hue = 0.0
        material_saturation = 0.0
    else:
        material_name = f"Wall Material {index}"
        material_hue = index * (0.07) % 1
        material_saturation = 0.5
        
    wall_material = bpy.data.materials.get(material_name)
    if wall_material is None:
        wall_material = create_default_material(material_name, material_hue, material_saturation)
    
    # Add material to mesh if not already existing
    if (not mesh.materials) or (not material_name in mesh.materials):
        mesh.materials.append(wall_material)
    
    return mesh.materials.find(material_name)


def main(wall_thickness=0.125, wall_height=2.5, fill_rims=True, material_per_face_ring=False, context=bpy.context):
    """
    Creates a new object containing a mesh with wall geometry based on the currently selected object.
    The reference object should contain a mesh with a 2D wireframe representation of the wall layout, without faces.
    Adds the new wall object to the same collections the reference object is linked to.
    Hides the reference object after execution and selects the newly created wall object.
    """
    scene = context.scene
    reference_object = context.active_object
    
    # Create new mesh linked to new object in scene
    new_mesh = bpy.data.meshes.new(reference_object.name + " Walls")  # add the new mesh
    generate_wall_mesh_data(reference_object, new_mesh, wall_thickness, fill_rims) # add the mesh data
    new_obj = bpy.data.objects.new(new_mesh.name, new_mesh) # add a new object containing the mesh

    # Link new object to all collections the old object is a part of
    for collection in reference_object.users_collection:
        collection.objects.link(new_obj)

    # Hide and de-select reference object
    reference_object.select_set(False)
    reference_object.hide_set(True)

    # Select new object
    context.view_layer.objects.active = new_obj
    new_obj.select_set(True)
    
    # Add default material to new object's mesh
    add_wall_material(0, new_obj.data)

    # Edit new object
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.select_all(action='SELECT')
    # Remove doubles
    bpy.ops.mesh.remove_doubles(threshold=0.0001)
    # Extrude upwards
    bpy.ops.mesh.extrude_region_move(TRANSFORM_OT_translate={"value":(0, 0, wall_height)})
    # Unwrap UV-Information for horizontal face loops (walls)
    bpy.ops.mesh.select_all(action='DESELECT')
    uv_unwrap_walls(bpy.context.edit_object, material_per_face_ring)
    bpy.ops.mesh.select_all(action='SELECT')
    bpy.ops.object.mode_set(mode='OBJECT')


class OBJECT_OT_create_walls(bpy.types.Operator):
    """
    WallBuilder "Create Walls" Operator
    """
    bl_idname = "wallbuilder.create_walls"
    bl_label = "Create Walls"
    bl_description = "Creates a new mesh object with wall geometry for a given wireframe blueprint."
    bl_options = {'REGISTER', 'UNDO'}
    
    wall_thickness: FloatProperty(name="Wall Thickness", default=0.25)
    wall_height: FloatProperty(name="Wall Height", default=2.5)
    fill_rims: BoolProperty(name="Fill Rims", default=True)
    material_per_face_ring: BoolProperty(name="Material per face ring", default=False)
    
    @classmethod
    def poll(cls, context):
        return (context.mode == 'OBJECT' and context.active_object is not None)
    
    def execute(self, context):
        # Execute main function
        main(
            self.wall_thickness/2, 
            self.wall_height, 
            self.fill_rims, 
            self.material_per_face_ring, 
            context=context
        )
        return {'FINISHED'}
    
    def invoke(self, context, event):
        try:
            # Execute main function
            main(
                self.wall_thickness/2, 
                self.wall_height, 
                self.fill_rims, 
                self.material_per_face_ring, 
                context=context
            )
            return {'FINISHED'}
        except Exception:
            self.report({'WARNING'}, "Failed to construct walls for reference object. Make sure the reference object contains a mesh with a 2D wireframe outline of the wall layout.")
            return {'CANCELLED'}


class OBJECT_MT_wallbuilder(bpy.types.Menu):
    """
    WallBuilder Menu
    """
    bl_idname = "OBJECT_MT_wallbuilder_menu"
    bl_label = "WallBuilder"
    
    def draw(self, context):
        layout = self.layout
        layout.operator(OBJECT_OT_create_walls.bl_idname)


def menu_draw(self, context):
    self.layout.operator_context = 'INVOKE_REGION_WIN'
    self.layout.menu(OBJECT_MT_wallbuilder.bl_idname)


def register():
    register_class(OBJECT_MT_wallbuilder)
    register_class(OBJECT_OT_create_walls)
    bpy.types.VIEW3D_MT_object.append(menu_draw)
    

def unregister():
    bpy.types.VIEW3D_MT_object.remove(menu_draw)
    unregister_class(OBJECT_MT_wallbuilder)
    unregister_class(OBJECT_OT_create_walls)


if __name__ == "__main__":
    # Unregister the operator class if it is already registered
    try:
        unregister()
    except RuntimeError:
        pass
    
    # (Re-)register the operator class
    register()
