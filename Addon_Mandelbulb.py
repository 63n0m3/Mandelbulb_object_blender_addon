# This is the code that generates customizeable 3d mandelbulb based on the Wikipedia formula( look at the picture in the docs for furmula)
# This addon can create grid of squares representing mandelbulb including vertex color map.
# Alternatively it can create a surface mesh of mandelbulb without vertex color map
# Depending on the options selected calculation takes different amount of time
# Resolution is a parameter representing number of voxel units in x,y,z axies
# Control this number as for example 128 resolution gives for example around 4 500 000 vertices( in form of voxel cubes )
# max_r is the parameter that wikipedia doesnt present. It is a cutoff parameter for determining length of calculated voxel position vector calculated x^2+y^2+z^2 will be compared to this value and abandoned if lower
# Parameters are named according to wikipedia's naming scheme. wikipedia.org/wiki/Mandelbulb
# Licence: MIT free use, distribute and profit with this copyright notice
# by Gen0me https://github.com/63n0m3/Mandelbulb_object_blender_addon/
# btc: bc1q2h9m4pzfdnnfu79hfnjlcnta0j4etjqr4n5gy0
bl_info = {
    "name": "Mandelbulb 3D",
    "author": "Gen0me",
    "version": (1, 0),
    "blender": (2, 80, 0),
    "location": "View3D > Add > Mesh > Mandelbulb 3D",
    "description": "Create customizeable 3d fractal",
    "warning": "",
    "doc_url": "github.com/63n0m3",
    "category": "Add Mesh",
}

import bpy, math, numpy
from bpy.types import Operator
from bpy.props import *
from bpy_extras.object_utils import AddObjectHelper, object_data_add
from mathutils import Vector
import bmesh

def add_object(self, context):

    verts = []
    edges = []
    faces = []
    tmp_loop_col = []
  #  interval = 0.03
 #   MBorder = 6
    verts_count = 0
    interval = self.size/self.resolution
    min = -self.size/2 + ( interval/2 )
    half_vox_size = interval/2 # Algorithm works this way that it iterates through every voxel and checks if it is inside mandelbulb
    for x in numpy.arange( min, min+self.size , interval ): # Those 3 loops iterate through every x,y,z of the voxels
        for y in numpy.arange( min, min+self.size , interval ):
            for z in numpy.arange( min, min+self.size , interval ):
                xfin = x # Asignment of the starting vector
                yfin = y
                zfin = z
                self.cx = x
                self.cy = y
                self.cz = z
                for i in range( 1, self.iteration+1 ): # This loop reiterates mandelbulb math updating new x,y,z vertor values
                    Rsq = xfin**2+yfin**2+zfin**2
                    R = Rsq**0.5 # Length of a vector
                    phi = math.atan(yfin/xfin)
                    theta = math.atan(math.sqrt(xfin**2+yfin**2)/zfin)
                    Rn = R**self.n
                    xfin = math.sin(self.n*theta) * math.cos(self.n*phi)
                    yfin = math.sin(self.n*theta) * math.sin(self.n*phi)
                    zfin = math.cos(self.n*theta)
                    xfin *= Rn + self.cx # Addition of the starting vector
                    yfin *= Rn + self.cy
                    zfin *= Rn + self.cz
                    if xfin**2 + yfin**2 + zfin**2 > self.max_r: # This is condition to check if voxel is inside mandelbulb
                        break
                    if i == self.iteration: # This checks if i-loop has gone the necessary number of iterations
                        verts.append([ x-half_vox_size, y-half_vox_size, z-half_vox_size ]) # Following is the creation of each seperate texel vertex and face data
                        verts.append([ x+half_vox_size, y-half_vox_size, z-half_vox_size ])
                        verts.append([ x+half_vox_size, y+half_vox_size, z-half_vox_size ])
                        verts.append([ x-half_vox_size, y+half_vox_size, z-half_vox_size ])
                        verts.append([ x-half_vox_size, y-half_vox_size, z+half_vox_size ])
                        verts.append([ x+half_vox_size, y-half_vox_size, z+half_vox_size ])
                        verts.append([ x+half_vox_size, y+half_vox_size, z+half_vox_size ])
                        verts.append([ x-half_vox_size, y+half_vox_size, z+half_vox_size ])

                        faces.append([ verts_count*8+3, verts_count*8+2, verts_count*8+1, verts_count*8 ])
                        faces.append([ verts_count*8+4, verts_count*8+5, verts_count*8+6, verts_count*8+7 ])
                        faces.append([ verts_count*8, verts_count*8+1, verts_count*8+5, verts_count*8+4 ])
                        faces.append([ verts_count*8+1, verts_count*8+2, verts_count*8+6, verts_count*8+5 ])
                        faces.append([ verts_count*8+2, verts_count*8+3, verts_count*8+7, verts_count*8+6 ])
                        faces.append([ verts_count*8, verts_count*8+4, verts_count*8+7, verts_count*8+3 ])

                        tmp_loop_col.append (xfin) # Vertex colors data
                        tmp_loop_col.append (yfin)
                        tmp_loop_col.append (zfin)
                        tmp_loop_col.append (Rn)

                        verts_count += 1


    mbrot3d = "Mandelbulb 3d" # Creation of the mandelbulb object inside the collection
    mbrot3dm = bpy.data.meshes.new(mbrot3d)
    mbrot3do = bpy.data.objects.new(mbrot3d, mbrot3dm)
    col = bpy.context.view_layer.active_layer_collection.collection
    col.objects.link(mbrot3do)
    bpy.context.view_layer.objects.active = mbrot3do
    mbrot3dm.from_pydata(verts, edges, faces)
    mod_rem = mbrot3do.modifiers.new('Remesh','REMESH')
    mod_rem.voxel_size = interval
    mod_rem.show_viewport = False
    if (self.remesh):
        bpy.ops.object.modifier_apply(modifier="Remesh")
        mod_rem.show_viewport = True

    bm = bmesh.new() # Following is the addition of the vertex colors to object using bmesh loops
    bm.from_mesh(mbrot3dm)
    cl = bm.loops.layers.color.new("visualization")
    count = 0
    for f in bm.faces:
        for loop in f.loops:
            tem = int(loop.vert.index/8)*4
            if tem <= len(tmp_loop_col):
                loop[cl] = (tmp_loop_col[tem], tmp_loop_col[tem+1], tmp_loop_col[tem+2], tmp_loop_col[tem+3])
                count += 1
    bm.to_mesh(mbrot3dm)

class OBJECT_OT_add_object(Operator, AddObjectHelper):
    """Create a new Mesh Object"""
    bl_idname = "mesh.add_object"
    bl_label = "Add Mesh Object"
    bl_options = {'REGISTER', 'UNDO'}

    #GUI box propherties
    size : FloatProperty(
        name = "Calculation Size",
        description = "Cube domain of calculation size",
        default = 2.5
    )
    n : FloatProperty(
        name = "n",
        description = "order",
        default = 5.0
    )
    iteration : IntProperty(
        name = "iteration",
        description = "iteration",
        default = 2
    )
    resolution : FloatProperty(
        name = "Resolution",
        description = "voxel count in each axis",
        default = 32.0
    )
    remesh : BoolProperty(
        name = "Apply remesh",
        description = "apply remesh to lower the output number of vertices",
    )
    max_r : FloatProperty(
        name = "max_r",
        description = "r above which calculated vectors are abandoned",
        default = 8.0
    )

    def execute(self, context):

        add_object(self, context)

        return {'FINISHED'}


# Registration

def add_object_button(self, context):
    self.layout.operator(
        OBJECT_OT_add_object.bl_idname,
        text="Mandelbulb fractal",
        icon='MONKEY')


# This allows you to right click on a button and link to documentation
def add_object_manual_map():
    url_manual_prefix = "https://docs.blender.org/manual/en/latest/"
    url_manual_mapping = (
        ("bpy.ops.mesh.add_object", "scene_layout/object/types.html"),
    )
    return url_manual_prefix, url_manual_mapping


def register():
    bpy.utils.register_class(OBJECT_OT_add_object)
    bpy.utils.register_manual_map(add_object_manual_map)
    bpy.types.VIEW3D_MT_mesh_add.append(add_object_button)


def unregister():
    bpy.utils.unregister_class(OBJECT_OT_add_object)
    bpy.utils.unregister_manual_map(add_object_manual_map)
    bpy.types.VIEW3D_MT_mesh_add.remove(add_object_button)


if __name__ == "__main__":
    register()
