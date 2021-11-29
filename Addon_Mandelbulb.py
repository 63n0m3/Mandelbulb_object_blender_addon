# This is the code that generates customizeable 3d mandelbulb based on the Wikipedia formula( look at the picture in the docs for furmula)
# This addon can create grid of squares representing mandelbulb including vertex color map.
# Alternatively it can create a surface mesh of mandelbulb without vertex color map
# Depending on the options selected calculation takes different amount of time
# Resolution is a parameter representing number of voxel units in x,y,z axies
# Control this number as for example 128 resolution gives for example around 4 500 000 vertices( in form of voxel cubes ). You can alternatively generate just vertices but you will lose color maps. And still you will be left with problem of generating mesh from them...
# In case you generate voxels their center is on the center of algorithms criterion, in case of vertices, vertices are in the centre of algorithm's criterion.
# max_r is the parameter that wikipedia doesnt present. It is a cutoff parameter for determining length of calculated voxel position vector calculated x^2+y^2+z^2 will be compared to this value and abandoned if lower
# Generated vertex colors(including alpha) are colored according to difference in distance of primary x,y,z and r=sqrt(x**2+y**2+z**2) to calculated by algorithm iteration resaults.
# Parameters are named according to wikipedia's naming scheme. wikipedia.org/wiki/Mandelbulb
# phase_theta and phase_phi give this skin flowing effect
# n/A, phase_theta/B, phase_phi/C, D parameters available only on Juliabulb and Quintic algorithms
# Delete internal is so slow, i do not recomment using it on bigger resolutions.
# Licence: MIT free use, distribute and profit with this copyright notice
# by Gen0me https://github.com/63n0m3/Mandelbulb_object_blender_addon/
# Updated 29.11.2021
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

def m_type_callback(self, context):
    return(
        ('JUL', 'Juliabulb', "n, phase_theta, phase_phi available"),
        ('CUB', 'Cubic', " "),
        ('QUAD', 'Quadratic', " "),
        ('QUIN', 'Quintic', "A, B, C, D available"),
        ('POW9', 'Power 9', " "),

    )

def add_object(self, context):

    verts = []
    edges = []
    faces = []
    tmp_loop_col = []
    voxels = []
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
                    if self.m_type == 'JUL':           # Juliabulb
                        Rsq = xfin**2+yfin**2+zfin**2
                        R = Rsq**0.5 # Length of a vector
                        phi = math.atan(yfin/xfin) + self.phase_phi
                        theta = math.atan(math.sqrt(xfin**2+yfin**2)/zfin) + self.phase_theta
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
                            voxels.append([ x, y, z ])
                            tmp_loop_col.append (xfin) # Vertex colors data
                            tmp_loop_col.append (yfin)
                            tmp_loop_col.append (zfin)
                            tmp_loop_col.append (Rn)
                            verts_count += 1
                    elif self.m_type == 'POW9':        # Power Nine
                        xfin = xfin**9 - 36*(xfin**7)*(yfin**2 + zfin**2) + 126*(xfin**5)*((yfin**2 + zfin**2)**2) - 84*(xfin**3)*((yfin**2 + zfin**2)**3) + 9*xfin*((yfin**2 + zfin**2)**4) + self.cx
                        yfin = yfin**9 - 36*(yfin**7)*(zfin**2 + xfin**2) + 126*(yfin**5)*((zfin**2 + xfin**2)**2) - 84*(yfin**3)*((zfin**2 + xfin**2)**3) + 9*yfin*((zfin**2 + xfin**2)**4) + self.cy
                        zfin = zfin**9 - 36*(zfin**7)*(xfin**2 + yfin**2) + 126*(zfin**5)*((xfin**2 + yfin**2)**2) - 84*(zfin**3)*((xfin**2 + yfin**2)**3) + 9*zfin*((xfin**2 + yfin**2)**4) + self.cz
                        Rsq = xfin**2+yfin**2+zfin**2
                        if Rsq > self.max_r:
                            break;
                        if i == self.iteration: # This checks if i-loop has gone the necessary number of iterations
                            voxels.append([ x, y, z ])
                            tmp_loop_col.append (xfin) # Vertex colors data
                            tmp_loop_col.append (yfin)
                            tmp_loop_col.append (zfin)
                            tmp_loop_col.append (Rsq**(1/2))
                            verts_count += 1
                    elif self.m_type == 'QUIN':      # Quintic
                        xfin = xfin**5 - 10*(xfin**3)*(yfin**2 + self.n*yfin*zfin + zfin**2) + 5*xfin*(yfin**4 + self.phase_theta*(yfin**3)*zfin + self.phase_phi*(yfin**2)*(zfin**2) + self.phase_theta*yfin*(zfin**3) + zfin**4) + self.param_D*(xfin**2)*yfin*zfin*(yfin+zfin) + self.cx
                        yfin = yfin**5 - 10*(yfin**3)*(zfin**2 + self.n*xfin*zfin + xfin**2) + 5*yfin*(zfin**4 + self.phase_theta*(zfin**3)*xfin + self.phase_phi*(zfin**2)*(xfin**2) + self.phase_theta*zfin*(xfin**3) + xfin**4) + self.param_D*(yfin**2)*zfin*xfin*(zfin+xfin) + self.cy
                        zfin = zfin**5 - 10*(zfin**3)*(xfin**2 + self.n*xfin*yfin + yfin**2) + 5*zfin*(xfin**4 + self.phase_theta*(xfin**3)*yfin + self.phase_phi*(xfin**2)*(yfin**2) + self.phase_theta*xfin*(yfin**3) + yfin**4) + self.param_D*(zfin**2)*xfin*yfin*(xfin+yfin) + self.cz
                        Rsq = xfin**2+yfin**2+zfin**2
                        if Rsq > self.max_r:
                            break;
                        if i == self.iteration: # This checks if i-loop has gone the necessary number of iterations
                            voxels.append([ x, y, z ])
                            tmp_loop_col.append (xfin) # Vertex colors data
                            tmp_loop_col.append (yfin)
                            tmp_loop_col.append (zfin)
                            tmp_loop_col.append (Rsq**(1/2))
                            verts_count += 1
                    elif self.m_type == 'CUB':      # Cubic
                        xfin = xfin**3 - 3*xfin*(yfin**2 + zfin**2) + self.cx
                        yfin = - yfin**3 + 3*yfin*(xfin**2) - yfin*(zfin**2) + self.cy
                        zfin = zfin**3 - 3*zfin*(xfin**2) + zfin*(yfin**2) + self.cz
                        Rsq = xfin**2+yfin**2+zfin**2
                        if Rsq > self.max_r:
                            break;
                        if i == self.iteration: # This checks if i-loop has gone the necessary number of iterations
                            voxels.append([ x, y, z ])
                            tmp_loop_col.append (xfin) # Vertex colors data
                            tmp_loop_col.append (yfin)
                            tmp_loop_col.append (zfin)
                            tmp_loop_col.append (Rsq**(1/2))
                            verts_count += 1
                    else: #eliif self.m_type == 'QUAD':      # Quadratic
                        xfin = xfin**2 - yfin**2 + self.cx
                        yfin = 2*xfin*zfin + self.cy
                        zfin = 2*xfin*yfin + self.cz
                        Rsq = xfin**2+yfin**2+zfin**2
                        if Rsq > self.max_r:
                            break;
                        if i == self.iteration: # This checks if i-loop has gone the necessary number of iterations
                            voxels.append([ x, y, z ])
                            tmp_loop_col.append (xfin) # Vertex colors data
                            tmp_loop_col.append (yfin)
                            tmp_loop_col.append (zfin)
                            tmp_loop_col.append (Rsq**(1/2))
                            verts_count += 1


    # so this is the algorithm that calculates Delete disconnected. it is so slow because it has to iterate every voxel by every voxel  multiple times.
    if self.del_disconnected == True:
        indices_connected = [[ -1 for x in range(6)] for y in range(verts_count)]
        #for v in range(0, verts_count):
            #for c in range(0, 6):
                #print(indices_connected[v][c])    # Starting false value for connected voxel indices
        for i in range(0, verts_count):     # This iterates every voxel by every voxel to check how many voxels are in each voxels surroundings
            indice_a_count = 0
            indice_b_count = 2
            indice_c_count = 4
            for j in range(0, verts_count):
                if (voxels[i][0] == voxels[j][0]):
                    if (voxels[i][1] == voxels[j][1]):
                        if (voxels[i][2] < voxels[j][2] - interval*0.9 and voxels[i][2] > voxels[j][2] - interval*1.1 ) or (voxels[i][2] > voxels[j][2] + interval*0.9 and voxels[i][2] < voxels[j][2] + interval*1.1):
                            indices_connected[i][indice_a_count] = j
                            indice_a_count += 1
                if (voxels[i][0] == voxels[j][0]):
                    if (voxels[i][2] == voxels[j][2]):
                        if (voxels[i][1] < voxels[j][1] - interval*0.9 and voxels[i][1] > voxels[j][1] - interval*1.1 ) or (voxels[i][1] > voxels[j][1] + interval*0.9 and voxels[i][1] < voxels[j][1] + interval*1.1):
                            indices_connected[i][indice_b_count] = j
                            indice_b_count += 1
                if (voxels[i][1] == voxels[j][1]):
                    if (voxels[i][2] == voxels[j][2]):
                        if (voxels[i][0] < voxels[j][0] - interval*0.9 and voxels[i][0] > voxels[j][0] - interval*1.1 ) or (voxels[i][0] > voxels[j][0] + interval*0.9 and voxels[i][0] < voxels[j][0] + interval*1.1):
                            indices_connected[i][indice_b_count] = j
                            indice_c_count += 1
        # Now inside indices_connected[n][0-5] should be indices of connected verts
        
        # So here comes the fast algorithm. Its criterion is the number of surrounding voxels
        delete_list_indices = []
        for q in range( 0, verts_count ):
            connected_num = 0
            for r in range( 0, 6 ):
                if indices_connected[q][r] != -1 :
                    connected_num += 1
            if connected_num <= self.del_disc_crit :
                delete_list_indices.append(q)
        delete_list_indices = sorted ( delete_list_indices, reverse = True )
        print (delete_list_indices)
        for u in range( 0, len( delete_list_indices )) :
            to_del = delete_list_indices[ u ]
            voxels.pop(to_del)
            tmp_loop_col.pop(4*to_del+3)
            tmp_loop_col.pop(4*to_del+2)
            tmp_loop_col.pop(4*to_del+1)
            tmp_loop_col.pop(4*to_del)
        verts_count -= len( delete_list_indices )        
                    
    # so this is the algorithm that calculates Delete internal. it is so slow because it has to iterate every voxel by every voxel  up to 3 times until hits are met.
    if self.del_internal == True:
        indices_internal = []       # this is a list for indices of internal voxels
        for i in range(0, verts_count):     # This loop iterates through every generated voxel and checks if it is internal voxel and deletes it if option is selected
            count_near_voxels_a = 0
            count_near_voxels_b = 0
            count_near_voxels_c = 0
            for j in range(0, verts_count):
                if (voxels[i][0] == voxels[j][0]):
                    if (voxels[i][1] == voxels[j][1]):
                        if (voxels[i][2] < voxels[j][2] - interval*0.9 and voxels[i][2] > voxels[j][2] - interval*1.1 ) or (voxels[i][2] > voxels[j][2] + interval*0.9 and voxels[i][2] < voxels[j][2] + interval*1.1):            
                            count_near_voxels_a += 1
                if (voxels[i][0] == voxels[j][0]):
                    if (voxels[i][2] == voxels[j][2]):
                        if (voxels[i][1] < voxels[j][1] - interval*0.9 and voxels[i][1] > voxels[j][1] - interval*1.1 ) or (voxels[i][1] > voxels[j][1] + interval*0.9 and voxels[i][1] < voxels[j][1] + interval*1.1):
                            count_near_voxels_b += 1
                if (voxels[i][1] == voxels[j][1]):
                    if (voxels[i][2] == voxels[j][2]):
                        if (voxels[i][0] < voxels[j][0] - interval*0.9 and voxels[i][0] > voxels[j][0] - interval*1.1 ) or (voxels[i][0] > voxels[j][0] + interval*0.9 and voxels[i][0] < voxels[j][0] + interval*1.1):
                            count_near_voxels_c += 1
            if count_near_voxels_a == 2 and count_near_voxels_b == 2 and count_near_voxels_c == 2:
                indices_internal.append(i)
        for i in range(len(indices_internal)-1, -1, -1):    # voxels with calculated indices are substracted
            to_del = indices_internal[i]
            voxels.pop(to_del)
            tmp_loop_col.pop(4*to_del+3)
            tmp_loop_col.pop(4*to_del+2)
            tmp_loop_col.pop(4*to_del+1)
            tmp_loop_col.pop(4*to_del)
        verts_count -= len(indices_internal)
        
    if self.cub_or_vert == True:    # Whether to generate cubes or just left voxel vertices
        for n in range( 0, verts_count ):
            x = voxels[n][0]
            y = voxels[n][1]
            z = voxels[n][2]

            verts.append([ x-half_vox_size, y-half_vox_size, z-half_vox_size ]) # Following is the creation of each seperate texel vertex and face data
            verts.append([ x+half_vox_size, y-half_vox_size, z-half_vox_size ])
            verts.append([ x+half_vox_size, y+half_vox_size, z-half_vox_size ])
            verts.append([ x-half_vox_size, y+half_vox_size, z-half_vox_size ])
            verts.append([ x-half_vox_size, y-half_vox_size, z+half_vox_size ])
            verts.append([ x+half_vox_size, y-half_vox_size, z+half_vox_size ])
            verts.append([ x+half_vox_size, y+half_vox_size, z+half_vox_size ])
            verts.append([ x-half_vox_size, y+half_vox_size, z+half_vox_size ])
            
            faces.append([ n*8+3, n*8+2, n*8+1, n*8 ])
            faces.append([ n*8+4, n*8+5, n*8+6, n*8+7 ])
            faces.append([ n*8, n*8+1, n*8+5, n*8+4 ])
            faces.append([ n*8+1, n*8+2, n*8+6, n*8+5 ])
            faces.append([ n*8+2, n*8+3, n*8+7, n*8+6 ])
            faces.append([ n*8, n*8+4, n*8+7, n*8+3 ])
            
    mbrot3d = "Mandelbulb 3d" # Creation of the mandelbulb object inside the collection
    mbrot3dm = bpy.data.meshes.new(mbrot3d)
    mbrot3do = bpy.data.objects.new(mbrot3d, mbrot3dm)
    col = bpy.context.view_layer.active_layer_collection.collection
    col.objects.link(mbrot3do)
    bpy.context.view_layer.objects.active = mbrot3do
    
    if self.cub_or_vert == True:
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
    else:
        mbrot3dm.from_pydata(voxels, edges, faces)

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
    m_type : EnumProperty(
        name = "Mandelbulb type:",
        items = m_type_callback,
    )
    n : FloatProperty(
        name = "n/A",
        description = "order",
        default = 5.0
    )
    phase_theta : FloatProperty(
        name = "Phase theta/B",
        description = "skin flow effect",
        default = 0.0
    )
    phase_phi : FloatProperty(
        name = "Phase phi/C",
        description = "skin flow effect",
        default = 0.0
    )
    param_D : FloatProperty(
        name = "D",
        description = "D in quintic formula",
        default = 0.0
    )
    iteration : IntProperty(
        name = "iteration",
        description = "iteration",
        default = 2
    )
    resolution : FloatProperty(
        name = "Resolution",
        description = "voxel count in each axis",
        default = 16.0
    )
    max_r : FloatProperty(
        name = "max_r",
        description = "r above which calculated vectors are abandoned",
        default = 8.0
    )
    cub_or_vert : BoolProperty(
        name = "Generate cubes/vertices",
        description = "tick to generate cubes, untick for vertices",
        default = True
    )
    del_internal : BoolProperty(
        name = "Delete internal",
        description = "deletes internal voxels",
        default = False
    )
    del_disconnected : BoolProperty(
        name = "Delete disconnected",
        description = "delete disconnected voxels, it is a fast, not accurate algorithm",
        default = False
    )
    del_disc_crit : IntProperty(
        name = "Delete disconnected size",
        description = "Criterion for deleting disconnected voxels. It is the max number of voxel's neighbouring voxels",
        default = 1
    )
    remesh : BoolProperty(
        name = "Apply remesh",
        description = "apply remesh to lower the output number of vertices",
        default = False
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

