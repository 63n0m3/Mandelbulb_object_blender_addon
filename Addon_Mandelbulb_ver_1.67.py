# Version 1.0
# works fully in python as blender script

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


# Version 2.0
# can work in 3 different ways:
# in python as python blender script or calling precompiled dll from blender python script that uses c++ calculations or OpenCL calculations( depending on option choosen by user ).
# Dll c++ is single-threaded. Same as python.
# Dll OpenCL works on a single GPU which you can choose( look at standard output: Window -> Toggle System_Console ).
# Version 2.0 has also automated material creation in nodes tree with visualization of vertex colors. ver 1.0 needs manual node creation( same vertex color data is available in version 1.0 ).
# Version 2.0 has more advanced mesh cleaning algorithms than version 1.0.
# Delete disconnected bunch grouped is a very accurate algorithm of cleaning but its complexity is overwhelming.
# Thin skin algorithms are based on comparing voxels to all six surrounding voxels that are face connected like in picture: "voxels calc order.jpg"
# Thick skin algorithms are based on comparing voxels to all 26 surrounding voxels that are face or vertex connected.
# Same Mandelbulb fractals in Version 1.0 and 2.0.
# Dll gets loaded and unloaded every time the script runs unless error ocurs( than it wont get unloaded ). This means slight kernel compilation delay at the beginning of OpenCL run.
# Updated 7.11.2023

# I do plan version 3.0. Version 2.0 works by calculating all voxels of mandelbulb, and afterwards cleaning the unnecessary voxels. Cleanning algorithms are more computionally intensive than generation algorithms by at least N number of elements( div by constant )( in thin skin not bunch grouped easy case. bunch grouped gets much worse )
# For bigger Mandelbulbs cleaning is necessary as tested blender object creation times are huge compared to generation times. There may be a faster way to do this.
# I am leaving timing data in standard output for user information. Python script algorithms are more integrated so generation and cleaning algorithms are bunched together without data rewriting phase.


# This python script when executed in blender will register new object in Add->Mesh-> Mandelbulb fractal.


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
import bmesh
import ctypes
import os
import time

dll_path = os.path.dirname(os.path.abspath(bpy.context.space_data.text.filepath)) + '\\Mandelbulb_Gen0me.dll'

def calc_type_callback(self, context):
    return(
        ('PYTON', 'Internal - Python', " "),
        ('DLL_CPU', 'Dll - Cpu Single Thread', " "),
        ('DLL_OPENCL', 'Dll - OpenCl', " "),
    )
def m_type_callback(self, context):
    return(
        ('JUL', 'Juliabulb', "n, phase_theta, phase_phi available"),
        ('CUB', 'Cubic', " "),
        ('QUAD', 'Quadratic', " "),
        ('QUIN', 'Quintic', " "),
        ('POW9', 'Power 9', " "),
    )
def di_type_callback(self, context):
    return(
        ('NO_DI', 'Dont delete internal', " "),
        ('THIN_SKIN', 'Thin skin', "Faster"),
        ('THICK_SKIN', 'Thick skin', "Slower"),
    )
def dd_type_callback(self, context):
    return(
        ('NO_DD', 'Dont delete disconnected', " "),
        ('THIN_SKIN_CONNECTED_DIRECTLY', 'Thin skin - connected directly', " "),
        ('THICK_SKIN_CONNECTED_DIRECTLY', 'Thick skin - connected directly', " "),
        ('THIN_SKIN_CONNECTED_IN_BUNCH', 'Thin skin - bunch grouped', " "),
        ('THICK_SKIN_CONNECTED_IN_BUNCH', 'Thick skin - bunch grouped', " "),
    )


def add_object(self, context):

    start_time = time.time()*1000

    tmp_loop_col = []
    voxels = []
    verts = []
    edges = []
    faces = []

    interval = self.size/self.resolution
    half_vox_size = interval/2

    if self.calc_type == 'PYTON':   # for python calculations
        print ("Generating mandelbulb: ")
        verts_count = 0
        min = -self.size/2 + ( interval/2 ) # Algorithm works this way that it iterates through every voxel and checks if it is inside mandelbulb
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
                            if i == self.iteration: # This checks if i-loop has gone the necessary number of iterations. If it went through it means the vertex is valid
                                voxels.append([ x, y, z ]) # Vertex coordinates
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
                            if i == self.iteration:
                                voxels.append([ x, y, z ])
                                tmp_loop_col.append (xfin)
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
                            if i == self.iteration:
                                voxels.append([ x, y, z ])
                                tmp_loop_col.append (xfin)
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
                            if i == self.iteration:
                                voxels.append([ x, y, z ])
                                tmp_loop_col.append (xfin)
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
                            if i == self.iteration:
                                voxels.append([ x, y, z ])
                                tmp_loop_col.append (xfin)
                                tmp_loop_col.append (yfin)
                                tmp_loop_col.append (zfin)
                                tmp_loop_col.append (Rsq**(1/2))
                                verts_count += 1

            # so this is the algorithm that calculates Delete disconnected. It is so slow because it has to iterate every voxel by every voxel multiple times.
        if self.dd_type == 'THIN_SKIN_CONNECTED_DIRECTLY' or self.di_type == 'THIN_SKIN':
            print ("And deleting internal/disconnected: ")
            indices_connected = [[ -1 for x in range(6)] for y in range(verts_count)]
            for i in range(0, verts_count):     # This iterates every voxel by every voxel to check how many voxels are in each voxels surroundings
                indice_count = 0
                for j in range(0, verts_count):
                    if (voxels[i][0] == voxels[j][0]):
                        if (voxels[i][1] == voxels[j][1]):
                            if (voxels[i][2] < voxels[j][2] - interval*0.9 and voxels[i][2] > voxels[j][2] - interval*1.1 ) or (voxels[i][2] > voxels[j][2] + interval*0.9 and voxels[i][2] < voxels[j][2] + interval*1.1):
                                indices_connected[i][indice_count] = j
                                indice_count += 1
                    if (voxels[i][0] == voxels[j][0]):
                        if (voxels[i][2] == voxels[j][2]):
                            if (voxels[i][1] < voxels[j][1] - interval*0.9 and voxels[i][1] > voxels[j][1] - interval*1.1 ) or (voxels[i][1] > voxels[j][1] + interval*0.9 and voxels[i][1] < voxels[j][1] + interval*1.1):
                                indices_connected[i][indice_count] = j
                                indice_count += 1
                    if (voxels[i][1] == voxels[j][1]):
                        if (voxels[i][2] == voxels[j][2]):
                            if (voxels[i][0] < voxels[j][0] - interval*0.9 and voxels[i][0] > voxels[j][0] - interval*1.1 ) or (voxels[i][0] > voxels[j][0] + interval*0.9 and voxels[i][0] < voxels[j][0] + interval*1.1):
                                indices_connected[i][indice_count] = j
                                indice_count += 1
            # Now inside indices_connected[n][0-5] should be indices of connected verts

            if self.dd_type == 'THIN_SKIN_CONNECTED_DIRECTLY':     # So here comes the fast algorithm. Its criterion is the number of surrounding voxels
                delete_list_indices = []
                for q in range(0, verts_count):
                    connected_num = 0
                    for r in range(0, 6):
                        if indices_connected[q][r] != -1:
                            connected_num += 1
                    if connected_num <= self.del_disc_crit:
                        delete_list_indices.append(q)
                delete_list_indices = sorted(delete_list_indices, reverse = True)
                for u in range(0, len(delete_list_indices)):
                    to_del = delete_list_indices[u]
                    voxels.pop(to_del)
                    tmp_loop_col.pop(4*to_del+3)
                    tmp_loop_col.pop(4*to_del+2)
                    tmp_loop_col.pop(4*to_del+1)
                    tmp_loop_col.pop(4*to_del)
                verts_count -= len(delete_list_indices)


            # so this is the algorithm that calculates Delete internal. It is so slow because it has to iterate every voxel by every voxel up to 3 times until hits are met.
            if self.di_type == 'THIN_SKIN':
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

        elif self.di_type == 'THICK_SKIN' or self.dd_type == 'THICK_SKIN_CONNECTED_DIRECTLY' or self.dd_type == 'THIN_SKIN_CONNECTED_IN_BUNCH' or self.dd_type == 'THICK_SKIN_CONNECTED_IN_BUNCH':
            print ("Bunch grouped and thick algorithm python versions do not work, there are CPU and OpenCL versions of those")

        if self.cub_or_vert == True:    # Whether to generate cubes or just leave voxel vertices
            for n in range(0, verts_count):
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

        print (round (time.time()*1000 - start_time, 2))

    else:   # if using dll
        compute_dll = ctypes.WinDLL( dll_path )
        gpu_compute_proto = ctypes.WINFUNCTYPE (
            ctypes.c_int,    #return
            ctypes.c_int,    #param 1
            ctypes.c_int,
            ctypes.c_int,
            ctypes.c_float,
            ctypes.c_int,
            ctypes.c_float,
            ctypes.c_float,
            ctypes.c_float,
            ctypes.c_float,
            ctypes.c_float)
        gpu_compute_params = ( 1, "p1", 0 ),( 1, "p2", 0 ),( 1, "p3", 0 ),( 1, "p4", 0 ),( 1, "p5", 0 ),( 1, "p6", 0 ),( 1, "p7", 0 ),( 1, "p8", 0 ),( 1, "p9", 0 ),( 1, "p10", 0 )
        gpu_compute_api = gpu_compute_proto(("_Z44Calculate_arrays_and_get_number_of_verts_GPUiiififffff", compute_dll), gpu_compute_params)

        cpu_compute_proto = ctypes.WINFUNCTYPE (
            ctypes.c_int,
            ctypes.c_int,
            ctypes.c_int,
            ctypes.c_float,
            ctypes.c_int,
            ctypes.c_float,
            ctypes.c_float,
            ctypes.c_float,
            ctypes.c_float,
            ctypes.c_float)
        cpu_compute_params = ( 1, "p1", 0 ),( 1, "p2", 0 ),( 1, "p3", 0 ),( 1, "p4", 0 ),( 1, "p5", 0 ),( 1, "p6", 0 ),( 1, "p7", 0 ),( 1, "p8", 0 ),( 1, "p9", 0 )
        cpu_compute_api = cpu_compute_proto(("_Z44Calculate_arrays_and_get_number_of_verts_CPUiififffff", compute_dll), cpu_compute_params)

        del_internal_disconnected_proto = ctypes.WINFUNCTYPE (
            ctypes.c_int,
            ctypes.c_int,
            ctypes.c_int,
            ctypes.c_int)
        del_internal_disconnected_params = ( 1, "p1", 0 ),( 1, "p2", 0 ),( 1, "p3", 0 )
        cpu_del_internal_disconnected_api = del_internal_disconnected_proto(("_Z32Delete_internal_disconnected_CPUiii", compute_dll), del_internal_disconnected_params)
        gpu_del_internal_disconnected_api = del_internal_disconnected_proto(("_Z32Delete_internal_disconnected_GPUiii", compute_dll), del_internal_disconnected_params)

        copy_data_proto = ctypes.WINFUNCTYPE (
            ctypes.c_int,
            ctypes.POINTER(ctypes.c_float) )
        copy_data_params = ( 1, "p1", 0 ),
        copy_data_gpu_api = copy_data_proto(("_Z18Export_voxel_d_GPUPf", compute_dll), copy_data_params)
        copy_data_cpu_api = copy_data_proto(("_Z18Export_voxel_d_CPUPf", compute_dll), copy_data_params)


        if self.m_type == 'JUL':     # no enum in ctypes
            m_typ = 0
        if self.m_type == 'CUB':
            m_typ = 1
        if self.m_type == 'QUAD':
            m_typ = 2
        if self.m_type == 'QUIN':
            m_typ = 3
        if self.m_type == 'POW9':
            m_typ = 4

        if self.di_type == 'NO_DI':
            di_typ = 0
        if self.di_type == 'THIN_SKIN':
            di_typ = 1
        if self.di_type == 'THICK_SKIN':
            di_typ = 2

        if self.dd_type == 'NO_DD':
            dd_typ = 0
        if self.dd_type == 'THIN_SKIN_CONNECTED_DIRECTLY':
            dd_typ = 1
        if self.dd_type == 'THICK_SKIN_CONNECTED_DIRECTLY':
            dd_typ = 2
        if self.dd_type == 'THIN_SKIN_CONNECTED_IN_BUNCH':
            dd_typ = 3
        if self.dd_type == 'THICK_SKIN_CONNECTED_IN_BUNCH':
            dd_typ = 4


        if self.calc_type == 'DLL_OPENCL':
            start_time = time.time()*1000
            print ("Generating mandelbulb: ")
            verts_count = gpu_compute_api( ctypes.c_int(self.ocl_device), ctypes.c_int(m_typ), ctypes.c_int(self.iteration), ctypes.c_float(self.size), ctypes.c_int(self.resolution), ctypes.c_float(self.n), ctypes.c_float(self.phase_theta), ctypes.c_float(self.phase_phi), ctypes.c_float(self.max_r), ctypes.c_float(self.param_D) )

            print (round (time.time()*1000 - start_time, 2))
            start_time = time.time()*1000

            if self.di_type != 'NO_DI' or self.dd_type != 'NO_DD':
                print ("Deleting internal/disconnected: ")
                verts_count = gpu_del_internal_disconnected_api( ctypes.c_int(di_typ), ctypes.c_int(dd_typ), ctypes.c_int(self.del_disc_crit))
                print (round (time.time()*1000 - start_time, 2))
                start_time = time.time()*1000

            print ("Writing data to python: ")

            if self.cub_or_vert == True:
                voxel_d = (ctypes.c_float*(7*verts_count))()
                copy_data_gpu_api(ctypes.cast(voxel_d, ctypes.POINTER(ctypes.c_float)))

                for n in range(0, verts_count):
                    x = voxel_d[7*n]
                    y = voxel_d[7*n+1]
                    z = voxel_d[7*n+2]

                    verts.append([ x-half_vox_size, y-half_vox_size, z-half_vox_size ])
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

                    tmp_loop_col.append(voxel_d[7*n+3])
                    tmp_loop_col.append(voxel_d[7*n+4])
                    tmp_loop_col.append(voxel_d[7*n+5])
                    tmp_loop_col.append(voxel_d[7*n+6])

            else:   # self.cub_or_vert == False
                voxel_d = (ctypes.c_float*(7*verts_count))()
                copy_data_gpu_api(ctypes.cast(voxel_d, ctypes.POINTER(ctypes.c_float)))

                for i in range(0, verts_count):                                       # move data from outside dll to known python-blender structures
                    voxels.append([ voxel_d[7*i], voxel_d[7*i+1], voxel_d[7*i+2] ])
                    tmp_loop_col.append (voxel_d[7*i+3])
                    tmp_loop_col.append (voxel_d[7*i+4])
                    tmp_loop_col.append (voxel_d[7*i+5])
                    tmp_loop_col.append (voxel_d[7*i+6])

            print (round (time.time()*1000 - start_time, 2))

        if self.calc_type == 'DLL_CPU':
            start_time = time.time()*1000
            print ("Generating mandelbulb: ")
            verts_count = cpu_compute_api(ctypes.c_int(m_typ), ctypes.c_int(self.iteration), ctypes.c_float(self.size), ctypes.c_int(self.resolution), ctypes.c_float(self.n), ctypes.c_float(self.phase_theta), ctypes.c_float(self.phase_phi), ctypes.c_float(self.max_r), ctypes.c_float(self.param_D))
            print (round (time.time()*1000 - start_time, 2))
            start_time = time.time()*1000

            if self.di_type != 'NO_DI' or self.dd_type != 'NO_DD':
                print ("Deleting internal/disconnected: ")
                verts_count = cpu_del_internal_disconnected_api(ctypes.c_int(di_typ), ctypes.c_int(dd_typ), ctypes.c_int(self.del_disc_crit))
                print (round (time.time()*1000 - start_time, 2))
                start_time = time.time()*1000
            print ("Writing data to python: ")
            voxel_d = (ctypes.c_float*(7*verts_count))()
            copy_data_cpu_api(ctypes.cast(voxel_d, ctypes.POINTER(ctypes.c_float)))

            for i in range(0, verts_count):
                tmp_loop_col.append (voxel_d[7*i+3])
                tmp_loop_col.append (voxel_d[7*i+4])
                tmp_loop_col.append (voxel_d[7*i+5])
                tmp_loop_col.append (voxel_d[7*i+6])

            if self.cub_or_vert == True:    # Whether to generate cubes or just leave voxel vertices
                for n in range(0, verts_count):
                    x = voxel_d[7*n]
                    y = voxel_d[7*n+1]
                    z = voxel_d[7*n+2]

                    verts.append([ x-half_vox_size, y-half_vox_size, z-half_vox_size ])
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

            else: # self.cub_or_vert == False
                for i in range(0, verts_count):
                    voxels.append([ voxel_d[7*i], voxel_d[7*i+1], voxel_d[7*i+2] ])

            print (round (time.time()*1000 - start_time, 2))

        libHandle = compute_dll._handle
        del compute_dll
        ctypes.windll.kernel32.FreeLibrary(libHandle)

    start_time = time.time()*1000
    print ("Creating object: ")

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
    mat = bpy.data.materials.new(name="Vertex_Colors_Mat")
    mat.use_nodes = True
    mbrot3do.data.materials.append(mat)
    mat_nodes = mat.node_tree.nodes
    mat_links = mat.node_tree.links
    snvc = mat_nodes.new("ShaderNodeVertexColor")
    snvc.location = (-250,220)
    snvc.layer_name = 'visualization'
    mat_links.new(snvc.outputs['Color'], mat_nodes['Principled BSDF'].inputs['Base Color'])

    print (round (time.time()*1000 - start_time, 2))
    print ("Done.")


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
    calc_type : EnumProperty(
        name = "Computation type:",
        items = calc_type_callback,
    )
    ocl_device : IntProperty(
        name = "OpenCL device:",
        description = "Check 'Window -> Toggle System Console'",
        default = 0
    )
    resolution : IntProperty(
        name = "Resolution",
        description = "voxel count in each axis",
        default = 16
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
    max_r : FloatProperty(
        name = "max_r",
        description = "r above which calculated vectors are abandoned",
        default = 8.0
    )
    di_type : EnumProperty(
        name = "Delete internal type:",
        items = di_type_callback,
    )
    dd_type : EnumProperty(
        name = "Delete disconnected type:",
        items = dd_type_callback,
    )
    del_disc_crit : IntProperty(
        name = "Delete disconnected size",
        description = "Criterion for deleting disconnected voxels. In bunch group it is the max number of voxels in a connected bunch. In connected directly it is the max number of voxel's neighbouring voxels",
        default = 1
    )
    cub_or_vert : BoolProperty(
        name = "Generate cubes/vertices",
        description = "tick to generate cubes, untick for vertices",
        default = True
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

