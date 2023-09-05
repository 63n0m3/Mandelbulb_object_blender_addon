///Add compiler option -c -DBUILD_DLL
///Add linker option -shared -o mandelbulb_opencl_63n0m3.dll
#include <windows.h>
#include <math.h>
#include <CL/opencl.hpp>
#include <vector>
#include <string>
#include <stdlib.h>             /// alloca
#include <iostream>

#define JUL 0
#define POW9 1
#define QUIN 2
#define CUB 3
#define QUAD 4

using namespace std;

vector <float> * results;
vector<cl::Platform> * platforms;
vector<cl::Device> * devices;
cl::Context * Cl_Context;
cl::Program * Prog;
cl::CommandQueue * Cl_Queue;
cl::Buffer * voxel_d;

/// returns valid voxel count and sets up results (voxel data array: x,y,z-coordinates,r,g,b,a-colours )
int Calculate_arrays_and_get_number_of_verts_CPU (int m_type, int iterations, float Size, int resolution, float n_order,
                                                  float phase_theta, float phase_phi, float max_r, float param_D){

    unsigned int voxel_d_size = pow(resolution,3)* 8;
    results = new vector <float>;
    results->resize( voxel_d_size, 0.0 );
    float interval = Size/(float)resolution;
    float min = -Size/2 + ( interval/2 );
    if (m_type==JUL){                     /// Juliabulb
        for (int x_i = 0; x_i<resolution; x_i++){
            float x = min + (interval*(float)x_i);
            for (int y_i = 0; y_i<resolution; y_i++){
                float y = min + (interval*(float)y_i);
                for (int z_i = 0; z_i<resolution; z_i++){                 /// Algorithm works this way that it iterates through every voxel and checks if it is inside mandelbulb
                    float z = min + (interval*(float)z_i);                  /// Those 3 loops / workgroups iterate through every x,y,z of the voxels
                    float xfin = x;                 /// Asignment of the starting vector
                    float yfin = y;
                    float zfin = z;
                    float cx = x;
                    float cy = y;
                    float cz = z;
                    for (int i = 1; i <= iterations; i++){                 /// This loop reiterates mandelbulb math updating new x,y,z vertor values
                        float Rsq = (xfin*xfin)+(yfin*yfin)+(zfin*zfin);
                        float R = sqrt(Rsq); /// Length of a vector
                        float phi = atan(yfin/xfin) + phase_phi;
                        float theta = atan(sqrt( (xfin*xfin)+(yfin*yfin))/zfin) + phase_theta;
                        float Rn = pow (R,n_order);
                        xfin = sin(n_order*theta) * cos(n_order*phi);
                        yfin = sin(n_order*theta) * sin(n_order*phi);
                        zfin = cos(n_order*theta);
                        xfin *= Rn + cx;        /// Addition of the starting vector
                        yfin *= Rn + cy;
                        zfin *= Rn + cz;
                        if (Rsq > max_r)                                    /// This is condition to check if voxel is inside mandelbulb
                            break;
                        if (i == iterations){                               /// This checks if i-loop has gone the necessary number of iterations. If it went through it means the vertex is valid
                            results->at( ( x_i + y_i*resolution + z_i*resolution*resolution  ) * 8 ) = x;           ///position
                            results->at( ( x_i + y_i*resolution + z_i*resolution*resolution  ) * 8 + 1 ) = y;
                            results->at( ( x_i + y_i*resolution + z_i*resolution*resolution  ) * 8 + 2 ) = z;
                            results->at( ( x_i + y_i*resolution + z_i*resolution*resolution  ) * 8 + 3 ) = xfin;    ///colors
                            results->at( ( x_i + y_i*resolution + z_i*resolution*resolution  ) * 8 + 4 ) = yfin;
                            results->at( ( x_i + y_i*resolution + z_i*resolution*resolution  ) * 8 + 5 ) = zfin;
                            results->at( ( x_i + y_i*resolution + z_i*resolution*resolution  ) * 8 + 6 ) = Rn;
                            results->at( ( x_i + y_i*resolution + z_i*resolution*resolution  ) * 8 + 7 ) = 1.0;     /// if it is a valid voxel
                        }
                    }
                }
            }
        }
    }
    else if (m_type == POW9){               /// Power Nine
        for (int x_i=0; x_i<resolution; x_i++){
            float x = min + (interval*(float)x_i);
            for (int y_i=0; y_i<resolution; y_i++){
                float y = min + (interval*(float)y_i);
                for (int z_i = 0; z_i<resolution; z_i++){
                    float z = min + (interval*(float)z_i);
                    float xfin = x;
                    float yfin = y;
                    float zfin = z;
                    float cx = x;
                    float cy = y;
                    float cz = z;
                    for (int i=1; i<=iterations; i++){
                        xfin = pow(xfin, 9) - 36*pow(xfin,7)*(yfin*yfin + zfin*zfin) + 126*pow(xfin,5)*pow((yfin*yfin + zfin*zfin),2) - 84*pow(xfin,3)*pow((yfin*yfin + zfin*zfin),3) + 9*xfin*pow((yfin*yfin + zfin*zfin),4) + cx;
                        yfin = pow(yfin, 9) - 36*pow(yfin,7)*(zfin*zfin + xfin*xfin) + 126*pow(yfin,5)*pow((zfin*zfin + xfin*xfin),2) - 84*pow(yfin,3)*pow((zfin*zfin + xfin*xfin),3) + 9*yfin*pow((zfin*zfin + xfin*xfin),4) + cy;
                        zfin = pow(zfin, 9) - 36*pow(zfin,7)*(xfin*xfin + yfin*yfin) + 126*pow(zfin,5)*pow((xfin*xfin + yfin*yfin),2) - 84*pow(zfin,3)*pow((xfin*xfin + yfin*yfin),3) + 9*zfin*pow((xfin*xfin + yfin*yfin),4) + cz;
                        float Rsq = (xfin*xfin)+(yfin*yfin)+(zfin*zfin);
                        if (Rsq > max_r)
                            break;
                        if (i == iterations){
                            results->at( ( x_i + y_i*resolution + z_i*resolution*resolution  ) * 8 ) = x;
                            results->at( ( x_i + y_i*resolution + z_i*resolution*resolution  ) * 8 + 1 ) = y;
                            results->at( ( x_i + y_i*resolution + z_i*resolution*resolution  ) * 8 + 2 ) = z;
                            results->at( ( x_i + y_i*resolution + z_i*resolution*resolution  ) * 8 + 3 ) = xfin;
                            results->at( ( x_i + y_i*resolution + z_i*resolution*resolution  ) * 8 + 4 ) = yfin;
                            results->at( ( x_i + y_i*resolution + z_i*resolution*resolution  ) * 8 + 5 ) = zfin;
                            results->at( ( x_i + y_i*resolution + z_i*resolution*resolution  ) * 8 + 6 ) = sqrt(Rsq);
                            results->at( ( x_i + y_i*resolution + z_i*resolution*resolution  ) * 8 + 7 ) = 1.0;
                        }
                    }
                }
            }
        }
    }
    else if (m_type == QUIN){               /// Quintic
        for (int x_i=0; x_i<resolution; x_i++){
            float x = min + (interval*(float)x_i);
            for (int y_i=0; y_i<resolution; y_i++){
                float y = min + (interval*(float)y_i);
                for (int z_i=0; z_i<resolution; z_i++){
                    float z = min + (interval*(float)z_i);
                    float xfin = x;
                    float yfin = y;
                    float zfin = z;
                    float cx = x;
                    float cy = y;
                    float cz = z;
                    for (int i = 1; i<=iterations; i++){
                        xfin = pow(xfin, 5) - 10*pow(xfin,3)*(yfin*yfin + n_order*yfin*zfin + zfin*zfin) + 5*xfin*(pow(yfin,4) + phase_theta*pow(yfin,3)*zfin + phase_phi*(yfin*yfin)*(zfin*zfin) + phase_theta*yfin*pow(zfin,3) + pow(zfin,4) ) + param_D*(xfin*xfin)*yfin*zfin*(yfin+zfin) + cx;
                        yfin = pow(yfin, 5) - 10*pow(yfin,3)*(zfin*zfin + n_order*xfin*zfin + xfin*xfin) + 5*yfin*(pow(zfin,4) + phase_theta*pow(zfin,3)*xfin + phase_phi*(zfin*zfin)*(xfin*xfin) + phase_theta*zfin*pow(xfin,3) + pow(xfin,4) ) + param_D*(yfin*yfin)*zfin*xfin*(zfin+xfin) + cy;
                        zfin = pow(zfin, 5) - 10*pow(zfin,3)*(xfin*xfin + n_order*xfin*yfin + yfin*yfin) + 5*zfin*(pow(xfin,4) + phase_theta*pow(xfin,3)*yfin + phase_phi*(xfin*xfin)*(yfin*yfin) + phase_theta*xfin*pow(yfin,3) + pow(yfin,4) ) + param_D*(zfin*zfin)*xfin*yfin*(xfin+yfin) + cz;
                        float Rsq = (xfin*xfin)+(yfin*yfin)+(zfin*zfin);
                        if (Rsq > max_r)
                            break;
                        if (i == iterations){
                            results->at( ( x_i + y_i*resolution + z_i*resolution*resolution  ) * 8 ) = x;
                            results->at( ( x_i + y_i*resolution + z_i*resolution*resolution  ) * 8 + 1 ) = y;
                            results->at( ( x_i + y_i*resolution + z_i*resolution*resolution  ) * 8 + 2 ) = z;
                            results->at( ( x_i + y_i*resolution + z_i*resolution*resolution  ) * 8 + 3 ) = xfin;
                            results->at( ( x_i + y_i*resolution + z_i*resolution*resolution  ) * 8 + 4 ) = yfin;
                            results->at( ( x_i + y_i*resolution + z_i*resolution*resolution  ) * 8 + 5 ) = zfin;
                            results->at( ( x_i + y_i*resolution + z_i*resolution*resolution  ) * 8 + 6 ) = sqrt(Rsq);
                            results->at( ( x_i + y_i*resolution + z_i*resolution*resolution  ) * 8 + 7 ) = 1.0;
                        }
                    }
                }
            }
        }
    }
    else if (m_type == CUB){                /// Cubic
        for (int x_i=0; x_i<resolution; x_i++){
            float x = min + (interval*(float)x_i);
            for (int y_i=0; y_i<resolution; y_i++){
                float y = min + (interval*(float)y_i);
                for (int z_i=0; z_i<resolution; z_i++){
                    float z = min + (interval*(float)z_i);
                    float xfin = x;
                    float yfin = y;
                    float zfin = z;
                    float cx = x;
                    float cy = y;
                    float cz = z;
                    for (int i=1; i<=iterations; i++){
                        xfin = pow(xfin, 3) - 3*xfin*( (yfin*yfin) + (zfin*zfin) ) + cx;
                        yfin = - pow(yfin, 3) + 3*yfin*xfin*xfin - yfin*zfin*zfin + cy;
                        zfin = pow(zfin, 3) - 3*zfin*xfin*xfin + zfin*yfin*yfin + cz;
                        float Rsq = (xfin*xfin)+(yfin*yfin)+(zfin*zfin);
                        if (Rsq > max_r)
                            break;
                        if (i == iterations){
                            results->at( ( x_i + y_i*resolution + z_i*resolution*resolution  ) * 8 ) = x;
                            results->at( ( x_i + y_i*resolution + z_i*resolution*resolution  ) * 8 + 1 ) = y;
                            results->at( ( x_i + y_i*resolution + z_i*resolution*resolution  ) * 8 + 2 ) = z;
                            results->at( ( x_i + y_i*resolution + z_i*resolution*resolution  ) * 8 + 3 ) = xfin;
                            results->at( ( x_i + y_i*resolution + z_i*resolution*resolution  ) * 8 + 4 ) = yfin;
                            results->at( ( x_i + y_i*resolution + z_i*resolution*resolution  ) * 8 + 5 ) = zfin;
                            results->at( ( x_i + y_i*resolution + z_i*resolution*resolution  ) * 8 + 6 ) = sqrt(Rsq);
                            results->at( ( x_i + y_i*resolution + z_i*resolution*resolution  ) * 8 + 7 ) = 1.0;
                        }
                    }
                }
            }
        }
    }
    else if (m_type == QUAD){               /// Quadratic
        for (int x_i=0; x_i<resolution; x_i++){
            float x = min + (interval*(float)x_i);
            for (int y_i=0; y_i<resolution; y_i++){
                float y = min + (interval*(float)y_i);
                for (int z_i = 0; z_i<resolution; z_i++){
                    float z = min + (interval*(float)z_i);
                    float xfin = x;
                    float yfin = y;
                    float zfin = z;
                    float cx = x;
                    float cy = y;
                    float cz = z;
                    for (int i=1; i<=iterations; i++){
                        xfin = (xfin*xfin) - (yfin*yfin) + cx;
                        yfin = 2*xfin*zfin + cy;
                        zfin = 2*xfin*yfin + cz;
                        float Rsq = (xfin*xfin)+(yfin*yfin)+(zfin*zfin);
                        if (Rsq > max_r)
                            break;
                        if (i == iterations){
                            results->at( ( x_i + y_i*resolution + z_i*resolution*resolution  ) * 8 ) = x;
                            results->at( ( x_i + y_i*resolution + z_i*resolution*resolution  ) * 8 + 1 ) = y;
                            results->at( ( x_i + y_i*resolution + z_i*resolution*resolution  ) * 8 + 2 ) = z;
                            results->at( ( x_i + y_i*resolution + z_i*resolution*resolution  ) * 8 + 3 ) = xfin;
                            results->at( ( x_i + y_i*resolution + z_i*resolution*resolution  ) * 8 + 4 ) = yfin;
                            results->at( ( x_i + y_i*resolution + z_i*resolution*resolution  ) * 8 + 5 ) = zfin;
                            results->at( ( x_i + y_i*resolution + z_i*resolution*resolution  ) * 8 + 6 ) = sqrt(Rsq);
                            results->at( ( x_i + y_i*resolution + z_i*resolution*resolution  ) * 8 + 7 ) = 1.0;
                        }
                    }
                }
            }
        }
    }
    int c_vox = 0;
    for (int vox_i=0; vox_i<pow(resolution,3); vox_i++){   /// consolidation of good results to get rid of empty voxel data. it puts voxel data close together disregarding invalid voxels
        if (results->at( vox_i*8 + 7 ) == 1.0){
            results->at( c_vox*7 ) = results->at( vox_i*8 );
            results->at( c_vox*7 + 1 ) = results->at( vox_i*8 + 1 );
            results->at( c_vox*7 + 2 ) = results->at( vox_i*8 + 2 );
            results->at( c_vox*7 + 3 ) = results->at( vox_i*8 + 3 );
            results->at( c_vox*7 + 4 ) = results->at( vox_i*8 + 4 );
            results->at( c_vox*7 + 5 ) = results->at( vox_i*8 + 5 );
            results->at( c_vox*7 + 6 ) = results->at( vox_i*8 + 6 );
            c_vox++;
        }
    }
    results->resize(c_vox*7);
    return c_vox;
}

int Export_voxel_d_CPU (float * copy_to){

    for (int i=0; i<results->size(); i++){
        copy_to[i] = results->at(i);
    }
    delete (results);
}

int Calculate_arrays_and_get_number_of_verts_GPU (int device_type, int m_type, int iterations, float Size, int resolution,
                                                  float n_order, float phase_theta, float phase_phi, float max_r, float param_D){

#include "Mandelbulb_gpu_kernels_1.6c.cpp"
    results = new vector <float>;
    platforms = new vector<cl::Platform>;
    devices = new vector<cl::Device>;
    cl::Platform::get(platforms);
    if (platforms->size()==0){
        cout<<"No OpenCL platforms found."<<endl;
        return 0;
    }
    for (int i=0; i<platforms->size(); i++){
        vector<cl::Device> devices_on_platform;
        platforms->at(i).getDevices( CL_DEVICE_TYPE_ALL, &devices_on_platform );
        for (int j=0; j<devices_on_platform.size(); j++)
            devices->push_back(devices_on_platform[j]);
    }
    cout<<"OpenCL devices:"<<endl;
    for (int i=0; i<devices->size(); i++)
        cout<<i<<") "<<devices->at(i).getInfo<CL_DEVICE_NAME>()<<endl;
    cout<<"Remember to choose OCL 3.0 device. For Nvidia GPU it means drivers after 465.xx"<<endl;
    if (device_type>=devices->size()){
        cout<<"Wrong device."<<endl;
        return 0;
    }
    cl::Device *Default_Device = &(devices->at(device_type));
    Cl_Context = new cl::Context({*Default_Device});
    cl::Program::Sources Prog_Sources;
    Prog_Sources.push_back(Jul_Compute.c_str());
    Prog_Sources.push_back(Pow9_Compute.c_str());
    Prog_Sources.push_back(Quin_Compute.c_str());
    Prog_Sources.push_back(Cub_Compute.c_str());
    Prog_Sources.push_back(Quad_Compute.c_str());
    Prog_Sources.push_back(Consol_Compute.c_str());
    Prog_Sources.push_back(Faces_Vertices_FB.c_str());
    Prog = new cl::Program(*Cl_Context, Prog_Sources);
    if (Prog->build(*Default_Device) != CL_SUCCESS){
        cout<<"Error building program"<<endl<<Prog->getBuildInfo<CL_PROGRAM_BUILD_LOG>(*Default_Device)<<endl;
        return 0;
    }
    else cout<<"Success building program."<<endl;
    Cl_Queue = new cl::CommandQueue(*Cl_Context, *Default_Device);
    unsigned int voxel_d_size = pow(resolution,3)* 8;
    results->resize(voxel_d_size, 0.0);
    voxel_d = new cl::Buffer (*Cl_Context, CL_MEM_READ_WRITE, sizeof(float)*voxel_d_size);
    Cl_Queue->enqueueWriteBuffer(*voxel_d, CL_TRUE, 0, sizeof(float)*voxel_d_size, (float *) results->data());
    float interval = Size/(float)resolution;
    float min = (-Size/2.0) + (interval/2.0);
    if (m_type == JUL){                     /// Juliabulb

        cl::Kernel Kernel_Julia=cl::Kernel(*Prog, "Kernel_Julia");
        Kernel_Julia.setArg(0, sizeof(int), &iterations);
        Kernel_Julia.setArg(1, sizeof(int), &resolution);
        Kernel_Julia.setArg(2, sizeof(float), &interval);
        Kernel_Julia.setArg(3, sizeof(float), &min);
        Kernel_Julia.setArg(4, sizeof(float), &n_order);
        Kernel_Julia.setArg(5, sizeof(float), &phase_theta);
        Kernel_Julia.setArg(6, sizeof(float), &phase_phi);
        Kernel_Julia.setArg(7, sizeof(float), &max_r);
        Kernel_Julia.setArg(8, sizeof(float), &param_D);
        Kernel_Julia.setArg(9, *voxel_d);
        Cl_Queue->enqueueNDRangeKernel (Kernel_Julia, cl::NullRange, cl::NDRange(resolution, resolution, resolution), cl::NullRange);
    }
    else if (m_type == POW9){               /// Power Nine

        cl::Kernel Kernel_Pow9=cl::Kernel(*Prog, "Kernel_Pow9");
        Kernel_Pow9.setArg(0, sizeof(int), &iterations);
        Kernel_Pow9.setArg(1, sizeof(int), &resolution);
        Kernel_Pow9.setArg(2, sizeof(float), &interval);
        Kernel_Pow9.setArg(3, sizeof(float), &min);
        Kernel_Pow9.setArg(4, sizeof(float), &n_order);
        Kernel_Pow9.setArg(5, sizeof(float), &phase_theta);
        Kernel_Pow9.setArg(6, sizeof(float), &phase_phi);
        Kernel_Pow9.setArg(7, sizeof(float), &max_r);
        Kernel_Pow9.setArg(8, sizeof(float), &param_D);
        Kernel_Pow9.setArg(9, *voxel_d);
        Cl_Queue->enqueueNDRangeKernel (Kernel_Pow9, cl::NullRange, cl::NDRange(resolution, resolution, resolution), cl::NullRange);
    }
    else if (m_type == QUIN){               /// Quintic

        cl::Kernel Kernel_Quin=cl::Kernel(*Prog, "Kernel_Quin");
        Kernel_Quin.setArg(0, sizeof(int), &iterations);
        Kernel_Quin.setArg(1, sizeof(int), &resolution);
        Kernel_Quin.setArg(2, sizeof(float), &interval);
        Kernel_Quin.setArg(3, sizeof(float), &min);
        Kernel_Quin.setArg(4, sizeof(float), &n_order);
        Kernel_Quin.setArg(5, sizeof(float), &phase_theta);
        Kernel_Quin.setArg(6, sizeof(float), &phase_phi);
        Kernel_Quin.setArg(7, sizeof(float), &max_r);
        Kernel_Quin.setArg(8, sizeof(float), &param_D);
        Kernel_Quin.setArg(9, *voxel_d);
        Cl_Queue->enqueueNDRangeKernel (Kernel_Quin, cl::NullRange, cl::NDRange(resolution, resolution, resolution), cl::NullRange);
    }
    else if (m_type == CUB){                /// Cubic

        cl::Kernel Kernel_Cub=cl::Kernel(*Prog, "Kernel_Cub");
        Kernel_Cub.setArg(0, sizeof(int), &iterations);
        Kernel_Cub.setArg(1, sizeof(int), &resolution);
        Kernel_Cub.setArg(2, sizeof(float), &interval);
        Kernel_Cub.setArg(3, sizeof(float), &min);
        Kernel_Cub.setArg(4, sizeof(float), &n_order);
        Kernel_Cub.setArg(5, sizeof(float), &phase_theta);
        Kernel_Cub.setArg(6, sizeof(float), &phase_phi);
        Kernel_Cub.setArg(7, sizeof(float), &max_r);
        Kernel_Cub.setArg(8, sizeof(float), &param_D);
        Kernel_Cub.setArg(9, *voxel_d);
        Cl_Queue->enqueueNDRangeKernel (Kernel_Cub, cl::NullRange, cl::NDRange(resolution, resolution, resolution), cl::NullRange);
    }
    else if (m_type == QUAD){               /// Quadratic

        cl::Kernel Kernel_Quad=cl::Kernel(*Prog, "Kernel_Quad");
        Kernel_Quad.setArg(0, sizeof(int), &iterations);
        Kernel_Quad.setArg(1, sizeof(int), &resolution);
        Kernel_Quad.setArg(2, sizeof(float), &interval);
        Kernel_Quad.setArg(3, sizeof(float), &min);
        Kernel_Quad.setArg(4, sizeof(float), &n_order);
        Kernel_Quad.setArg(5, sizeof(float), &phase_theta);
        Kernel_Quad.setArg(6, sizeof(float), &phase_phi);
        Kernel_Quad.setArg(7, sizeof(float), &max_r);
        Kernel_Quad.setArg(8, sizeof(float), &param_D);
        Kernel_Quad.setArg(9, *voxel_d);
        Cl_Queue->enqueueNDRangeKernel (Kernel_Quad, cl::NullRange, cl::NDRange(resolution, resolution, resolution), cl::NullRange);
    }
    int c_vox = voxel_d_size/8;
    cl::Buffer voxel_d_s(*Cl_Context, CL_MEM_READ_WRITE, sizeof(int));
    Cl_Queue->enqueueWriteBuffer(voxel_d_s, CL_TRUE, 0, sizeof(int), &c_vox);
    cl::Kernel Kernel_Consol=cl::Kernel(*Prog, "Kernel_Consol");       /// consolidation of good results to get rid of empty voxel data. it puts voxel data close together disregarding invalid voxels
    Kernel_Consol.setArg(0, voxel_d_s);
    Kernel_Consol.setArg(1, *voxel_d);
    Cl_Queue->enqueueNDRangeKernel(Kernel_Consol, cl::NullRange, cl::NDRange(1), cl::NDRange(1));
    Cl_Queue->enqueueReadBuffer(voxel_d_s, CL_TRUE, 0, sizeof(int), &c_vox);
    results->resize(c_vox*7);
    Cl_Queue->enqueueReadBuffer(*voxel_d, CL_TRUE, 0, sizeof(float)*c_vox*7, results->data());
    return c_vox;
}


int Export_voxel_d(float * copy_to){

    for (int i=0; i<results->size(); i++){
        copy_to[i] = results->at(i);
    }
    Cl_Queue->finish();
    delete (voxel_d);
    delete (Cl_Queue);
    delete (Prog);
    delete (Cl_Context);
    delete (devices);
    delete (platforms);
    delete (results);
}

int Export_Faces_Vertices_d (float half_vox_size, float * copy_to_verts, unsigned int * copy_to_faces, float * copy_to_col_loop){

    int Size_f_v = 24*(results->size()/7);
    unsigned int fac[Size_f_v];
    float ver[Size_f_v];
    cl::Buffer faces_d(*Cl_Context, CL_MEM_READ_WRITE, sizeof(unsigned int)*Size_f_v);
    cl::Buffer verts_d(*Cl_Context, CL_MEM_READ_WRITE, sizeof(float)*Size_f_v);
    cl::Kernel Faces_Vertices_FB = cl::Kernel(*Prog, "Faces_Vertices_FB");
    Faces_Vertices_FB.setArg(0, sizeof(float), &half_vox_size);
    Faces_Vertices_FB.setArg(1, *voxel_d);
    Faces_Vertices_FB.setArg(2, faces_d);
    Faces_Vertices_FB.setArg(3, verts_d);
    Cl_Queue->enqueueNDRangeKernel(Faces_Vertices_FB, cl::NullRange, cl::NDRange(results->size()/7), cl::NullRange);
    Cl_Queue->enqueueReadBuffer(faces_d, CL_TRUE, 0, sizeof(unsigned int)*Size_f_v, fac);
    Cl_Queue->enqueueReadBuffer(verts_d, CL_TRUE, 0, sizeof(float)*Size_f_v, ver);
    float vox_buff[results->size()];
    Cl_Queue->enqueueReadBuffer(*voxel_d, CL_TRUE, 0, sizeof(float)*results->size(), vox_buff);
    Cl_Queue->finish();
    for (int i=0; i<results->size()/7; i++){
        copy_to_col_loop[4*i] = results->at(7*i+3);
        copy_to_col_loop[4*i+1] = results->at(7*i+4);
        copy_to_col_loop[4*i+2] = results->at(7*i+5);
        copy_to_col_loop[4*i+3] = results->at(7*i+6);
        for (int j=0; j<24; j++){
            copy_to_faces[24*i+j] = fac[24*i+j];
            copy_to_verts[24*i+j] = ver[24*i+j];
        }
    }
    delete (voxel_d);
    delete (Cl_Queue);
    delete (Prog);
    delete (Cl_Context);
    delete (devices);
    delete (platforms);
    delete (results);
}
