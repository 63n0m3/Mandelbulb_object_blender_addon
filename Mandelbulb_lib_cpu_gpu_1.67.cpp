
/** This file is compilable to usefull DLL that provides mandelbulb generation mechanisms
    for external sources including blender Mandelbulb generation addon
    by Gen0me https://github.com/63n0m3/Mandelbulb_object_blender_addon/
    Licence: MIT free use, distribute and profit with this copyright notice
    Updated 7.11.2023
    btc: bc1q2h9m4pzfdnnfu79hfnjlcnta0j4etjqr4n5gy0
    More information you can find on github in Addon_Mandelbulb**.py comments and wikipedia.org/wiki/Mandelbulb

    Following is the information about compilation and accessing the dll functions:
    You will need to have working OpenCL with search directories pointing to OpenCL-SDK\include for #include <CL/opencl.hpp>
    Add to linker resources: OpenCL.dll
    Add compiler option -c -DBUILD_DLL
    Add linker option -shared -o Mandelbulb_Gen0me.dll
    This should build a dll.

    There are two ways to run following functions. Both ways are independent of each other.

    For single threaded C++ execute in order:
    Calculate_arrays_and_get_number_of_verts_CPU - returns voxel count
    Delete_internal_disconnected_CPU - optional, returns voxel count
    Export_voxel_d_CPU

    For OpenCL execute in order:
    Calculate_arrays_and_get_number_of_verts_GPU - returns voxel count
    Delete_internal_disconnected_GPU - optional, returns voxel count
    Export_voxel_d_GPU

    More information what values represent you will find in python script and in wikipedia page.
    device_type - is an OpenCL device number, probably 0 if you have single GPU
    m_type - is a mandelbulb type - see defines
    di_type, dd_type - are types of cleaning algorithms used: delete internal, delete disconnected - see defines
    Thin skin algorithms are based on comparing voxels to all six surrounding voxels that are face connected like in picture: "voxels calc order.jpg"
    Thick skin algorithms are based on comparing voxels to all 26 surrounding voxels that are face or vertex connected.
    del_disc_crit - is the number of surrounding or in bunch voxels used for those algorithms.
    copy_to - is a pointer to data to be copied to as 3 float coordinate values x,y,z, following 4 float colour values of vertex colour data.
    size available in copy_to should be equal to 7*sizeof(float)*voxel count.
*/

#include <windows.h>
#include <math.h>
#include <CL/opencl.hpp>
#include <vector>
#include <string>
#include <iostream>

///m_type
#define JUL 0
#define CUB 1
#define QUAD 2
#define QUIN 3
#define POW9 4
///di_type
#define NO_DI 0
#define THIN_SKIN 1
#define THICK_SKIN 2
///dd_type
#define NO_DD 0
#define THIN_SKIN_CONNECTED_DIRECTLY 1
#define THICK_SKIN_CONNECTED_DIRECTLY 2
#define THIN_SKIN_CONNECTED_IN_BUNCH 3
#define THICK_SKIN_CONNECTED_IN_BUNCH 4

using namespace std;

float interval;
int voxel_count;
vector <float> * results;
vector<cl::Platform> * platforms;
vector<cl::Device> * devices;
cl::Context * Cl_Context;
cl::Program * Prog;
cl::CommandQueue * Cl_Queue;
cl::Buffer * voxel_d;
cl::Buffer * voxel_d_s;
cl::Kernel * Kernel_Consol_8_to_8;

/// returns valid voxel count and sets up results (voxel data array: x,y,z-coordinates,r,g,b,a-colours )
int Calculate_arrays_and_get_number_of_verts_CPU (int m_type, int iterations, float Size, int resolution, float n_order,
                                                  float phase_theta, float phase_phi, float max_r, float param_D){

    unsigned int voxel_d_size = pow(resolution,3)* 8;
    results = new vector <float>;
    results->resize (voxel_d_size, 0.0);
    interval = Size/(float)resolution;
    float min = -Size/2 + (interval/2);
    if (m_type == JUL){                     /// Juliabulb
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
                for (int z_i=0; z_i<resolution; z_i++){
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
            results->at( c_vox*8 ) = results->at( vox_i*8 );
            results->at( c_vox*8 + 1 ) = results->at( vox_i*8 + 1 );
            results->at( c_vox*8 + 2 ) = results->at( vox_i*8 + 2 );
            results->at( c_vox*8 + 3 ) = results->at( vox_i*8 + 3 );
            results->at( c_vox*8 + 4 ) = results->at( vox_i*8 + 4 );
            results->at( c_vox*8 + 5 ) = results->at( vox_i*8 + 5 );
            results->at( c_vox*8 + 6 ) = results->at( vox_i*8 + 6 );
            results->at( c_vox*8 + 7 ) = 1.0;
            c_vox++;
        }
    }
    results->resize(c_vox*8);
    return c_vox;
}



bool In_vec_vec_int (vector <vector<int>> v_v_int, int int_v){
    for (int i=0; i<v_v_int.size(); i++){
        for (int j=0; j<v_v_int[i].size(); j++){
            if (v_v_int[i][j] == int_v) return true;
        }
    }
    return false;
}
bool In_vec_int (vector<int> v_int, int int_v){
    for (int i=0; i<v_int.size(); i++){
        if (v_int[i] == int_v) return true;
    }
    return false;
}

/// returns valid voxel count and re-sets up results
int Delete_internal_disconnected_CPU (int di_type, int dd_type, int del_disc_crit){

    int verts_count = results->size()/8;
    if (di_type == NO_DI && dd_type == NO_DD)
        return verts_count;

    int skin_bunch;

    if (dd_type == THIN_SKIN_CONNECTED_IN_BUNCH)
        skin_bunch = 6;
    else if (dd_type == THICK_SKIN_CONNECTED_IN_BUNCH)
        skin_bunch = 26;
    else skin_bunch = 0;

    vector <vector <int>> indices_connected (verts_count , vector <int> (skin_bunch, -1));

    for (int i=0; i<verts_count; i++){     /// This loop iterates through every generated voxel and checks how many voxels connect to this voxel
        int count_near_voxels_thick = 0;
        int count_near_voxels_thin = 0;
/// Thin skin algorithms
        if (di_type == THIN_SKIN || dd_type == THIN_SKIN_CONNECTED_DIRECTLY || dd_type == THIN_SKIN_CONNECTED_IN_BUNCH){
            for (int j=0; j<verts_count; j++){

                if (results->at(8*i) == results->at(8*j)){
                    if (results->at(8*i+1) == results->at(8*j+1)){
                        if ( ( ( results->at(8*i+2) < results->at(8*j+2) - interval*0.9 ) && ( results->at(8*i+2) > results->at(8*j+2) - interval*1.1 ) ) || ( ( results->at(8*i+2) > results->at(8*j+2) + interval*0.9 ) && results->at(8*i+2) < results->at(8*j+2) + interval*1.1 ) ){
                            if (dd_type == THIN_SKIN_CONNECTED_IN_BUNCH)
                                indices_connected[i][count_near_voxels_thin] = j;
                            count_near_voxels_thin += 1;
                        }
                    }
                }
                if (results->at(8*i) == results->at(8*j)){
                    if (results->at(8*i+2) == results->at(8*j+2)){
                        if ( ( ( results->at(8*i+1) < results->at(8*j+1) - interval*0.9 ) && ( results->at(8*i+1) > results->at(8*j+1) - interval*1.1 ) ) || ( ( results->at(8*i+1) > results->at(8*j+1) + interval*0.9 ) && results->at(8*i+1) < results->at(8*j+1) + interval*1.1 ) ){
                            if (dd_type == THIN_SKIN_CONNECTED_IN_BUNCH)
                                indices_connected[i][count_near_voxels_thin] = j;
                            count_near_voxels_thin += 1;
                        }
                    }
                }
                if (results->at(8*i+1) == results->at(8*j+1)){
                    if (results->at(8*i+2) == results->at(8*j+2)){
                        if ( ( ( results->at(8*i) < results->at(8*j) - interval*0.9 ) && ( results->at(8*i) > results->at(8*j) - interval*1.1 ) ) || ( ( results->at(8*i) > results->at(8*j) + interval*0.9 ) && results->at(8*i) < results->at(8*j) + interval*1.1 ) ){
                            if (dd_type == THIN_SKIN_CONNECTED_IN_BUNCH)
                                indices_connected[i][count_near_voxels_thin] = j;
                            count_near_voxels_thin += 1;
                        }
                    }
                }
            }
        }
/// Thick skin algorithms
        if (di_type == THICK_SKIN || dd_type == THICK_SKIN_CONNECTED_DIRECTLY || dd_type == THICK_SKIN_CONNECTED_IN_BUNCH){
            for (int j=0; j<verts_count; j++)
                if ( ( results->at(8*i+1) > results->at(8*j+1) - interval*1.1 ) && ( results->at(8*i+1) < results->at(8*j+1) + interval*1.1 ) )
                    if ( ( results->at(8*i+2) > results->at(8*j+2) - interval*1.1 ) && ( results->at(8*i+2) < results->at(8*j+2) + interval*1.1 ) )
                        if ( ( results->at(8*i) > results->at(8*j) - interval*1.1 ) && ( results->at(8*i) < results->at(8*j) + interval*1.1 ) ){
                            if ( ( results->at(8*i) == results->at(8*j) ) && ( results->at(8*i+1) == results->at(8*j+1) ) && ( results->at(8*i+2) == results->at(8*j+2) ) )
                                continue;
                            if (dd_type == THICK_SKIN_CONNECTED_IN_BUNCH)
                                indices_connected[i][count_near_voxels_thick] = j;
                            count_near_voxels_thick += 1;
                        }
        }

        if (di_type == THIN_SKIN)
            if (count_near_voxels_thin == 6)     /// setting internal voxels invalid if there are 6 voxels connecting to this voxel (1 for every face of a cube) in thin skin
                results->at(8*i+7) = 0.0;

        if (di_type == THICK_SKIN)
            if (count_near_voxels_thick == 26)
                results->at(8*i+7) = 0.0;

        if (dd_type == THIN_SKIN_CONNECTED_DIRECTLY){  /// So here comes the fast algorithm. Its criterion is the number of surrounding voxels
            if (count_near_voxels_thin <= del_disc_crit)
                results->at(8*i+7) = 0.0;
        }
        if (dd_type == THICK_SKIN_CONNECTED_DIRECTLY){
            if (count_near_voxels_thick <= del_disc_crit)
                results->at(8*i+7) = 0.0;
        }
    }

    /// Now inside indices_connected[n][0-5 or 0-26] should be indices of connected verts (for bunch)

    if (dd_type==THIN_SKIN_CONNECTED_IN_BUNCH || dd_type==THICK_SKIN_CONNECTED_IN_BUNCH){             /// And this is the slow algorithm. It groups voxels and checks how many are connected in each group
                                                                                                    /// criterion is nubber of voxels in each group
        vector <vector <int>> voxel_group;

        int group_counter = 0;

        for (int i=0; i<verts_count; i++){

            if (In_vec_vec_int(voxel_group, i) == false){       ///this adds indice to voxel group if indice was not grouped yet
                voxel_group.push_back(vector<int>());
                voxel_group[group_counter].push_back(i);
            }
            else continue;

            unsigned int i_group=0;
            while (i_group < voxel_group[group_counter].size()){                   ///this adds next indices to voxel group from the list of connected voxels to checked voxel group
                for (int j=0; j<skin_bunch; j++){
                    if (indices_connected [ voxel_group[group_counter] [i_group] ] [j] != -1){
                        if (In_vec_int (voxel_group[group_counter], indices_connected [ voxel_group[group_counter] [i_group] ] [j]) == false){
                            voxel_group[group_counter].push_back (indices_connected [ voxel_group[group_counter] [i_group] ] [j]);
                        }
                    }
                }
                i_group++;
            }
            group_counter++;
        }

        for (int i=0; i<voxel_group.size(); i++)               /// set voxels for deletion
            if (voxel_group[i].size() <= del_disc_crit)
                for (int j=0; j<voxel_group[i].size(); j++)
                    results->at(8*(voxel_group[i][j])+7) = 0.0;
    }

    int c_vox = 0;
    for (int vox_i=0; vox_i<results->size()/8; vox_i++){   /// consolidation of results
        if (results->at( vox_i*8 + 7 ) == 1.0){
            results->at( c_vox*8 ) = results->at( vox_i*8 );
            results->at( c_vox*8 + 1 ) = results->at( vox_i*8 + 1 );
            results->at( c_vox*8 + 2 ) = results->at( vox_i*8 + 2 );
            results->at( c_vox*8 + 3 ) = results->at( vox_i*8 + 3 );
            results->at( c_vox*8 + 4 ) = results->at( vox_i*8 + 4 );
            results->at( c_vox*8 + 5 ) = results->at( vox_i*8 + 5 );
            results->at( c_vox*8 + 6 ) = results->at( vox_i*8 + 6 );
            results->at( c_vox*8 + 7 ) = 1.0;
            c_vox++;
        }
    }

    results->resize(c_vox*8);

    return c_vox;
}


int Export_voxel_d_CPU (float * copy_to){

    for (int i=0; i<results->size()/8; i++){
        copy_to[7*i] = results->at(8*i);
        copy_to[7*i+1] = results->at(8*i+1);
        copy_to[7*i+2] = results->at(8*i+2);
        copy_to[7*i+3] = results->at(8*i+3);
        copy_to[7*i+4] = results->at(8*i+4);
        copy_to[7*i+5] = results->at(8*i+5);
        copy_to[7*i+6] = results->at(8*i+6);
    }
    delete (results);
}

int Calculate_arrays_and_get_number_of_verts_GPU (int device_type, int m_type, int iterations, float Size, int resolution,
                                                  float n_order, float phase_theta, float phase_phi, float max_r, float param_D){

#include "Mandelbulb_gpu_kernels_1.67.cpp"
    platforms = new vector<cl::Platform>;
    devices = new vector<cl::Device>;
    cl::Platform::get(platforms);
    if (platforms->size() == 0){
        cout<<"No OpenCL platforms found."<<endl;
        return 0;
    }
    for (int i=0; i<platforms->size(); i++){
        vector<cl::Device> devices_on_platform;
        platforms->at(i).getDevices (CL_DEVICE_TYPE_ALL, &devices_on_platform);
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

    cout<<"MAX WORKGROUP SIZE: "<<Default_Device->getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>()<<endl;

    Cl_Context = new cl::Context({*Default_Device});
    cl::Program::Sources Prog_Sources;
    Prog_Sources.push_back(Jul_Compute.c_str());
    Prog_Sources.push_back(Pow9_Compute.c_str());
    Prog_Sources.push_back(Quin_Compute.c_str());
    Prog_Sources.push_back(Cub_Compute.c_str());
    Prog_Sources.push_back(Quad_Compute.c_str());
    Prog_Sources.push_back(Consol_Compute_8_to_7.c_str());
    Prog_Sources.push_back(Consol_Compute_8_to_8.c_str());
    Prog_Sources.push_back(Delete_Internal_Disconnected_Compute.c_str());
    Prog_Sources.push_back(Delete_Disconnected_Bunch_Compute.c_str());
    Prog_Sources.push_back(Delete_Internal_Compute.c_str());    /// depreciated
    Prog = new cl::Program(*Cl_Context, Prog_Sources);
    if (Prog->build(*Default_Device) != CL_SUCCESS){
        cout<<"Error building program"<<endl<<Prog->getBuildInfo<CL_PROGRAM_BUILD_LOG>(*Default_Device)<<endl;
        return 0;
    }
    else cout<<"Success building program."<<endl;
    Cl_Queue = new cl::CommandQueue (*Cl_Context, *Default_Device);
    unsigned int voxel_d_size = pow(resolution,3) * 8;
    voxel_d = new cl::Buffer (*Cl_Context, CL_MEM_READ_WRITE, sizeof(float)*voxel_d_size);
    Cl_Queue->enqueueFillBuffer (*voxel_d, (float)0.0, 0, sizeof(float)*voxel_d_size);

    Cl_Queue->finish();
    interval = Size/(float)resolution;
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

    Cl_Queue->finish();
    voxel_count = voxel_d_size/8;
    voxel_d_s = new cl::Buffer (*Cl_Context, CL_MEM_READ_WRITE, sizeof(int));
    Cl_Queue->enqueueWriteBuffer (*voxel_d_s, CL_TRUE, 0, sizeof(int), &voxel_count);
    Cl_Queue->finish();
    Kernel_Consol_8_to_8 = new cl::Kernel (*Prog, "Kernel_Consol_8_to_8");       /// consolidation of good results to get rid of empty voxel data. it puts voxel data close together disregarding invalid voxels
    Kernel_Consol_8_to_8->setArg(0, *voxel_d_s);
    Kernel_Consol_8_to_8->setArg(1, *voxel_d);

    int error = Cl_Queue->enqueueNDRangeKernel (*Kernel_Consol_8_to_8, cl::NullRange, cl::NDRange(1), cl::NDRange(1));
    if (error == CL_SUCCESS) cout<<"CL_SUCCESS\n";
    else if (error == CL_OUT_OF_RESOURCES) cout<<"CL_OUT_OF_RESOURCES\n";
    else cout<<"CL_ERROR\n";

    Cl_Queue->finish();

    error = Cl_Queue->enqueueReadBuffer (*voxel_d_s, CL_TRUE, 0, sizeof(int), &voxel_count);
    if (error == CL_SUCCESS) cout<<"CL_SUCCESS\n";
    else if (error == CL_OUT_OF_RESOURCES) cout<<"CL_OUT_OF_RESOURCES\n";
    else cout<<"CL_ERROR\n";

    Cl_Queue->finish();
    return voxel_count;
}


int Delete_internal_disconnected_GPU (int di_type, int dd_type, int del_disc_crit){

    int error;

    if (di_type == NO_DI && dd_type==NO_DD)
        goto DIDGPU;

    int skin_bunch;

    if (dd_type == THIN_SKIN_CONNECTED_IN_BUNCH)
        skin_bunch = 6;
    else if (dd_type == THICK_SKIN_CONNECTED_IN_BUNCH)
        skin_bunch = 26;
    else skin_bunch = 0;

    {
        cl::Buffer indices (*Cl_Context, CL_MEM_READ_WRITE, sizeof(int)*skin_bunch*voxel_count);
        Cl_Queue->enqueueFillBuffer (indices, (int)-1, 0, sizeof(int)*skin_bunch*voxel_count);

        cl::Kernel Kernel_Delete_Internal_Disconnected = cl::Kernel (*Prog, "Kernel_Delete_Internal_Disconnected");
        Kernel_Delete_Internal_Disconnected.setArg(0, *voxel_d_s);
        Kernel_Delete_Internal_Disconnected.setArg(1, *voxel_d);
        Kernel_Delete_Internal_Disconnected.setArg(2, sizeof(float), &interval);
        Kernel_Delete_Internal_Disconnected.setArg(3, sizeof(int), &di_type);
        Kernel_Delete_Internal_Disconnected.setArg(4, sizeof(int), &dd_type);
        Kernel_Delete_Internal_Disconnected.setArg(5, sizeof(int), &del_disc_crit);
        Kernel_Delete_Internal_Disconnected.setArg(6, indices);
        Kernel_Delete_Internal_Disconnected.setArg(7, sizeof(int), &skin_bunch);
        Cl_Queue->enqueueNDRangeKernel (Kernel_Delete_Internal_Disconnected, cl::NullRange, cl::NDRange (voxel_count), cl::NullRange);

        if (dd_type == THIN_SKIN_CONNECTED_IN_BUNCH || dd_type == THICK_SKIN_CONNECTED_IN_BUNCH){
            cl::Buffer groups (*Cl_Context, CL_MEM_READ_WRITE, sizeof(int)*voxel_count*voxel_count);
            cl::Kernel Kernel_Delete_Disconnected_Bunch = cl::Kernel (*Prog, "Kernel_Delete_Disconnected_Bunch");
            Kernel_Delete_Disconnected_Bunch.setArg(0, *voxel_d_s);
            Kernel_Delete_Disconnected_Bunch.setArg(1, *voxel_d);
            Kernel_Delete_Disconnected_Bunch.setArg(2, sizeof(int), &del_disc_crit);
            Kernel_Delete_Disconnected_Bunch.setArg(3, indices);
            Kernel_Delete_Disconnected_Bunch.setArg(4, sizeof(int), &skin_bunch);
            Kernel_Delete_Disconnected_Bunch.setArg(5, groups);
            Cl_Queue->enqueueNDRangeKernel (Kernel_Delete_Disconnected_Bunch, cl::NullRange, cl::NDRange (voxel_count), cl::NullRange);
        }
    }

    Kernel_Consol_8_to_8->setArg(0, *voxel_d_s);
    Kernel_Consol_8_to_8->setArg(1, *voxel_d);

    error = Cl_Queue->enqueueNDRangeKernel(*Kernel_Consol_8_to_8, cl::NullRange, cl::NDRange(1), cl::NDRange(1));
    if (error == CL_SUCCESS) cout<<"CL_SUCCESS\n";
    else if (error == CL_OUT_OF_RESOURCES) cout<<"CL_OUT_OF_RESOURCES\n";
    else cout<<"CL_ERROR\n";

    Cl_Queue->finish();

DIDGPU:

    error = Cl_Queue->enqueueReadBuffer(*voxel_d_s, CL_TRUE, 0, sizeof(int), &voxel_count);
    if (error == CL_SUCCESS) cout<<"CL_SUCCESS\n";
    else if (error == CL_OUT_OF_RESOURCES) cout<<"CL_OUT_OF_RESOURCES\n";
    else cout<<"CL_ERROR\n";
    Cl_Queue->finish();

    return voxel_count;
}



int Delete_internal_GPU (){       /// this is depreciated algorithm. it deletes internal according to thin skin algorithm
    cl::Kernel Kernel_Delete_Internal = cl::Kernel (*Prog, "Kernel_Delete_Internal");
    Kernel_Delete_Internal.setArg(0, *voxel_d_s);
    Kernel_Delete_Internal.setArg(1, *voxel_d);
    Kernel_Delete_Internal.setArg(2, sizeof(float), &interval);
    Cl_Queue->enqueueNDRangeKernel (Kernel_Delete_Internal, cl::NullRange, cl::NDRange (voxel_count), cl::NullRange);
    Cl_Queue->finish();
    Kernel_Consol_8_to_8->setArg(0, *voxel_d_s);
    Kernel_Consol_8_to_8->setArg(1, *voxel_d);
    Cl_Queue->enqueueNDRangeKernel (*Kernel_Consol_8_to_8, cl::NullRange, cl::NDRange(1), cl::NDRange(1));
    Cl_Queue->finish();
    Cl_Queue->enqueueReadBuffer (*voxel_d_s, CL_TRUE, 0, sizeof(int), &voxel_count);
    Cl_Queue->finish();
    return voxel_count;
}

int Export_voxel_d_GPU (float * copy_to){
    cl::Kernel Kernel_Consol_8_to_7 = cl::Kernel (*Prog, "Kernel_Consol_8_to_7");       /// consolidation 8 component voxel data into 7 component
    Kernel_Consol_8_to_7.setArg(0, *voxel_d_s);
    Kernel_Consol_8_to_7.setArg(1, *voxel_d);
    Cl_Queue->enqueueNDRangeKernel (Kernel_Consol_8_to_7, cl::NullRange, cl::NDRange(1), cl::NDRange(1));
    Cl_Queue->finish();
    Cl_Queue->enqueueReadBuffer (*voxel_d_s, CL_TRUE, 0, sizeof(int), &voxel_count);
    Cl_Queue->finish();
    Cl_Queue->enqueueReadBuffer (*voxel_d, CL_TRUE, 0, sizeof(float)*voxel_count*7, copy_to);
    Cl_Queue->finish();

    delete (Kernel_Consol_8_to_8);
    delete (voxel_d_s);
    delete (voxel_d);
    delete (Cl_Queue);
    delete (Prog);
    delete (Cl_Context);
    delete (devices);
    delete (platforms);
}


main(){




  ///      This works
        Calculate_arrays_and_get_number_of_verts_GPU (0, 3, 2, 2.5, 70,
                                                  5, 0, 0, 8, 0);
        int voxel_count = Delete_internal_disconnected_GPU (1, 1, 2);
        float * voxel_d = (float*)malloc (7*voxel_count*sizeof(float));
        Export_voxel_d_GPU (voxel_d);


  /// This breaks on 673
  /*      Calculate_arrays_and_get_number_of_verts_GPU (0, 3, 2, 2.5, 100,
                                                  5, 0, 0, 8, 0);
        int voxel_count = Delete_internal_disconnected_GPU (1, 1, 2);
        float * voxel_d = (float*)malloc (7*voxel_count*sizeof(float));
        Export_voxel_d_GPU (voxel_d);
    */
  /// This breaks on 608
  /*      Calculate_arrays_and_get_number_of_verts_GPU (0, 3, 2, 2.5, 200,
                                                  5, 0, 0, 8, 0);
        float * voxel_d = (float*)malloc (7*voxel_count*sizeof(float));
        Export_voxel_d_GPU (voxel_d);
    */
}
