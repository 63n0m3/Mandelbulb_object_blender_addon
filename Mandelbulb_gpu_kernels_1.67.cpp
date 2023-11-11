string Jul_Compute = R""""(
void kernel Kernel_Julia (const int iterations, const int resolution, const float interval, const float min, const float n_order,
                          const float phase_theta, const float phase_phi, const float max_r, const float param_D, global float *voxel_d){

    float x = min + (interval*(float)get_global_id(0));         /// Algorithm works this way that it iterates through every voxel and checks if it is inside mandelbulb
    float y = min + (interval*(float)get_global_id(1));         /// Those 3 workgroups iterate through every x,y,z of the voxels
    float z = min + (interval*(float)get_global_id(2));
    float xfin = x;                 /// Asignment of the starting vector
    float yfin = y;
    float zfin = z;
    float cx = x;
    float cy = y;
    float cz = z;
    for (int i=1; i<=iterations; i++){                     /// This loop reiterates mandelbulb math updating new x,y,z vertor values
        float Rsq = (xfin*xfin)+(yfin*yfin)+(zfin*zfin);
        float R = sqrt(Rsq); /// Length of a vector
        float phi = atan(yfin/xfin) + phase_phi;
        float theta = atan(sqrt( (xfin*xfin)+(yfin*yfin) )/zfin) + phase_theta;
        float Rn = pow (R,n_order);
        xfin = sin(n_order*theta) * cos(n_order*phi);
        yfin = sin(n_order*theta) * sin(n_order*phi);
        zfin = cos(n_order*theta);
        xfin *= Rn + cx;        /// Addition of the starting vector
        yfin *= Rn + cy;
        zfin *= Rn + cz;
        if (Rsq > max_r)                                        /// This is condition to check if voxel is inside mandelbulb
            break;
        if (i == iterations){                                   /// This checks if i-loop has gone the necessary number of iterations. If it went through it means the vertex is valid
            voxel_d[ ( get_global_id(0) + get_global_id(1)*resolution + get_global_id(2)*resolution*resolution  ) * 8 ] = x;           ///position
            voxel_d[ ( get_global_id(0) + get_global_id(1)*resolution + get_global_id(2)*resolution*resolution  ) * 8 + 1 ] = y;
            voxel_d[ ( get_global_id(0) + get_global_id(1)*resolution + get_global_id(2)*resolution*resolution  ) * 8 + 2 ] = z;
            voxel_d[ ( get_global_id(0) + get_global_id(1)*resolution + get_global_id(2)*resolution*resolution  ) * 8 + 3 ] = xfin;    ///colors
            voxel_d[ ( get_global_id(0) + get_global_id(1)*resolution + get_global_id(2)*resolution*resolution  ) * 8 + 4 ] = yfin;
            voxel_d[ ( get_global_id(0) + get_global_id(1)*resolution + get_global_id(2)*resolution*resolution  ) * 8 + 5 ] = zfin;
            voxel_d[ ( get_global_id(0) + get_global_id(1)*resolution + get_global_id(2)*resolution*resolution  ) * 8 + 6 ] = Rn;
            voxel_d[ ( get_global_id(0) + get_global_id(1)*resolution + get_global_id(2)*resolution*resolution  ) * 8 + 7 ] = 1.0;     /// if it is a valid voxel
        }
    }
}
)"""";
string Pow9_Compute = R""""(
void kernel Kernel_Pow9 (const int iterations, const int resolution, const float interval, const float min, const float n_order,
                        const float phase_theta, const float phase_phi, const float max_r, const float param_D, global float *voxel_d){

        float x = min + (interval*(float)get_global_id(0));
        float y = min + (interval*(float)get_global_id(1));
        float z = min + (interval*(float)get_global_id(2));
        float xfin = x;
        float yfin = y;
        float zfin = z;
        float cx = x;
        float cy = y;
        float cz = z;
        for (int i=1; i<=iterations; i++){
            xfin = pow (xfin, 9) - 36*pow(xfin,7)*(yfin*yfin + zfin*zfin) + 126*pow(xfin,5)*pow((yfin*yfin + zfin*zfin),2) - 84*pow(xfin,3)*pow((yfin*yfin + zfin*zfin),3) + 9*xfin*pow((yfin*yfin + zfin*zfin),4) + cx;
            yfin = pow (yfin, 9) - 36*pow(yfin,7)*(zfin*zfin + xfin*xfin) + 126*pow(yfin,5)*pow((zfin*zfin + xfin*xfin),2) - 84*pow(yfin,3)*pow((zfin*zfin + xfin*xfin),3) + 9*yfin*pow((zfin*zfin + xfin*xfin),4) + cy;
            zfin = pow (zfin, 9) - 36*pow(zfin,7)*(xfin*xfin + yfin*yfin) + 126*pow(zfin,5)*pow((xfin*xfin + yfin*yfin),2) - 84*pow(zfin,3)*pow((xfin*xfin + yfin*yfin),3) + 9*zfin*pow((xfin*xfin + yfin*yfin),4) + cz;
            float Rsq = (xfin*xfin)+(yfin*yfin)+(zfin*zfin);
            if (Rsq > max_r)
                break;
            if (i == iterations){
                voxel_d[ ( get_global_id(0) + get_global_id(1)*resolution + get_global_id(2)*resolution*resolution  ) * 8 ] = x;
                voxel_d[ ( get_global_id(0) + get_global_id(1)*resolution + get_global_id(2)*resolution*resolution  ) * 8 + 1 ] = y;
                voxel_d[ ( get_global_id(0) + get_global_id(1)*resolution + get_global_id(2)*resolution*resolution  ) * 8 + 2 ] = z;
                voxel_d[ ( get_global_id(0) + get_global_id(1)*resolution + get_global_id(2)*resolution*resolution  ) * 8 + 3 ] = xfin;
                voxel_d[ ( get_global_id(0) + get_global_id(1)*resolution + get_global_id(2)*resolution*resolution  ) * 8 + 4 ] = yfin;
                voxel_d[ ( get_global_id(0) + get_global_id(1)*resolution + get_global_id(2)*resolution*resolution  ) * 8 + 5 ] = zfin;
                voxel_d[ ( get_global_id(0) + get_global_id(1)*resolution + get_global_id(2)*resolution*resolution  ) * 8 + 6 ] = sqrt(Rsq);
                voxel_d[ ( get_global_id(0) + get_global_id(1)*resolution + get_global_id(2)*resolution*resolution  ) * 8 + 7 ] = 1.0;
            }
        }
}
)"""";
string Quin_Compute = R""""(
void kernel Kernel_Quin (const int iterations, const int resolution, const float interval, const float min, const float n_order,
                         const float phase_theta, const float phase_phi, const float max_r, const float param_D, global float *voxel_d){

        float x = min + (interval*(float)get_global_id(0));
        float y = min + (interval*(float)get_global_id(1));
        float z = min + (interval*(float)get_global_id(2));
        float xfin = x;
        float yfin = y;
        float zfin = z;
        float cx = x;
        float cy = y;
        float cz = z;
        for(int i = 1; i <= iterations; i++){
            xfin = pow (xfin, 5) - 10*pow(xfin,3)*(yfin*yfin + n_order*yfin*zfin + zfin*zfin) + 5*xfin*(pow(yfin,4) + phase_theta*pow(yfin,3)*zfin + phase_phi*(yfin*yfin)*(zfin*zfin) + phase_theta*yfin*pow(zfin,3) + pow(zfin,4) ) + param_D*(xfin*xfin)*yfin*zfin*(yfin+zfin) + cx;
            yfin = pow (yfin, 5) - 10*pow(yfin,3)*(zfin*zfin + n_order*xfin*zfin + xfin*xfin) + 5*yfin*(pow(zfin,4) + phase_theta*pow(zfin,3)*xfin + phase_phi*(zfin*zfin)*(xfin*xfin) + phase_theta*zfin*pow(xfin,3) + pow(xfin,4) ) + param_D*(yfin*yfin)*zfin*xfin*(zfin+xfin) + cy;
            zfin = pow (zfin, 5) - 10*pow(zfin,3)*(xfin*xfin + n_order*xfin*yfin + yfin*yfin) + 5*zfin*(pow(xfin,4) + phase_theta*pow(xfin,3)*yfin + phase_phi*(xfin*xfin)*(yfin*yfin) + phase_theta*xfin*pow(yfin,3) + pow(yfin,4) ) + param_D*(zfin*zfin)*xfin*yfin*(xfin+yfin) + cz;
            float Rsq = (xfin*xfin)+(yfin*yfin)+(zfin*zfin);
            if (Rsq > max_r)
                break;
            if (i == iterations){
                voxel_d[ ( get_global_id(0) + get_global_id(1)*resolution + get_global_id(2)*resolution*resolution  ) * 8 ] = x;
                voxel_d[ ( get_global_id(0) + get_global_id(1)*resolution + get_global_id(2)*resolution*resolution  ) * 8 + 1 ] = y;
                voxel_d[ ( get_global_id(0) + get_global_id(1)*resolution + get_global_id(2)*resolution*resolution  ) * 8 + 2 ] = z;
                voxel_d[ ( get_global_id(0) + get_global_id(1)*resolution + get_global_id(2)*resolution*resolution  ) * 8 + 3 ] = xfin;
                voxel_d[ ( get_global_id(0) + get_global_id(1)*resolution + get_global_id(2)*resolution*resolution  ) * 8 + 4 ] = yfin;
                voxel_d[ ( get_global_id(0) + get_global_id(1)*resolution + get_global_id(2)*resolution*resolution  ) * 8 + 5 ] = zfin;
                voxel_d[ ( get_global_id(0) + get_global_id(1)*resolution + get_global_id(2)*resolution*resolution  ) * 8 + 6 ] = sqrt(Rsq);
                voxel_d[ ( get_global_id(0) + get_global_id(1)*resolution + get_global_id(2)*resolution*resolution  ) * 8 + 7 ] = 1.0;
            }
        }
}
)"""";
string Cub_Compute = R""""(
void kernel Kernel_Cub (const int iterations, const int resolution, const float interval, const float min, const float n_order,
                        const float phase_theta, const float phase_phi, const float max_r, const float param_D, global float *voxel_d){

        float x = min + (interval*(float)get_global_id(0));
        float y = min + (interval*(float)get_global_id(1));
        float z = min + (interval*(float)get_global_id(2));
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
                voxel_d[ ( get_global_id(0) + get_global_id(1)*resolution + get_global_id(2)*resolution*resolution  ) * 8 ] = x;
                voxel_d[ ( get_global_id(0) + get_global_id(1)*resolution + get_global_id(2)*resolution*resolution  ) * 8 + 1 ] = y;
                voxel_d[ ( get_global_id(0) + get_global_id(1)*resolution + get_global_id(2)*resolution*resolution  ) * 8 + 2 ] = z;
                voxel_d[ ( get_global_id(0) + get_global_id(1)*resolution + get_global_id(2)*resolution*resolution  ) * 8 + 3 ] = xfin;
                voxel_d[ ( get_global_id(0) + get_global_id(1)*resolution + get_global_id(2)*resolution*resolution  ) * 8 + 4 ] = yfin;
                voxel_d[ ( get_global_id(0) + get_global_id(1)*resolution + get_global_id(2)*resolution*resolution  ) * 8 + 5 ] = zfin;
                voxel_d[ ( get_global_id(0) + get_global_id(1)*resolution + get_global_id(2)*resolution*resolution  ) * 8 + 6 ] = sqrt(Rsq);
                voxel_d[ ( get_global_id(0) + get_global_id(1)*resolution + get_global_id(2)*resolution*resolution  ) * 8 + 7 ] = 1.0;
            }
        }
}
)"""";
string Quad_Compute = R""""(
void kernel Kernel_Quad (const int iterations, const int resolution, const float interval, const float min, const float n_order,
                         const float phase_theta, const float phase_phi, const float max_r, const float param_D, global float *voxel_d){

        float x = min + (interval*(float)get_global_id(0));
        float y = min + (interval*(float)get_global_id(1));
        float z = min + (interval*(float)get_global_id(2));
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
                voxel_d[ ( get_global_id(0) + get_global_id(1)*resolution + get_global_id(2)*resolution*resolution  ) * 8 ] = x;
                voxel_d[ ( get_global_id(0) + get_global_id(1)*resolution + get_global_id(2)*resolution*resolution  ) * 8 + 1 ] = y;
                voxel_d[ ( get_global_id(0) + get_global_id(1)*resolution + get_global_id(2)*resolution*resolution  ) * 8 + 2 ] = z;
                voxel_d[ ( get_global_id(0) + get_global_id(1)*resolution + get_global_id(2)*resolution*resolution  ) * 8 + 3 ] = xfin;
                voxel_d[ ( get_global_id(0) + get_global_id(1)*resolution + get_global_id(2)*resolution*resolution  ) * 8 + 4 ] = yfin;
                voxel_d[ ( get_global_id(0) + get_global_id(1)*resolution + get_global_id(2)*resolution*resolution  ) * 8 + 5 ] = zfin;
                voxel_d[ ( get_global_id(0) + get_global_id(1)*resolution + get_global_id(2)*resolution*resolution  ) * 8 + 6 ] = sqrt(Rsq);
                voxel_d[ ( get_global_id(0) + get_global_id(1)*resolution + get_global_id(2)*resolution*resolution  ) * 8 + 7 ] = 1.0;
            }
        }
}
)"""";

string Consol_Compute_8_to_8 = R""""(
void kernel Kernel_Consol_8_to_8 (global int *c_vox, global float *voxel_d){

    if (get_global_id(0) == 0){
        int vox_d_n = 0;
        for (int vox_i=0; vox_i<*c_vox; vox_i++){
            if (voxel_d[ vox_i*8 + 7 ] == 1.0){
                voxel_d[ vox_d_n*8 ] = voxel_d[ vox_i*8 ];
                voxel_d[ vox_d_n*8 + 1 ] = voxel_d[ vox_i*8 + 1 ];
                voxel_d[ vox_d_n*8 + 2 ] = voxel_d[ vox_i*8 + 2 ];
                voxel_d[ vox_d_n*8 + 3 ] = voxel_d[ vox_i*8 + 3 ];
                voxel_d[ vox_d_n*8 + 4 ] = voxel_d[ vox_i*8 + 4 ];
                voxel_d[ vox_d_n*8 + 5 ] = voxel_d[ vox_i*8 + 5 ];
                voxel_d[ vox_d_n*8 + 6 ] = voxel_d[ vox_i*8 + 6 ];
                voxel_d[ vox_d_n*8 + 7 ] = 1.0;
                vox_d_n++;
            }
        }
        *c_vox = vox_d_n;
    }
}
)"""";

string Consol_Compute_8_to_7 = R""""(
void kernel Kernel_Consol_8_to_7 (global int *c_vox, global float *voxel_d){

    if (get_global_id(0) == 0){
        int vox_d_n = 0;
        for (int vox_i=0; vox_i<*c_vox; vox_i++){
            if (voxel_d[ vox_i*8 + 7 ] == 1.0){
                voxel_d[ vox_d_n*7 ] = voxel_d[ vox_i*8 ];
                voxel_d[ vox_d_n*7 + 1 ] = voxel_d[ vox_i*8 + 1 ];
                voxel_d[ vox_d_n*7 + 2 ] = voxel_d[ vox_i*8 + 2 ];
                voxel_d[ vox_d_n*7 + 3 ] = voxel_d[ vox_i*8 + 3 ];
                voxel_d[ vox_d_n*7 + 4 ] = voxel_d[ vox_i*8 + 4 ];
                voxel_d[ vox_d_n*7 + 5 ] = voxel_d[ vox_i*8 + 5 ];
                voxel_d[ vox_d_n*7 + 6 ] = voxel_d[ vox_i*8 + 6 ];
                vox_d_n++;
            }
        }
        *c_vox = vox_d_n;
    }
}
)"""";

string Delete_Internal_Compute = R""""(
void kernel Kernel_Delete_Internal (global int *c_vox, global float *voxel_d, const float interval){
                                                /// This iterates through every generated voxel and checks if it is internal voxel
        int i = get_global_id(0);
        int count_near_voxels = 0;

        for (int j=0; j<*c_vox; j++){
            if (voxel_d[8*i] == voxel_d[8*j]){
                if (voxel_d[8*i+1] == voxel_d[8*j+1]){
                    if ( ( ( voxel_d[8*i+2] < voxel_d[8*j+2] - interval*0.9 ) && ( voxel_d[8*i+2] > voxel_d[8*j+2] - interval*1.1 ) ) || ( ( voxel_d[8*i+2] > voxel_d[8*j+2] + interval*0.9 ) && voxel_d[8*i+2] < voxel_d[8*j+2] + interval*1.1 ) ){
                        count_near_voxels += 1;
                        if (count_near_voxels == 6)                                 /// invalidate it
                            voxel_d[i*8 + 7] = 0.0;
                    }
                }
            }
            if (voxel_d[8*i] == voxel_d[8*j]){
                if (voxel_d[8*i+2] == voxel_d[8*j+2]){
                    if ( ( ( voxel_d[8*i+1] < voxel_d[8*j+1] - interval*0.9 ) && ( voxel_d[8*i+1] > voxel_d[8*j+1] - interval*1.1 ) ) || ( ( voxel_d[8*i+1] > voxel_d[8*j+1] + interval*0.9 ) && voxel_d[8*i+1] < voxel_d[8*j+1] + interval*1.1 ) )
                        count_near_voxels += 1;
                        if (count_near_voxels == 6)                                 /// invalidate it
                            voxel_d[i*8 + 7] = 0.0;
                }
            }
            if (voxel_d[8*i+1] == voxel_d[8*j+1]){
                if (voxel_d[8*i+2] == voxel_d[8*j+2]){
                    if ( ( ( voxel_d[8*i] < voxel_d[8*j] - interval*0.9 ) && ( voxel_d[8*i] > voxel_d[8*j] - interval*1.1 ) ) || ( ( voxel_d[8*i] > voxel_d[8*j] + interval*0.9 ) && voxel_d[8*i] < voxel_d[8*j] + interval*1.1 ) )
                        count_near_voxels += 1;
                        if (count_near_voxels == 6)                                 /// invalidate it
                            voxel_d[i*8 + 7] = 0.0;
                }
            }
        }
}
)"""";



string Delete_Internal_Disconnected_Compute = R""""(
void kernel Kernel_Delete_Internal_Disconnected (const global int *c_vox, global float *voxel_d, const float interval, const int di_type, const int dd_type, const int del_disc_crit, global int *indices_connected, const int skin_bunch){

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

    int i = get_global_id(0);

    int count_near_voxels_thick = 0;
    int count_near_voxels_thin = 0;
/// Thin skin algorithms
    if (di_type == THIN_SKIN || dd_type == THIN_SKIN_CONNECTED_DIRECTLY || dd_type == THIN_SKIN_CONNECTED_IN_BUNCH){
        for (int j=0; j<*c_vox; j++){
            if (voxel_d[8*i] == voxel_d[8*j]){
                if (voxel_d[8*i+1] == voxel_d[8*j+1]){
                    if ( ( ( voxel_d[8*i+2] < voxel_d[8*j+2] - interval*0.9 ) && ( voxel_d[8*i+2] > voxel_d[8*j+2] - interval*1.1 ) ) || ( ( voxel_d[8*i+2] > voxel_d[8*j+2] + interval*0.9 ) && voxel_d[8*i+2] < voxel_d[8*j+2] + interval*1.1 ) ){
                        if (dd_type == THIN_SKIN_CONNECTED_IN_BUNCH)
                            indices_connected[i*skin_bunch + count_near_voxels_thin] = j;
                        count_near_voxels_thin += 1;
                    }
                }
            }
            if (voxel_d[8*i] == voxel_d[8*j]){
                if (voxel_d[8*i+2] == voxel_d[8*j+2]){
                    if ( ( ( voxel_d[8*i+1] < voxel_d[8*j+1] - interval*0.9 ) && ( voxel_d[8*i+1] > voxel_d[8*j+1] - interval*1.1 ) ) || ( ( voxel_d[8*i+1] > voxel_d[8*j+1] + interval*0.9 ) && voxel_d[8*i+1] < voxel_d[8*j+1] + interval*1.1 ) ){
                        if (dd_type == THIN_SKIN_CONNECTED_IN_BUNCH)
                            indices_connected[i*skin_bunch + count_near_voxels_thin] = j;
                        count_near_voxels_thin += 1;
                    }
                }
            }
            if (voxel_d[8*i+1] == voxel_d[8*j+1]){
                if (voxel_d[8*i+2] == voxel_d[8*j+2]){
                    if ( ( ( voxel_d[8*i] < voxel_d[8*j] - interval*0.9 ) && ( voxel_d[8*i] > voxel_d[8*j] - interval*1.1 ) ) || ( ( voxel_d[8*i] > voxel_d[8*j] + interval*0.9 ) && voxel_d[8*i] < voxel_d[8*j] + interval*1.1 ) ){
                        if (dd_type == THIN_SKIN_CONNECTED_IN_BUNCH)
                            indices_connected[i*skin_bunch + count_near_voxels_thin] = j;
                        count_near_voxels_thin += 1;
                    }
                }
            }
        }
    }
/// Thick skin algorithms
    if (di_type == THICK_SKIN || dd_type == THICK_SKIN_CONNECTED_DIRECTLY || dd_type == THICK_SKIN_CONNECTED_IN_BUNCH){
        for (int j=0; j<*c_vox; j++){
            if ( ( voxel_d[8*i+1] > voxel_d[8*j+1] - interval*1.1 ) && ( voxel_d[8*i+1] < voxel_d[8*j+1] + interval*1.1 ) ){
                if ( ( voxel_d[8*i+2] > voxel_d[8*j+2] - interval*1.1 ) && ( voxel_d[8*i+2] < voxel_d[8*j+2] + interval*1.1 ) ){
                    if ( ( voxel_d[8*i] > voxel_d[8*j] - interval*1.1 ) && ( voxel_d[8*i] < voxel_d[8*j] + interval*1.1 ) ){
                        if ( ( voxel_d[8*i] == voxel_d[8*j] ) && ( voxel_d[8*i+1] == voxel_d[8*j+1] ) && ( voxel_d[8*i+2] == voxel_d[8*j+2] ) )
                            continue;
                        if (dd_type == THICK_SKIN_CONNECTED_IN_BUNCH)
                            indices_connected[i*skin_bunch + count_near_voxels_thick] = j;
                        count_near_voxels_thick += 1;
                    }
                }
            }
        }
    }

    if (di_type == THIN_SKIN){
        if (count_near_voxels_thin == 6)     /// setting internal voxels invalid if there are 6 voxels connecting to this voxel (1 for every face of a cube) in thin skin
            voxel_d[8*i+7] = 0.0;
    }
    if (di_type == THICK_SKIN){
        if (count_near_voxels_thick == 26)
            voxel_d[8*i+7] = 0.0;
    }
    if (dd_type == THIN_SKIN_CONNECTED_DIRECTLY){  /// So here comes the fast algorithm. Its criterion is the number of surrounding voxels
        if (count_near_voxels_thin <= del_disc_crit)
            voxel_d[8*i+7] = 0.0;
    }
    if (dd_type == THICK_SKIN_CONNECTED_DIRECTLY){
        if (count_near_voxels_thick <= del_disc_crit)
            voxel_d[8*i+7] = 0.0;
    }
}
)"""";
/// Now inside indices_connected[n][0-5 or 0-26] should be indices of connected verts (for bunch)

string Delete_Disconnected_Bunch_Compute = R""""(
void kernel Kernel_Delete_Disconnected_Bunch (const global int *c_vox, global float *voxel_d, const int del_disc_crit, const global int *indices_connected, const int skin_bunch, global int *voxel_group){

                 /// And this is the slow algorithm. It groups voxels and checks how many are connected in each group

        int voxel_group_size = 1;
        int global_id = get_global_id(0);
        voxel_group[ (*c_vox) * global_id ] = global_id;
        int element = 0;

        for (int j=0; j<skin_bunch; j++){
            if ( indices_connected [ voxel_group[ (*c_vox) * global_id + element ] * skin_bunch + j ] == -1 ){
                element += 1;
                j=0;
                if (element == voxel_group_size)
                    break;
            }
            bool result = false;
            for (int i=0; i<voxel_group_size; i++){
                if ( indices_connected [ voxel_group[(*c_vox) * global_id + element] * skin_bunch + j ] == voxel_group[(*c_vox) * global_id + i] )
                    result = true;
            }
            if (result == false){
                voxel_group[ (*c_vox) * global_id + voxel_group_size] = indices_connected [ voxel_group[(*c_vox) * global_id + element] * skin_bunch + j ];
                voxel_group_size += 1;
            }
        }
        if (voxel_group_size <= del_disc_crit){
            for (int j=0; j<voxel_group_size; j++){
                voxel_d[ voxel_group[(*c_vox) * global_id + j] * 8 + 7 ] = 0.0;       /// mark voxels for deletion  /// change this to atomic!!!
            }
        }
}
)"""";
