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
void kernel Kernel_Pow9(const int iterations, const int resolution, const float interval, const float min, const float n_order,
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
        for( int i = 1; i <= iterations; i++ ){
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
string Consol_Compute = R""""(
void kernel Kernel_Consol (global int *c_vox, global float *voxel_d){

    if ( get_local_id(0) == 0){
        if ( get_local_id(1) == 0){
            if ( get_local_id(2) == 0){
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
    }
}
)"""";
string Faces_Vertices_FB = R""""(
void kernel Faces_Vertices_FB (const float half_vox_size, global float *voxel_d, global unsigned int *faces, global float *verts) {

    unsigned int n = get_global_id(0);
    int x = voxel_d[7*n];
    int y = voxel_d[7*n+1];
    int z = voxel_d[7*n+2];
    verts[24*n] = x-half_vox_size;
    verts[24*n+1] = y-half_vox_size;
    verts[24*n+2] = z-half_vox_size;
    verts[24*n+3] = x+half_vox_size;
    verts[24*n+4] = y-half_vox_size;
    verts[24*n+5] = z-half_vox_size;
    verts[24*n+6] = x+half_vox_size;
    verts[24*n+7] = y+half_vox_size;
    verts[24*n+8] = z-half_vox_size;
    verts[24*n+9] = x-half_vox_size;
    verts[24*n+10] = y+half_vox_size;
    verts[24*n+11] = z-half_vox_size;
    verts[24*n+12] = x-half_vox_size;
    verts[24*n+13] = y-half_vox_size;
    verts[24*n+14] = z+half_vox_size;
    verts[24*n+15] = x+half_vox_size;
    verts[24*n+16] = y-half_vox_size;
    verts[24*n+17] = z+half_vox_size;
    verts[24*n+18] = x+half_vox_size;
    verts[24*n+19] = y+half_vox_size;
    verts[24*n+20] = z+half_vox_size;
    verts[24*n+21] = x-half_vox_size;
    verts[24*n+22] = y+half_vox_size;
    verts[24*n+23] = z+half_vox_size;
    faces[24*n] = n*8+3;
    faces[24*n+1] = n*8+2;
    faces[24*n+2] = n*8+1;
    faces[24*n+3] = n*8;
    faces[24*n+4] = n*8+4;
    faces[24*n+5] = n*8+5;
    faces[24*n+6] = n*8+6;
    faces[24*n+7] = n*8+7;
    faces[24*n+8] = n*8;
    faces[24*n+9] = n*8+1;
    faces[24*n+10] = n*8+5;
    faces[24*n+11] = n*8+4;
    faces[24*n+12] = n*8+1;
    faces[24*n+13] = n*8+2;
    faces[24*n+14] = n*8+6;
    faces[24*n+15] = n*8+5;
    faces[24*n+16] = n*8+2;
    faces[24*n+17] = n*8+3;
    faces[24*n+18] = n*8+7;
    faces[24*n+19] = n*8+6;
    faces[24*n+20] = n*8;
    faces[24*n+21] = n*8+4;
    faces[24*n+22] = n*8+7;
    faces[24*n+23] = n*8+3;
}
)"""";
