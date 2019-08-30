#include "operators.h"
#include "split_op.h"
#include "fileIO.h"
#include "kernels.h"
#include "dynamic.h"

void laplacian(Grid &par, double2 *data, double2* out, int xDim, int yDim,
               int zDim, double dx, double dy, double dz){

    dim3 grid = par.grid;
    dim3 threads = par.threads;
    int gsize = xDim * yDim * zDim;

    double2 *temp_derivative;
    cudaHandleError( cudaMalloc((void **) &temp_derivative, sizeof(double2)*gsize) );
    derive<<<grid, threads>>>(data, temp_derivative, 1, gsize, dx);
    cudaCheckError();
    derive<<<grid, threads>>>(temp_derivative, temp_derivative, 1, gsize, dx);
    cudaCheckError();

    copy<<<grid, threads>>>(temp_derivative, out);
    cudaCheckError();

    derive<<<grid, threads>>>(data, temp_derivative, xDim, gsize, dy);
    cudaCheckError();
    derive<<<grid, threads>>>(temp_derivative, temp_derivative,
                              xDim, gsize, dy);
    cudaCheckError();

    sum<<<grid, threads>>>(temp_derivative, out, out);
    cudaCheckError();

    derive<<<grid, threads>>>(data, temp_derivative, xDim*yDim, gsize, dz);
    cudaCheckError();
    derive<<<grid, threads>>>(temp_derivative, temp_derivative,
                              xDim*yDim, gsize, dz);
    cudaCheckError();

    sum<<<grid, threads>>>(temp_derivative, out, out);
    cudaCheckError();

    cudaHandleError( cudaFree(temp_derivative) );

}

void laplacian(Grid &par, double2 *data, double2* out, int xDim, int yDim,
               double dx, double dy){


    dim3 grid = par.grid;
    dim3 threads = par.threads;
    int gsize = xDim * yDim;

    double2 *temp_derivative;
    cudaHandleError( cudaMalloc((void **) &temp_derivative, sizeof(double2)*gsize) );
    derive<<<grid, threads>>>(data, temp_derivative, 1, gsize, dx);
    cudaCheckError();
    derive<<<grid, threads>>>(temp_derivative, temp_derivative, 1, gsize, dx);
    cudaCheckError();

    copy<<<grid, threads>>>(temp_derivative, out);
    cudaCheckError();

    derive<<<grid, threads>>>(data, temp_derivative, xDim, gsize, dy);
    cudaCheckError();
    derive<<<grid, threads>>>(temp_derivative, temp_derivative,
                              xDim, gsize, dy);
    cudaCheckError();

    sum<<<grid, threads>>>(temp_derivative, out, out);
    cudaCheckError();

    cudaHandleError( cudaFree(temp_derivative) );
}

void laplacian(Grid &par, double2 *data, double2* out, int xDim, double dx){

    dim3 grid = par.grid;
    dim3 threads = par.threads;
    int gsize = xDim;

    derive<<<grid, threads>>>(data, out, 1, gsize, dx);
    cudaCheckError();
    derive<<<grid, threads>>>(out, out, 1, gsize, dx);
    cudaCheckError();

}

double sign(double x){
    if (x < 0){
        return -1.0;
    }
    else if (x == 0){
        return 0.0;
    }
    else{
        return 1.0;
    }
}

// Function to take the curl of Ax and Ay in 2d
// note: This is on the cpu, there should be a GPU version too.
double *curl2d(Grid &par, double *Ax, double *Ay){
    int xDim = par.ival("xDim");
    int yDim = par.ival("yDim");

    int size = sizeof(double) * xDim * yDim;
    double *curl;
    curl = (double *)calloc(size, sizeof(double));

    int index;

    // Note: To take the curl, we need a change in x and y to create a dx or dy
    //       For this reason, we have added yDim to y and 1 to x
    for (int i = 0; i < yDim-1; i++){
        for (int j = 0; j < xDim-1; j++){
            index = j + xDim * i;
            curl[index] = (Ay[index] - Ay[index+1]) 
                           - (Ax[index] - Ax[index+yDim]);
        }
    }

    return curl;
}

double *curl3d_r(Grid &par, double *Bx, double *By){
    int xDim = par.ival("xDim");
    int yDim = par.ival("yDim");
    int zDim = par.ival("zDim");

    int size = sizeof(double) * xDim * yDim * zDim;
    double *curl;
    curl = (double *)malloc(size);

    for (int i = 0; i < xDim*yDim*zDim; ++i){
        curl[i] = sqrt(Bx[i]*Bx[i] + By[i] * By[i]);
    }

    return curl;
}

double *curl3d_phi(Grid &par, double *Bx, double *By){
    int xDim = par.ival("xDim");
    int yDim = par.ival("yDim");
    int zDim = par.ival("zDim");

    int size = sizeof(double) * xDim * yDim * zDim;
    double *curl;
    curl = (double *)malloc(size);

    for (int i = 0; i < xDim*yDim*zDim; ++i){
        curl[i] = atan2(By[i], Bx[i])+M_PI;
    }

    return curl;
}


// Function to take the curl of Ax and Ay in 2d
// note: This is on the cpu, there should be a GPU version too.
// Not complete yet!
double *curl3d_x(Grid &par, double *Ax, double *Ay, double *Az){
    int xDim = par.ival("xDim");
    int yDim = par.ival("yDim");
    int zDim = par.ival("zDim");

    int size = sizeof(double) * xDim * yDim * zDim;
    double *curl;
    curl = (double *)calloc(size, sizeof(double));

    int index;

    // Note: To take the curl, we need a change in x and y to create a dx or dy
    //       For this reason, we have added yDim to y and 1 to x
    for (int i = 0; i < zDim-1; i++){
        for (int j = 0; j < yDim-1; j++){
            for (int k = 0; k < xDim-1; k++){
                index = k + xDim * j + xDim * yDim * i;
                curl[index] = (Az[index] - Az[index + zDim])
                              -(Ay[index] - Ay[index + zDim*yDim]);
            }
        }
    }

    return curl;
}

// Function to take the curl of Ax and Ay in 2d
// note: This is on the cpu, there should be a GPU version too.
// Not complete yet!
double *curl3d_y(Grid &par, double *Ax, double *Ay, double *Az){
    int xDim = par.ival("xDim");
    int yDim = par.ival("yDim");
    int zDim = par.ival("zDim");

    int size = sizeof(double) * xDim * yDim * zDim;
    double *curl;
    curl = (double *)calloc(size, sizeof(double));

    int index;

    // Note: To take the curl, we need a change in x and y to create a dx or dy
    //       For this reason, we have added yDim to y and 1 to x
    for (int i = 0; i < zDim-1; i++){
        for (int j = 0; j < yDim-1; j++){
            for (int k = 0; k < xDim - 1; k++){
                index = k + xDim * j + xDim * yDim * i;
                curl[index] = -(Az[index] - Az[index + 1])
                              -(Ax[index] - Ax[index + xDim*yDim]);
            }
        }
    }

    return curl;
}

// Function to take the curl of Ax and Ay in 2d
// note: This is on the cpu, there should be a GPU version too.
// Not complete yet!
double *curl3d_z(Grid &par, double *Ax, double *Ay, double *Az){
    int xDim = par.ival("xDim");
    int yDim = par.ival("yDim");
    int zDim = par.ival("zDim");

    int size = sizeof(double) * xDim * yDim * zDim;
    double *curl;
    curl = (double *)calloc(size, sizeof(double));

    int index;

    // Note: To take the curl, we need a change in x and y to create a dx or dy
    //       For this reason, we have added yDim to y and 1 to x
    for (int i = 0; i < zDim-1; i++){
        for (int j = 0; j < yDim-1; j++){
            for (int k = 0; k < xDim-1; k++){
                index = k + xDim * j + xDim * yDim * i;
                curl[index] = (Ay[index] - Ay[index + 1])
                              -(Ax[index] - Ax[index + xDim]);
            }
        }
    }

    return curl;
}

// Function to check whether a file exists
std::string filecheck(std::string filename){

    struct stat buffer = {0};
    if (stat(filename.c_str(), &buffer) == -1){
        std::cout << "File " << filename << " does not exist!" << '\n';
        std::cout << "Please select a new file:" << '\n'; 
        std::cin >> filename; 
        filename = filecheck(filename);
    } 

    return filename;
}

/*----------------------------------------------------------------------------//
* GPU KERNELS
*-----------------------------------------------------------------------------*/

// Function to generate momentum grids
void generate_p_space(Grid &par){

    int dimnum = par.ival("dimnum");
    int xDim = par.ival("xDim");
    int yDim = par.ival("yDim");
    int zDim = par.ival("zDim");
    double xMax = par.dval("xMax");
    double yMax = 0;
    if (dimnum > 1){
        yMax = par.dval("yMax");
    }
    double zMax = 0;
    if (dimnum == 3){
        zMax = par.dval("zMax");
    }

    double pxMax = par.dval("pxMax");
    double pyMax = 0;
    if (dimnum > 1){
        pyMax = par.dval("pyMax");
    }
    double pzMax = 0;
    if (dimnum == 3){
        pzMax = par.dval("pzMax");
    }

    double dx = par.dval("dx");
    double dy = 0;
    if (dimnum > 1){
        dy = par.dval("dy");
    }
    double dz = 0;
    if (dimnum == 3){
        dz = par.dval("dz");
    }

    double dpx = par.dval("dpx");
    double dpy = 0;
    if (dimnum > 1){
        dpy = par.dval("dpy");
    }
    double dpz = 0;
    if (dimnum == 3){
        dpz = par.dval("dpz");
    }

    double *x, *y, *z, *px, *py, *pz,
           *x_gpu, *y_gpu, *z_gpu, 
           *px_gpu, *py_gpu, *pz_gpu;

    x = (double *) malloc(sizeof(double) * xDim);
    y = (double *) malloc(sizeof(double) * yDim);
    z = (double *) malloc(sizeof(double) * zDim);
    px = (double *) malloc(sizeof(double) * xDim);
    py = (double *) malloc(sizeof(double) * yDim);
    pz = (double *) malloc(sizeof(double) * zDim);

    if (dimnum == 2){

        for(int i=0; i<xDim/2; ++i){
            x[i] = -xMax + i*dx;
            x[i + (xDim/2)] = i*dx;

            px[i] = i*dpx;
            px[i + (xDim/2)] = -pxMax + i*dpx;

        }
        for(int i=0; i<yDim/2; ++i){
            y[i] = -yMax + i*dy;
            y[i + (yDim/2)] = i*dy;

            py[i] = i*dpy;
            py[i + (yDim/2)] = -pyMax + i*dpy;

        }

        for(int i = 0; i < zDim; ++i){
            z[i] = 0;
            pz[i] = 0;
        }

    }
    else if(dimnum == 3){
        for(int i=0; i<xDim/2; ++i){
            x[i] = -xMax + i*dx;
            x[i + (xDim/2)] = i*dx;

            px[i] = i*dpx;
            px[i + (xDim/2)] = -pxMax + i*dpx;

        }
        for(int i=0; i<yDim/2; ++i){
            y[i] = -yMax + i*dy;
            y[i + (yDim/2)] = i*dy;

            py[i] = i*dpy;
            py[i + (yDim/2)] = -pyMax + i*dpy;

        }
        for(int i=0; i<zDim/2; ++i){
            z[i] = -zMax + i*dz;
            z[i + (zDim/2)] = i*dz;

            pz[i] = i*dpz;
            pz[i + (zDim/2)] = -pzMax + i*dpz;

        }

    }
    else if (dimnum == 1){
        for(int i=0; i<xDim/2; ++i){
            x[i] = -xMax + i*dx;
            x[i + (xDim/2)] = i*dx;

            px[i] = i*dpx;
            px[i + (xDim/2)] = -pxMax + i*dpx;

        }

        for(int i = 0; i < zDim; ++i){
            z[i] = 0;
            pz[i] = 0;
            y[i] = 0;
            py[i] = 0;
        }

    }
    par.store("x",x);
    par.store("y",y);
    par.store("z",z);
    par.store("px",px);
    par.store("py",py);
    par.store("pz",pz);

    // Now move these items to the gpu
    cudaHandleError( cudaMalloc((void**) &x_gpu, sizeof(double) * xDim) );
    cudaHandleError( cudaMalloc((void**) &y_gpu, sizeof(double) * yDim) );
    cudaHandleError( cudaMalloc((void**) &z_gpu, sizeof(double) * zDim) );
    cudaHandleError( cudaMalloc((void**) &px_gpu, sizeof(double) * xDim) );
    cudaHandleError( cudaMalloc((void**) &py_gpu, sizeof(double) * yDim) );
    cudaHandleError( cudaMalloc((void**) &pz_gpu, sizeof(double) * zDim) );

    cudaHandleError( cudaMemcpy(x_gpu, x, sizeof(double)*xDim, cudaMemcpyHostToDevice) );
    cudaHandleError( cudaMemcpy(y_gpu, y, sizeof(double)*yDim, cudaMemcpyHostToDevice) );
    cudaHandleError( cudaMemcpy(z_gpu, z, sizeof(double)*zDim, cudaMemcpyHostToDevice) );
    cudaHandleError( cudaMemcpy(px_gpu, px, sizeof(double)*xDim, cudaMemcpyHostToDevice) );
    cudaHandleError( cudaMemcpy(py_gpu, py, sizeof(double)*yDim, cudaMemcpyHostToDevice) );
    cudaHandleError( cudaMemcpy(pz_gpu, pz, sizeof(double)*zDim, cudaMemcpyHostToDevice) );

    par.store("x_gpu",x_gpu);
    par.store("y_gpu",y_gpu);
    par.store("z_gpu",z_gpu);
    par.store("px_gpu",px_gpu);
    par.store("py_gpu",py_gpu);
    par.store("pz_gpu",pz_gpu);
}

// This function is basically a wrapper to call the appropriate K kernel
void generate_K(Grid &par){

    // For k, we need xp, yp, and zp. These will also be used in generating 
    // pAxyz parameters, so it should already be stored in par.
    double *px_gpu = par.dsval("px_gpu");
    double *py_gpu = par.dsval("py_gpu");
    double *pz_gpu = par.dsval("pz_gpu");
    double gSize = par.ival("gSize");
    double mass = par.dval("mass");
    int wfc_num = par.ival("wfc_num");

    // Creating K to work with
    std::vector<double *> K(wfc_num), K_gpu(wfc_num);
    for (int w = 0; w < wfc_num; ++w){
        K[w] = (double*)malloc(sizeof(double)*gSize);
        cudaHandleError(cudaMalloc((void**) &K_gpu[w], sizeof(double)*gSize));

        simple_K<<<par.grid, par.threads>>>(px_gpu, py_gpu, pz_gpu,
                                            mass, K_gpu[w]);

        cudaHandleError(cudaMemcpy(K[w], K_gpu[w], sizeof(double)*gSize,
                        cudaMemcpyDeviceToHost));
    }
    par.store("K",K);
    par.store("K_gpu",K_gpu);
}

// Simple kernel for generating K
__global__ void simple_K(double *xp, double *yp, double *zp, double mass,
                         double *K){

    unsigned int gid = getGid3d3d();
    unsigned int xid = blockDim.x*blockIdx.x + threadIdx.x;
    unsigned int yid = blockDim.y*blockIdx.y + threadIdx.y;
    unsigned int zid = blockDim.z*blockIdx.z + threadIdx.z;
    K[gid] = (HBAR*HBAR/(2*mass))*(xp[xid]*xp[xid] + yp[yid]*yp[yid]
                                  + zp[zid]*zp[zid]);
}

// Function to generate game fields
void generate_gauge(Grid &par){
    int gSize = par.ival("gSize");
    int dimnum = par.ival("dimnum");
    int wfc_num = par.ival("wfc_num");

    std::vector<double *> Ax(wfc_num), Ay(wfc_num), Az(wfc_num),
                          Ax_gpu(wfc_num), Ay_gpu(wfc_num), Az_gpu(wfc_num);
    double *x_gpu = par.dsval("x_gpu");
    double *y_gpu = par.dsval("y_gpu");
    double *z_gpu = par.dsval("z_gpu");

    double xMax = par.dval("xMax");
    double yMax = par.dval("yMax");
    double zMax = 1;
    if (dimnum == 3){
        zMax = par.dval("zMax");
    }
    double omegaX = par.dval("omegaX");
    double omegaY = par.dval("omegaY");
    double omegaZ;
    if (dimnum == 3){
        omegaZ = par.dval("omegaZ");
    }
    double omega = par.dval("omega");
    double fudge = par.dval("fudge");
    for (int w = 0; w < wfc_num;++w){

        Ax[w] = (double *)malloc(sizeof(double)*gSize);
        Ay[w] = (double *)malloc(sizeof(double)*gSize);
        Az[w] = (double *)malloc(sizeof(double)*gSize);

        cudaHandleError(cudaMalloc((void**) &Ax_gpu[w], sizeof(double)*gSize));
        cudaHandleError(cudaMalloc((void**) &Ay_gpu[w], sizeof(double)*gSize));
        cudaHandleError(cudaMalloc((void**) &Az_gpu[w], sizeof(double)*gSize));
        if (par.is_ast_gpu("Ax")){
            double dx = par.dval("dx");
            double dy = par.dval("dy");
            double dz = par.dval("dz");
            double xMax = par.dval("xMax");
            double yMax = par.dval("yMax");
            double zMax = 0;
            if (dimnum == 3){ 
                zMax = par.dval("zMax");
            }

            EqnNode_gpu *eqn = par.astval("Ax");

            find_field<<<par.grid, par.threads>>>(Ax_gpu[w], dx, dy, dz, 
                                                  xMax, yMax, zMax, 0, eqn);
            cudaCheckError();
        }
        else{
            par.Ax_fn<<<par.grid, par.threads>>>(x_gpu, y_gpu, z_gpu, 
                                                  xMax, yMax, zMax, 
                                                  omegaX, omegaY, omegaZ, 
                                                  omega, fudge, Ax_gpu[w]);
            cudaCheckError();
        }
        if (par.is_ast_gpu("Ay")){
            double dx = par.dval("dx");
            double dy = par.dval("dy");
            double dz = par.dval("dz");
            double xMax = par.dval("xMax");
            double yMax = par.dval("yMax");
            double zMax = 0;
            if (dimnum == 3){
                zMax = par.dval("zMax");
            }

            EqnNode_gpu *eqn = par.astval("Ay");

            find_field<<<par.grid, par.threads>>>(Ay_gpu[w], dx, dy, dz,
                                                  xMax, yMax, zMax, 0, eqn);
            cudaCheckError();
        }
        else{
            par.Ay_fn<<<par.grid, par.threads>>>(x_gpu, y_gpu, z_gpu, 
                                                  xMax, yMax, zMax, 
                                                  omegaX, omegaY, omegaZ, 
                                                  omega, fudge, Ay_gpu[w]);
            cudaCheckError();
        }
        if (dimnum == 3){
            if (par.is_ast_gpu("Az")){
                double dx = par.dval("dx");
                double dy = par.dval("dy");
                double dz = par.dval("dz");

                double xMax = par.dval("xMax");
                double yMax = par.dval("yMax");
                double zMax = 0;
                if (dimnum == 3){
                    zMax = par.dval("zMax");
                }

                EqnNode_gpu *eqn = par.astval("Az");

                find_field<<<par.grid, par.threads>>>(Az_gpu[w], dx, dy, dz,
                                                      xMax, yMax, zMax,  
                                                      0, eqn);
                cudaCheckError();
            }
            else{
                par.Az_fn<<<par.grid, par.threads>>>(x_gpu, y_gpu, z_gpu, 
                                                      xMax, yMax, zMax, 
                                                      omegaX, omegaY,
                                                      omegaZ, omega, fudge,
                                                      Az_gpu[w]);
                cudaCheckError();
            }
        }
        else{
            kconstant_A<<<par.grid, par.threads>>>(x_gpu, y_gpu, z_gpu, 
                                                    xMax, yMax, zMax, 
                                                    omegaX, omegaY, omegaZ, 
                                                    omega, fudge, Az_gpu[w]);
            cudaCheckError();
        }
        cudaHandleError(cudaMemcpy(Ax[w], Ax_gpu[w],
                   sizeof(double)*gSize,cudaMemcpyDeviceToHost));
        cudaHandleError(cudaMemcpy(Ay[w], Ay_gpu[w],
                   sizeof(double)*gSize,cudaMemcpyDeviceToHost));
        cudaHandleError(cudaMemcpy(Az[w], Az_gpu[w],
                   sizeof(double)*gSize,cudaMemcpyDeviceToHost));
    }

    par.store("Ax", Ax);
    par.store("Ay", Ay);
    par.store("Az", Az);

    par.store("Ax_gpu", Ax_gpu);
    par.store("Ay_gpu", Ay_gpu);
    par.store("Az_gpu", Az_gpu);

}

// constant Kernel A
__global__ void kconstant_A(double *x, double *y, double *z,
                            double xMax, double yMax, double zMax,
                            double omegaX, double omegaY, double omegaZ,
                            double omega, double fudge, double *A){
    int gid = getGid3d3d();
    A[gid] = 0;        
}

// Kernel for simple rotational case, Ax
__global__ void krotation_Ax(double *x, double *y, double *z,
                             double xMax, double yMax, double zMax,
                             double omegaX, double omegaY, double omegaZ,
                             double omega, double fudge, double *A){
    int gid = getGid3d3d();
    int yid = blockDim.y*blockIdx.y + threadIdx.y;
    A[gid] = -y[yid] * omega * omegaX;
}

// Kernel for simple rotational case, Ay
__global__ void krotation_Ay(double *x, double *y, double *z,
                             double xMax, double yMax, double zMax,
                             double omegaX, double omegaY, double omegaZ,
                             double omega, double fudge, double *A){
    int gid = getGid3d3d();
    int xid = blockDim.x*blockIdx.x + threadIdx.x;
    A[gid] = x[xid] * omega * omegaY;
}

// Kernel for simple rotational case, Ax
__global__ void kring_rotation_Ax(double *x, double *y, double *z,
                                  double xMax, double yMax, double zMax,
                                  double omegaX, double omegaY, double omegaZ,
                                  double omega, double fudge, double *A){
    int gid = getGid3d3d();
    int xid = blockDim.x*blockIdx.x + threadIdx.x;
    int yid = blockDim.y*blockIdx.y + threadIdx.y;
    int zid = blockDim.z*blockIdx.z + threadIdx.z;
    double theta = atan2(y[yid],x[xid]);
    A[gid] = (z[zid]+zMax)*cos(theta)*omega*omegaX;
}

// Kernel for simple rotational case, Ay
__global__ void kring_rotation_Ay(double *x, double *y, double *z,
                                  double xMax, double yMax, double zMax,
                                  double omegaX, double omegaY, double omegaZ,
                                  double omega, double fudge, double *A){
    int gid = getGid3d3d();
    int xid = blockDim.x*blockIdx.x + threadIdx.x;
    int yid = blockDim.y*blockIdx.y + threadIdx.y;
    int zid = blockDim.z*blockIdx.z + threadIdx.z;
    double theta = atan2(y[yid],x[xid]);
    A[gid] = (z[zid]+zMax)*sin(theta)*omega*omegaX;
}

// Kernel for simple rotational case, Az
__global__ void kring_rotation_Az(double *x, double *y, double *z,
                                  double xMax, double yMax, double zMax,
                                  double omegaX, double omegaY, double omegaZ,
                                  double omega, double fudge, double *A){
    int gid = getGid3d3d();
    int xid = blockDim.x*blockIdx.x + threadIdx.x;
    int yid = blockDim.y*blockIdx.y + threadIdx.y;
    double r = sqrt(x[xid]*x[xid] + y[yid]*y[yid]);
    A[gid] = r*omega*omegaX;
}


// kernel for a simple vortex ring
__global__ void kring_Az(double *x, double *y, double *z,
                         double xMax, double yMax, double zMax,
                         double omegaX, double omegaY, double omegaZ,
                         double omega, double fudge, double *A){
    int gid = getGid3d3d();
    int xid = blockDim.x*blockIdx.x + threadIdx.x;
    int yid = blockDim.y*blockIdx.y + threadIdx.y;

    double rad = sqrt(x[xid]*x[xid] + y[yid]*y[yid]);

    A[gid] = omega * exp(-rad*rad / (0.0001*xMax)) * 0.01;
}

// testing kernel Ax
__global__ void ktest_Ax(double *x, double *y, double *z,
                         double xMax, double yMax, double zMax,
                         double omegaX, double omegaY, double omegaZ,
                         double omega, double fudge, double *A){
    int gid = getGid3d3d();
    int yid = blockDim.y*blockIdx.y + threadIdx.y;
    A[gid] = (sin(y[yid] * 100000)+1) * yMax * omega;
}

// testing kernel Ay
__global__ void ktest_Ay(double *x, double *y, double *z,
                         double xMax, double yMax, double zMax,
                         double omegaX, double omegaY, double omegaZ,
                         double omega, double fudge, double *A){
    int gid = getGid3d3d();
    A[gid] = 0;
}

// function to generate V
void generate_fields(Grid &par){

    generate_p_space(par);
    generate_K(par);
    std::cout << "generating gauge fields...\n";
    if (par.bval("read_a")) {
        FileIO::loadA(par);
    } else {
        generate_gauge(par);
    }

    int gSize = par.ival("gSize");
    int dimnum = par.ival("dimnum");
    int winding = par.dval("winding");
    int wfc_num = par.ival("wfc_num");

    bool energy_calc = par.bval("energy_calc");

    double dt = par.dval("dt");
    double gdt = par.dval("gdt");
    double *x_gpu = par.dsval("x_gpu");
    double *y_gpu = par.dsval("y_gpu");
    double *z_gpu = par.dsval("z_gpu");
    double *px_gpu = par.dsval("px_gpu");
    double *py_gpu = par.dsval("py_gpu");
    double *pz_gpu = par.dsval("pz_gpu");
    std::vector<double *> Ax_gpu = par.dsvecval("Ax_gpu");
    std::vector<double *> Ay_gpu = par.dsvecval("Ay_gpu");
    std::vector<double *> Az_gpu = par.dsvecval("Az_gpu");
    std::vector<double *> K_gpu = par.dsvecval("K_gpu");

    // Creating items list for kernels

    double *items, *items_gpu;
    int item_size = 18;
    items = (double*)malloc(sizeof(double)*item_size);
    cudaHandleError( cudaMalloc((void**) &items_gpu, sizeof(double)*item_size) );

    for (int i = 0; i < item_size; ++i){
        items[i] = 0;
    }
    items[0] = par.dval("xMax");
    items[1] = par.dval("yMax");
    if (dimnum == 3){
        items[2] = par.dval("zMax");
    }

    items[3] = par.dval("omegaX");
    items[4] = par.dval("omegaY");
    if (dimnum == 3){
        items[5] = par.dval("omegaZ");
    }

    items[6] = par.dval("x0_shift");
    items[7] = par.dval("y0_shift");
    if (dimnum == 3){
        items[8] = par.dval("z0_shift");
    }
    else{
        items[8] = 0.0;
    }

    items[9] = par.dval("mass");
    items[10] = par.dval("gammaY");
    items[11] = 1.0; // For gammaZ
    items[12] = par.dval("fudge");
    items[13] = 0.0; // For time

    items[14] = par.dval("Rxy");

    items[15] = par.dval("a0x");
    items[16] = par.dval("a0y");
    if (dimnum == 3){
        items[17] = par.dval("a0z");
    }
    else{
        items[17] = 1.0;
    }

    cudaHandleError( cudaMemcpy(items_gpu, items, sizeof(double)*item_size,
                                cudaMemcpyHostToDevice) );

    double fudge = par.dval("fudge");

    // Generating V

    std::vector<double *> V(wfc_num), V_gpu(wfc_num);

    // Generating wfc
    std::vector<double2 *> wfc_array(wfc_num), wfc_gpu_array(wfc_num);
    std::vector<double *> phi(wfc_num), phi_gpu(wfc_num);

    // generating aux fields.
    std::vector<double2 *> GV(wfc_num), EV(wfc_num), GK(wfc_num), EK(wfc_num);
    std::vector<double2 *> GV_gpu(wfc_num), EV_gpu(wfc_num), GK_gpu(wfc_num),
                           EK_gpu(wfc_num);
    std::vector<double2 *> GpAx(wfc_num), GpAy(wfc_num), GpAz(wfc_num),
                           EpAx(wfc_num), EpAy(wfc_num), EpAz(wfc_num);
    std::vector<double2 *> GpAx_gpu(wfc_num), GpAy_gpu(wfc_num), 
                           GpAz_gpu(wfc_num), EpAx_gpu(wfc_num),
                           EpAy_gpu(wfc_num), EpAz_gpu(wfc_num);
    std::vector<double *> pAx(wfc_num), pAy(wfc_num), pAz(wfc_num);
    std::vector<double *> pAx_gpu(wfc_num), pAy_gpu(wfc_num), pAz_gpu(wfc_num);
    
    std::cout << "iterating through all wavefunctions...\n";
    for (int w = 0; w < wfc_num; ++w){

        V[w] = (double *)malloc(sizeof(double)*gSize);

        cudaHandleError(cudaMalloc((void **) &V_gpu[w], sizeof(double)*gSize));

        if (par.is_ast_gpu("V")){
            double dx = par.dval("dx");
            double dy = par.dval("dy");
            double dz = par.dval("dz");

            double xMax = par.dval("xMax");
            double yMax = par.dval("yMax");
            double zMax = 0;
            if (dimnum == 3){ 
                zMax = par.dval("zMax");
            }

            EqnNode_gpu *eqn = par.astval("V");
            find_field<<<par.grid, par.threads>>>(V_gpu[w], dx, dy, dz, 
                                                  xMax, yMax, zMax, 0, eqn);
            cudaCheckError();
        }
        else{
            par.V_fn<<<par.grid, par.threads>>>(x_gpu, y_gpu, z_gpu, items_gpu,
                                                Ax_gpu[w], Ay_gpu[w],
                                                Az_gpu[w], V_gpu[w]);
            cudaCheckError();
        }

        cudaHandleError(cudaMemcpy(V[w], V_gpu[w], sizeof(double)*gSize,
                        cudaMemcpyDeviceToHost));

        wfc_array[w] = (double2 *)malloc(sizeof(double2)*gSize);
        phi[w] = (double *)malloc(sizeof(double)*gSize);

        cudaHandleError(cudaMalloc((void**) &wfc_gpu_array[w],
                        sizeof(double2)*gSize));
        cudaHandleError(cudaMalloc((void**) &phi_gpu[w], sizeof(double)*gSize));

        if (par.bval("read_wfc")){
            wfc_array = par.d2svecval("wfc_array");
            cudaHandleError(cudaMemcpy(wfc_gpu_array[w], wfc_array[w],
                            sizeof(double2)*gSize, cudaMemcpyHostToDevice));
        }
        else{
            par.wfc_fn<<<par.grid, par.threads>>>(x_gpu, y_gpu, z_gpu,
                                                  items_gpu, winding,
                                                  phi_gpu[w], wfc_gpu_array[w]);
            cudaHandleError(cudaMemcpy(wfc_array[w], wfc_gpu_array[w],
                            sizeof(double2)*gSize, cudaMemcpyDeviceToHost));
        }
    
        cudaHandleError(cudaMemcpy(phi[w], phi_gpu[w], sizeof(double)*gSize,
                        cudaMemcpyDeviceToHost));

        GV[w] = (double2 *)malloc(sizeof(double2)*gSize);
        EV[w] = (double2 *)malloc(sizeof(double2)*gSize);
        GK[w] = (double2 *)malloc(sizeof(double2)*gSize);
        EK[w] = (double2 *)malloc(sizeof(double2)*gSize);
    
        GpAx[w] = (double2 *)malloc(sizeof(double2)*gSize);
        EpAx[w] = (double2 *)malloc(sizeof(double2)*gSize);
        GpAy[w] = (double2 *)malloc(sizeof(double2)*gSize);
        EpAy[w] = (double2 *)malloc(sizeof(double2)*gSize);
        GpAz[w] = (double2 *)malloc(sizeof(double2)*gSize);
        EpAz[w] = (double2 *)malloc(sizeof(double2)*gSize);
    
        pAx[w] = (double *)malloc(sizeof(double)*gSize);
        pAy[w] = (double *)malloc(sizeof(double)*gSize);
        pAz[w] = (double *)malloc(sizeof(double)*gSize);
    
        cudaHandleError(cudaMalloc((void**) &GV_gpu[w], sizeof(double2)*gSize));
        cudaHandleError(cudaMalloc((void**) &EV_gpu[w], sizeof(double2)*gSize));
        cudaHandleError(cudaMalloc((void**) &GK_gpu[w], sizeof(double2)*gSize));
        cudaHandleError(cudaMalloc((void**) &EK_gpu[w], sizeof(double2)*gSize));
    
        cudaHandleError(cudaMalloc((void**) &GpAx_gpu[w],
                        sizeof(double2)*gSize));
        cudaHandleError(cudaMalloc((void**) &EpAx_gpu[w],
                        sizeof(double2)*gSize));
        cudaHandleError(cudaMalloc((void**) &GpAy_gpu[w],
                        sizeof(double2)*gSize));
        cudaHandleError(cudaMalloc((void**) &EpAy_gpu[w],
                        sizeof(double2)*gSize));
        cudaHandleError(cudaMalloc((void**) &GpAz_gpu[w],
                        sizeof(double2)*gSize));
        cudaHandleError(cudaMalloc((void**) &EpAz_gpu[w],
                        sizeof(double2)*gSize));
    
        cudaHandleError(cudaMalloc((void**) &pAx_gpu[w], sizeof(double)*gSize));
        cudaHandleError(cudaMalloc((void**) &pAy_gpu[w], sizeof(double)*gSize));
        cudaHandleError(cudaMalloc((void**) &pAz_gpu[w], sizeof(double)*gSize));
    
        aux_fields<<<par.grid, par.threads>>>(V_gpu[w], K_gpu[w], gdt, dt,
                                              Ax_gpu[w], Ay_gpu[w], Az_gpu[w],
                                              px_gpu, py_gpu, pz_gpu,
                                              pAx_gpu[w], pAy_gpu[w],
                                              pAz_gpu[w], GV_gpu[w], EV_gpu[w],
                                              GK_gpu[w], EK_gpu[w], GpAx_gpu[w],
                                              GpAy_gpu[w], GpAz_gpu[w],
                                              EpAx_gpu[w], EpAy_gpu[w],
                                              EpAz_gpu[w]);
        cudaCheckError();

        cudaHandleError(cudaMemcpy(GV[w], GV_gpu[w], sizeof(double2)*gSize,
                        cudaMemcpyDeviceToHost));
        cudaHandleError(cudaMemcpy(EV[w], EV_gpu[w], sizeof(double2)*gSize,
                        cudaMemcpyDeviceToHost));
        cudaHandleError(cudaMemcpy(GK[w], GK_gpu[w], sizeof(double2)*gSize,
                        cudaMemcpyDeviceToHost));
        cudaHandleError(cudaMemcpy(EK[w], EK_gpu[w], sizeof(double2)*gSize,
                        cudaMemcpyDeviceToHost));

        cudaHandleError(cudaMemcpy(GpAx[w], GpAx_gpu[w], sizeof(double2)*gSize,
                        cudaMemcpyDeviceToHost));
        cudaHandleError(cudaMemcpy(EpAx[w], EpAx_gpu[w], sizeof(double2)*gSize,
                        cudaMemcpyDeviceToHost));
        cudaHandleError(cudaMemcpy(GpAy[w], GpAy_gpu[w], sizeof(double2)*gSize,
                        cudaMemcpyDeviceToHost));
        cudaHandleError(cudaMemcpy(EpAy[w], EpAy_gpu[w], sizeof(double2)*gSize,
                        cudaMemcpyDeviceToHost));
        cudaHandleError(cudaMemcpy(GpAz[w], GpAz_gpu[w], sizeof(double2)*gSize,
                        cudaMemcpyDeviceToHost));
        cudaHandleError(cudaMemcpy(EpAz[w], EpAz_gpu[w], sizeof(double2)*gSize,
                        cudaMemcpyDeviceToHost));

        cudaHandleError(cudaMemcpy(pAx[w], pAx_gpu[w], sizeof(double)*gSize,
                        cudaMemcpyDeviceToHost));
        cudaHandleError(cudaMemcpy(pAy[w], pAy_gpu[w], sizeof(double)*gSize,
                        cudaMemcpyDeviceToHost));
        cudaHandleError(cudaMemcpy(pAz[w], pAz_gpu[w], sizeof(double)*gSize,
                        cudaMemcpyDeviceToHost));

        // Storing variables
        //cudaHandleError(cudaFree(phi_gpu));
        cudaHandleError(cudaFree(GV_gpu[w]));
        cudaHandleError(cudaFree(EV_gpu[w]));
        cudaHandleError(cudaFree(GK_gpu[w]));
        cudaHandleError(cudaFree(EK_gpu[w]));

        cudaHandleError(cudaFree(pAx_gpu[w]));
        cudaHandleError(cudaFree(pAy_gpu[w]));
        cudaHandleError(cudaFree(pAz_gpu[w]));

        cudaHandleError(cudaFree(GpAx_gpu[w]));
        cudaHandleError(cudaFree(GpAy_gpu[w]));
        cudaHandleError(cudaFree(GpAz_gpu[w]));

        cudaHandleError(cudaFree(EpAx_gpu[w]));
        cudaHandleError(cudaFree(EpAy_gpu[w]));
        cudaHandleError(cudaFree(EpAz_gpu[w]));

        if (!energy_calc){
            cudaHandleError(cudaFree(K_gpu[w]));
            cudaHandleError(cudaFree(V_gpu[w]));
            cudaHandleError(cudaFree(Ax_gpu[w]));
            cudaHandleError(cudaFree(Ay_gpu[w]));
            cudaHandleError(cudaFree(Az_gpu[w]));
        }
        else{
            par.store("V_gpu",V_gpu);
        }
        std::cout << "wavefunction " << w << " set!\n";
    }
    cudaHandleError(cudaFree(items_gpu));
    cudaHandleError(cudaFree(x_gpu));
    cudaHandleError(cudaFree(y_gpu));
    cudaHandleError(cudaFree(z_gpu));

    cudaHandleError(cudaFree(px_gpu));
    cudaHandleError(cudaFree(py_gpu));
    cudaHandleError(cudaFree(pz_gpu));

    par.store("V",V);
    par.store("items", items);
    //par.store("items_gpu", items_gpu);
    par.store("wfc_array", wfc_array);
    par.store("wfc_gpu_array", wfc_gpu_array);
    par.store("Phi", phi);
    par.store("Phi_gpu", phi_gpu);

    par.store("GV",GV);
    par.store("EV",EV);
    par.store("GK",GK);
    par.store("EK",EK);
    //par.store("GV_gpu",GV_gpu);
    //par.store("EV_gpu",EV_gpu);
    //par.store("GK_gpu",GK_gpu);
    //par.store("EK_gpu",EK_gpu);

    par.store("GpAx",GpAx);
    par.store("EpAx",EpAx);
    par.store("GpAy",GpAy);
    par.store("EpAy",EpAy);
    par.store("GpAz",GpAz);
    par.store("EpAz",EpAz);

    par.store("pAx",pAx);
    par.store("pAy",pAy);
    par.store("pAz",pAz);
    //par.store("pAx_gpu",pAx_gpu);
    //par.store("pAy_gpu",pAy_gpu);
    //par.store("pAz_gpu",pAz_gpu);

}

__global__ void kharmonic_V(double *x, double *y, double *z, double* items,
                            double *Ax, double *Ay, double *Az, double *V){

    int gid = getGid3d3d();
    int xid = blockDim.x*blockIdx.x + threadIdx.x;
    int yid = blockDim.y*blockIdx.y + threadIdx.y;
    int zid = blockDim.z*blockIdx.z + threadIdx.z;

    double V_x = items[3]*(x[xid]+items[6]);
    double V_y = items[10]*items[4]*(y[yid]+items[7]);
    double V_z = items[11]*items[5]*(z[zid]+items[8]);

    V[gid] = 0.5*items[9]*((V_x*V_x + V_y*V_y + V_z*V_z)
             + (Ax[gid]*Ax[gid] + Ay[gid]*Ay[gid] + Az[gid]*Az[gid]));
}

// kernel for simple 3d torus trapping potential
__global__ void ktorus_V(double *x, double *y, double *z, double* items,
                         double *Ax, double *Ay, double *Az, double *V){

    int gid = getGid3d3d();
    int xid = blockDim.x*blockIdx.x + threadIdx.x;
    int yid = blockDim.y*blockIdx.y + threadIdx.y;
    int zid = blockDim.z*blockIdx.z + threadIdx.z;

    double rad = sqrt((x[xid] - items[6]) * (x[xid] - items[6])
                      + (y[yid] - items[7]) * (y[yid] - items[7])) 
                      - 0.5*items[0];
    double omegaR = (items[3]*items[3] + items[4]*items[4]);
    double V_tot = (2*items[5]*items[5]*(z[zid] - items[8])*(z[zid] - items[8])
                    + omegaR*(rad*rad + items[12]*rad*z[zid]));
    V[gid] = 0.5*items[9]*(V_tot
                           + Ax[gid]*Ax[gid]
                           + Ay[gid]*Ay[gid]
                           + Az[gid]*Az[gid]);
}

__global__ void kstd_wfc(double *x, double *y, double *z, double *items,
                         double winding, double *phi, double2 *wfc_array){

    int gid = getGid3d3d();
    int xid = blockDim.x*blockIdx.x + threadIdx.x;
    int yid = blockDim.y*blockIdx.y + threadIdx.y;
    int zid = blockDim.z*blockIdx.z + threadIdx.z;

    phi[gid] = -fmod(winding*atan2(y[yid], x[xid]),2*PI);

    wfc_array[gid].x = exp(-(x[xid]*x[xid]/(items[14]*items[14]*items[15]*items[15]) 
                     + y[yid]*y[yid]/(items[14]*items[14]*items[16]*items[16]) 
                     + z[zid]*z[zid]/(items[14]*items[14]*items[17]*items[17])))
                     * cos(phi[gid]);
    wfc_array[gid].y = -exp(-(x[xid]*x[xid]/(items[14]*items[14]*items[15]*items[15]) 
                     + y[yid]*y[yid]/(items[14]*items[14]*items[16]*items[16]) 
                     + z[zid]*z[zid]/(items[14]*items[14]*items[17]*items[17])))
                     * sin(phi[gid]);

}

__global__ void ktorus_wfc(double *x, double *y, double *z, double *items,
                           double winding, double *phi, double2 *wfc_array){

    int gid = getGid3d3d();
    int xid = blockDim.x*blockIdx.x + threadIdx.x;
    int yid = blockDim.y*blockIdx.y + threadIdx.y;
    int zid = blockDim.z*blockIdx.z + threadIdx.z;

    double rad = sqrt((x[xid] - items[6]) * (x[xid] - items[6])
                      + (y[yid] - items[7]) * (y[yid] - items[7])) 
                      - 0.5*items[0];

    wfc_array[gid].x = exp(-( pow((rad)/(items[14]*items[15]*0.5),2) +
                   pow((z[zid])/(items[14]*items[17]*0.5),2) ) );
    wfc_array[gid].y = 0.0;
}

__global__ void aux_fields(double *V, double *K, double gdt, double dt,
                           double* Ax, double *Ay, double* Az,
                           double *px, double *py, double *pz,
                           double* pAx, double* pAy, double* pAz,
                           double2* GV, double2* EV, double2* GK, double2* EK,
                           double2* GpAx, double2* GpAy, double2* GpAz,
                           double2* EpAx, double2* EpAy, double2* EpAz){
    int gid = getGid3d3d();
    int xid = blockDim.x*blockIdx.x + threadIdx.x;
    int yid = blockDim.y*blockIdx.y + threadIdx.y;
    int zid = blockDim.z*blockIdx.z + threadIdx.z;

    GV[gid].x = exp(-V[gid]*(gdt/(2*HBAR)));
    GK[gid].x = exp(-K[gid]*(gdt/HBAR));
    GV[gid].y = 0.0;
    GK[gid].y = 0.0;

    // Ax and Ay will be calculated here but are used only for
    // debugging. They may be needed later for magnetic field calc

    pAx[gid] = Ax[gid] * px[xid];
    pAy[gid] = Ay[gid] * py[yid];
    pAz[gid] = Az[gid] * pz[zid];

    GpAx[gid].x = exp(-pAx[gid]*gdt);
    GpAx[gid].y = 0;
    GpAy[gid].x = exp(-pAy[gid]*gdt);
    GpAy[gid].y = 0;
    GpAz[gid].x = exp(-pAz[gid]*gdt);
    GpAz[gid].y = 0;

    EV[gid].x=cos(V[gid]*(dt/(2*HBAR)));
    EV[gid].y=sin(V[gid]*(dt/(2*HBAR)));
    EK[gid].x=cos(K[gid]*(dt/HBAR));
    EK[gid].y=sin(K[gid]*(dt/HBAR));

    EpAz[gid].x=cos(pAz[gid]*dt);
    EpAz[gid].y=sin(pAz[gid]*dt);
    EpAy[gid].x=cos(pAy[gid]*dt);
    EpAy[gid].y=sin(pAy[gid]*dt);
    EpAx[gid].x=cos(pAx[gid]*dt);
    EpAx[gid].y=sin(pAx[gid]*dt);
}

// Function to generate grids and treads for 2d and 3d cases
void generate_grid(Grid& par){

    int dimnum = par.ival("dimnum");
    int xDim = par.ival("xDim");
    int yDim = par.ival("yDim");
    int zDim = par.ival("zDim");
    int xD = 1, yD = 1, zD = 1;
    int max_threads = 256;
    if (xDim < max_threads){
        max_threads = xDim;
    }

    if (dimnum == 2){
        if (xDim <= max_threads){
            par.threads.x = xDim;
            par.threads.y = 1;
            par.threads.z = 1;
    
            xD = 1;
            yD = yDim;
            zD = 1;
        }
        else{
            int count = 0;
            int dim_tmp = xDim;
            while (dim_tmp > max_threads){
                count++;
                dim_tmp /= 2;
            }
    
            std::cout << "count is: " << count << '\n';
    
            par.threads.x = dim_tmp;
            par.threads.y = 1;
            par.threads.z = 1;
            xD = pow(2,count);
            yD = yDim;
            zD = 1;
        }

    }
    else if (dimnum == 3){

        if (xDim <= max_threads){
            par.threads.x = xDim;
            par.threads.y = 1;
            par.threads.z = 1;
    
            xD = 1;
            yD = yDim;
            zD = zDim;
        }
        else{
            int count = 0;
            int dim_tmp = xDim;
            while (dim_tmp > max_threads){
                count++;
                dim_tmp /= 2;
            }
    
            std::cout << "count is: " << count << '\n';
    
            par.threads.x = dim_tmp;
            par.threads.y = 1;
            par.threads.z = 1;
            xD = pow(2,count);
            yD = yDim;
            zD = zDim;
        }
    
    }
    else if (dimnum == 1){
        par.threads.x = xDim;
    }
    par.grid.x=xD;
    par.grid.y=yD;
    par.grid.z=zD;

    std::cout << "threads in x are: " << par.threads.x << '\n';
    std::cout << "dimensions are: " << par.grid.x << '\t' 
                                    << par.grid.y << '\t' 
                                    << par.grid.z << '\n';

}
