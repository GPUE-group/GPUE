#include "../include/operators.h"
#include "../include/kernels.h"
#include "../include/dynamic.h"

void laplacian(Grid &par, double2 *data, double2* out, int xDim, int yDim,
               int zDim, double dx, double dy, double dz){

    dim3 grid = par.grid;
    dim3 threads = par.threads;
    int gsize = xDim * yDim * zDim;

    double2 *temp_derivative;
    cudaMalloc((void **) &temp_derivative, sizeof(double2)*gsize);
    derive<<<grid, threads>>>(data, temp_derivative, 1, gsize, dx);
    derive<<<grid, threads>>>(temp_derivative, temp_derivative, 1, gsize, dx);

    copy<<<grid, threads>>>(temp_derivative, out);

    derive<<<grid, threads>>>(data, temp_derivative, xDim, gsize, dy);
    derive<<<grid, threads>>>(temp_derivative, temp_derivative,
                              xDim, gsize, dy);

    sum<<<grid, threads>>>(temp_derivative, out, out);

    derive<<<grid, threads>>>(data, temp_derivative, xDim*yDim, gsize, dz);
    derive<<<grid, threads>>>(temp_derivative, temp_derivative,
                              xDim*yDim, gsize, dz);

    sum<<<grid, threads>>>(temp_derivative, out, out);

    cudaFree(temp_derivative);

}

void laplacian(Grid &par, double2 *data, double2* out, int xDim, int yDim,
               double dx, double dy){


    dim3 grid = par.grid;
    dim3 threads = par.threads;
    int gsize = xDim * yDim;

    double2 *temp_derivative;
    cudaMalloc((void **) &temp_derivative, sizeof(double2)*gsize);
    derive<<<grid, threads>>>(data, temp_derivative, 1, gsize, dx);
    derive<<<grid, threads>>>(temp_derivative, temp_derivative, 1, gsize, dx);

    copy<<<grid, threads>>>(temp_derivative, out);

    derive<<<grid, threads>>>(data, temp_derivative, xDim, gsize, dy);
    derive<<<grid, threads>>>(temp_derivative, temp_derivative,
                              xDim, gsize, dy);

    sum<<<grid, threads>>>(temp_derivative, out, out);

    cudaFree(temp_derivative);
}

void laplacian(Grid &par, double2 *data, double2* out, int xDim, double dx){

    dim3 grid = par.grid;
    dim3 threads = par.threads;
    int gsize = xDim;

    derive<<<grid, threads>>>(data, out, 1, gsize, dx);
    derive<<<grid, threads>>>(out, out, 1, gsize, dx);

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
    curl = (double *)malloc(size);

    int index;

    // Note: To take the curl, we need a change in x and y to create a dx or dy
    //       For this reason, we have added yDim to y and 1 to x
    for (int i = 0; i < xDim; i++){
        for (int j = 0; j < yDim-1; j++){
            index = j + yDim * i;
            curl[index] = (Ay[index] - Ay[index+xDim]) 
                          - (Ax[index] - Ax[index+1]);
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
        curl[i] = atan2(By[i], Bx[i]);
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
    curl = (double *)malloc(size);

    int index;

    // Note: To take the curl, we need a change in x and y to create a dx or dy
    //       For this reason, we have added yDim to y and 1 to x
    for (int i = 0; i < xDim; i++){
        for (int j = 0; j < yDim-1; j++){
            for (int k = 0; k < zDim - 1; k++){
                index = k + zDim * j + zDim * yDim * i;
                curl[index] = (Az[index] - Az[index + xDim])
                              -(Ay[index] - Ay[index + xDim*yDim]);
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
    curl = (double *)malloc(size);

    int index;

    // Note: To take the curl, we need a change in x and y to create a dx or dy
    //       For this reason, we have added yDim to y and 1 to x
    for (int i = 0; i < xDim-1; i++){
        for (int j = 0; j < yDim; j++){
            for (int k = 0; k < zDim - 1; k++){
                index = k + zDim * j + zDim * yDim * i;
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
    curl = (double *)malloc(size);

    int index;

    // Note: To take the curl, we need a change in x and y to create a dx or dy
    //       For this reason, we have added yDim to y and 1 to x
    for (int i = 0; i < xDim - 1; i++){
        for (int j = 0; j < yDim-1; j++){
            for (int k = 0; k < zDim; k++){
                index = k + zDim * j + zDim * yDim * i;
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

// Function to read Ax from file.
// Note that this comes with a special method in init...
void file_A(std::string filename, double *A, double omega){
    std::fstream infile(filename, std::ios_base::in);

    double inval;
    int count = 0;
    while (infile >> inval){
        A[count] = omega*inval;
        count++;
    }
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
    cudaMalloc((void**) &x_gpu, sizeof(double) * xDim);
    cudaMalloc((void**) &y_gpu, sizeof(double) * yDim);
    cudaMalloc((void**) &z_gpu, sizeof(double) * zDim);
    cudaMalloc((void**) &px_gpu, sizeof(double) * xDim);
    cudaMalloc((void**) &py_gpu, sizeof(double) * yDim);
    cudaMalloc((void**) &pz_gpu, sizeof(double) * zDim);

    cudaMemcpy(x_gpu, x, sizeof(double)*xDim, cudaMemcpyHostToDevice);
    cudaMemcpy(y_gpu, y, sizeof(double)*yDim, cudaMemcpyHostToDevice);
    cudaMemcpy(z_gpu, z, sizeof(double)*zDim, cudaMemcpyHostToDevice);
    cudaMemcpy(px_gpu, px, sizeof(double)*xDim, cudaMemcpyHostToDevice);
    cudaMemcpy(py_gpu, py, sizeof(double)*yDim, cudaMemcpyHostToDevice);
    cudaMemcpy(pz_gpu, pz, sizeof(double)*zDim, cudaMemcpyHostToDevice);

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

    // Creating K to work with
    double *K, *K_gpu;
    K = (double*)malloc(sizeof(double)*gSize);
    cudaMalloc((void**) &K_gpu, sizeof(double)*gSize);

    simple_K<<<par.grid, par.threads>>>(px_gpu, py_gpu, pz_gpu, mass, K_gpu);

    cudaMemcpy(K, K_gpu, sizeof(double)*gSize, cudaMemcpyDeviceToHost);
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

    double *Ax, *Ay, *Az, *Ax_gpu, *Ay_gpu, *Az_gpu;
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

    Ax = (double *)malloc(sizeof(double)*gSize);
    Ay = (double *)malloc(sizeof(double)*gSize);
    Az = (double *)malloc(sizeof(double)*gSize);

    cudaMalloc((void**) &Ax_gpu, sizeof(double)*gSize);
    cudaMalloc((void**) &Ay_gpu, sizeof(double)*gSize);
    cudaMalloc((void**) &Az_gpu, sizeof(double)*gSize);

    if (par.Afn == "file"){
        file_A(par.Axfile, Ax, omega);
        cudaMemcpy(Ax_gpu, Ax, sizeof(double)*gSize, cudaMemcpyHostToDevice);

        if (dimnum > 1){
            file_A(par.Ayfile, Ay, omega);
            cudaMemcpy(Ay_gpu,Ay,sizeof(double)*gSize,cudaMemcpyHostToDevice);
        }

        if (dimnum == 3){
            file_A(par.Azfile, Az, omega);
            cudaMemcpy(Az_gpu,Az,sizeof(double)*gSize,cudaMemcpyHostToDevice);
        }

        std::cout << "finished reading Ax / Ay / Az from file" << '\n';
    }
    else{
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

            find_field<<<par.grid, par.threads>>>(Ax_gpu, dx, dy, dz, 
                                                  xMax, yMax, zMax, 0, eqn);
        }
        else{
            par.Ax_fn<<<par.grid, par.threads>>>(x_gpu, y_gpu, z_gpu, 
                                                 xMax, yMax, zMax, 
                                                 omegaX, omegaY, omegaZ, 
                                                 omega, fudge, Ax_gpu);
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

            find_field<<<par.grid, par.threads>>>(Ay_gpu, dx, dy, dz,
                                                  xMax, yMax, zMax , 0, eqn);
        }
        else{
            par.Ay_fn<<<par.grid, par.threads>>>(x_gpu, y_gpu, z_gpu, 
                                                 xMax, yMax, zMax, 
                                                 omegaX, omegaY, omegaZ, 
                                                 omega, fudge, Ay_gpu);
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

                find_field<<<par.grid, par.threads>>>(Az_gpu, dx, dy, dz,
                                                      xMax, yMax, zMax,  
                                                      0, eqn);
            }
            else{
                par.Az_fn<<<par.grid, par.threads>>>(x_gpu, y_gpu, z_gpu, 
                                                     xMax, yMax, zMax, 
                                                     omegaX, omegaY, omegaZ, 
                                                     omega, fudge, Az_gpu);
            }
        }
        else{
            kconstant_A<<<par.grid, par.threads>>>(x_gpu, y_gpu, z_gpu, 
                                                   xMax, yMax, zMax, 
                                                   omegaX, omegaY, omegaZ, 
                                                   omega, fudge, Az_gpu);
        }
    }
    cudaMemcpy(Ax, Ax_gpu, sizeof(double)*gSize,cudaMemcpyDeviceToHost);
    cudaMemcpy(Ay, Ay_gpu, sizeof(double)*gSize,cudaMemcpyDeviceToHost);
    cudaMemcpy(Az, Az_gpu, sizeof(double)*gSize,cudaMemcpyDeviceToHost);

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
    generate_gauge(par);

    int gSize = par.ival("gSize");
    int dimnum = par.ival("dimnum");
    int winding = par.dval("winding");

    bool energy_calc = par.bval("energy_calc");

    double dt = par.dval("dt");
    double gdt = par.dval("gdt");
    double *x_gpu = par.dsval("x_gpu");
    double *y_gpu = par.dsval("y_gpu");
    double *z_gpu = par.dsval("z_gpu");
    double *px_gpu = par.dsval("px_gpu");
    double *py_gpu = par.dsval("py_gpu");
    double *pz_gpu = par.dsval("pz_gpu");
    double *Ax_gpu = par.dsval("Ax_gpu");
    double *Ay_gpu = par.dsval("Ay_gpu");
    double *Az_gpu = par.dsval("Az_gpu");
    double *K_gpu = par.dsval("K_gpu");

    // Creating items list for kernels

    double *items, *items_gpu;
    int item_size = 18;
    items = (double*)malloc(sizeof(double)*item_size);
    cudaMalloc((void**) &items_gpu, sizeof(double)*item_size);

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

    cudaMemcpy(items_gpu, items, sizeof(double)*item_size,
               cudaMemcpyHostToDevice);

    double fudge = par.dval("fudge");

    // Generating V

    double *V, *V_gpu;

    V = (double *)malloc(sizeof(double)*gSize);

    cudaMalloc((void **) &V_gpu, sizeof(double)*gSize);

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
        find_field<<<par.grid, par.threads>>>(V_gpu, dx, dy, dz, 
                                              xMax, yMax, zMax, 0, eqn);
    }
    else{
        par.V_fn<<<par.grid, par.threads>>>(x_gpu, y_gpu, z_gpu, items_gpu,
                                            Ax_gpu, Ay_gpu, Az_gpu, V_gpu);
    }

    cudaMemcpy(V, V_gpu, sizeof(double)*gSize, cudaMemcpyDeviceToHost);

    // Generating wfc

    double2 *wfc, *wfc_gpu;
    double *phi, *phi_gpu;

    wfc = (double2 *)malloc(sizeof(double2)*gSize);
    phi = (double *)malloc(sizeof(double)*gSize);

    cudaMalloc((void**) &wfc_gpu, sizeof(double2)*gSize);
    cudaMalloc((void**) &phi_gpu, sizeof(double)*gSize);

    if (par.bval("read_wfc")){
        wfc = par.cufftDoubleComplexval("wfc");
        cudaMemcpy(wfc_gpu, wfc, sizeof(double2)*gSize, cudaMemcpyHostToDevice);
    }
    else{
        par.wfc_fn<<<par.grid, par.threads>>>(x_gpu, y_gpu, z_gpu, items_gpu,
                                              winding, phi_gpu, wfc_gpu);
        cudaMemcpy(wfc, wfc_gpu, sizeof(double2)*gSize, cudaMemcpyDeviceToHost);
    }

    cudaMemcpy(phi, phi_gpu, sizeof(double)*gSize, cudaMemcpyDeviceToHost);

    // generating aux fields.
    double2 *GV, *EV, *GK, *EK;
    double2 *GV_gpu, *EV_gpu, *GK_gpu, *EK_gpu;
    double2 *GpAx, *GpAy, *GpAz, *EpAx, *EpAy, *EpAz;
    double2 *GpAx_gpu, *GpAy_gpu, *GpAz_gpu, *EpAx_gpu, *EpAy_gpu, *EpAz_gpu;
    double *pAx, *pAy, *pAz;
    double *pAx_gpu, *pAy_gpu, *pAz_gpu;

    GV = (double2 *)malloc(sizeof(double2)*gSize);
    EV = (double2 *)malloc(sizeof(double2)*gSize);
    GK = (double2 *)malloc(sizeof(double2)*gSize);
    EK = (double2 *)malloc(sizeof(double2)*gSize);

    GpAx = (double2 *)malloc(sizeof(double2)*gSize);
    EpAx = (double2 *)malloc(sizeof(double2)*gSize);
    GpAy = (double2 *)malloc(sizeof(double2)*gSize);
    EpAy = (double2 *)malloc(sizeof(double2)*gSize);
    GpAz = (double2 *)malloc(sizeof(double2)*gSize);
    EpAz = (double2 *)malloc(sizeof(double2)*gSize);

    pAx = (double *)malloc(sizeof(double)*gSize);
    pAy = (double *)malloc(sizeof(double)*gSize);
    pAz = (double *)malloc(sizeof(double)*gSize);

    cudaMalloc((void**) &GV_gpu, sizeof(double2)*gSize);
    cudaMalloc((void**) &EV_gpu, sizeof(double2)*gSize);
    cudaMalloc((void**) &GK_gpu, sizeof(double2)*gSize);
    cudaMalloc((void**) &EK_gpu, sizeof(double2)*gSize);

    cudaMalloc((void**) &GpAx_gpu, sizeof(double2)*gSize);
    cudaMalloc((void**) &EpAx_gpu, sizeof(double2)*gSize);
    cudaMalloc((void**) &GpAy_gpu, sizeof(double2)*gSize);
    cudaMalloc((void**) &EpAy_gpu, sizeof(double2)*gSize);
    cudaMalloc((void**) &GpAz_gpu, sizeof(double2)*gSize);
    cudaMalloc((void**) &EpAz_gpu, sizeof(double2)*gSize);

    cudaMalloc((void**) &pAx_gpu, sizeof(double)*gSize);
    cudaMalloc((void**) &pAy_gpu, sizeof(double)*gSize);
    cudaMalloc((void**) &pAz_gpu, sizeof(double)*gSize);

    aux_fields<<<par.grid, par.threads>>>(V_gpu, K_gpu, gdt, dt,
                                          Ax_gpu, Ay_gpu, Az_gpu,
                                          px_gpu, py_gpu, pz_gpu,
                                          pAx_gpu, pAy_gpu, pAz_gpu,
                                          GV_gpu, EV_gpu, GK_gpu, EK_gpu,
                                          GpAx_gpu, GpAy_gpu, GpAz_gpu,
                                          EpAx_gpu, EpAy_gpu, EpAz_gpu);

    cudaMemcpy(GV, GV_gpu, sizeof(double2)*gSize, cudaMemcpyDeviceToHost);
    cudaMemcpy(EV, EV_gpu, sizeof(double2)*gSize, cudaMemcpyDeviceToHost);
    cudaMemcpy(GK, GK_gpu, sizeof(double2)*gSize, cudaMemcpyDeviceToHost);
    cudaMemcpy(EK, EK_gpu, sizeof(double2)*gSize, cudaMemcpyDeviceToHost);

    cudaMemcpy(GpAx, GpAx_gpu, sizeof(double2)*gSize, cudaMemcpyDeviceToHost);
    cudaMemcpy(EpAx, EpAx_gpu, sizeof(double2)*gSize, cudaMemcpyDeviceToHost);
    cudaMemcpy(GpAy, GpAy_gpu, sizeof(double2)*gSize, cudaMemcpyDeviceToHost);
    cudaMemcpy(EpAy, EpAy_gpu, sizeof(double2)*gSize, cudaMemcpyDeviceToHost);
    cudaMemcpy(GpAz, GpAz_gpu, sizeof(double2)*gSize, cudaMemcpyDeviceToHost);
    cudaMemcpy(EpAz, EpAz_gpu, sizeof(double2)*gSize, cudaMemcpyDeviceToHost);

    cudaMemcpy(pAx, pAx_gpu, sizeof(double)*gSize, cudaMemcpyDeviceToHost);
    cudaMemcpy(pAy, pAy_gpu, sizeof(double)*gSize, cudaMemcpyDeviceToHost);
    cudaMemcpy(pAz, pAz_gpu, sizeof(double)*gSize, cudaMemcpyDeviceToHost);

    // Storing variables
    cudaFree(items_gpu);
    //cudaFree(phi_gpu);
    cudaFree(GV_gpu);
    cudaFree(EV_gpu);
    cudaFree(GK_gpu);
    cudaFree(EK_gpu);

    cudaFree(pAx_gpu);
    cudaFree(pAy_gpu);
    cudaFree(pAz_gpu);

    cudaFree(GpAx_gpu);
    cudaFree(GpAy_gpu);
    cudaFree(GpAz_gpu);

    cudaFree(EpAx_gpu);
    cudaFree(EpAy_gpu);
    cudaFree(EpAz_gpu);

    cudaFree(x_gpu);
    cudaFree(y_gpu);
    cudaFree(z_gpu);

    cudaFree(px_gpu);
    cudaFree(py_gpu);
    cudaFree(pz_gpu);

    if (!energy_calc){
        cudaFree(K_gpu);
        cudaFree(V_gpu);
        cudaFree(Ax_gpu);
        cudaFree(Ay_gpu);
        cudaFree(Az_gpu);
    }
    else{
        par.store("V_gpu",V_gpu);
    }

    par.store("V",V);
    par.store("items", items);
    //par.store("items_gpu", items_gpu);
    par.store("wfc", wfc);
    par.store("wfc_gpu", wfc_gpu);
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
                      - 0.4*items[0];
    double omegaR = (items[3]*items[3] + items[4]*items[4]);
    double V_tot = (2*items[5]*items[5]*(z[zid] - items[8])*(z[zid] - items[8])
                    + omegaR*(rad*rad + items[12]*rad*z[zid]));
    V[gid] = 0.5*items[9]*(V_tot
                           + Ax[gid]*Ax[gid]
                           + Ay[gid]*Ay[gid]
                           + Az[gid]*Az[gid]);
}

__global__ void kstd_wfc(double *x, double *y, double *z, double *items,
                         double winding, double *phi, double2 *wfc){

    int gid = getGid3d3d();
    int xid = blockDim.x*blockIdx.x + threadIdx.x;
    int yid = blockDim.y*blockIdx.y + threadIdx.y;
    int zid = blockDim.z*blockIdx.z + threadIdx.z;

    phi[gid] = fmod(winding*atan2(y[yid], x[xid]),2*PI);

    wfc[gid].x = exp(-(x[xid]*x[xid]/(items[14]*items[14]*items[15]*items[15]) 
                     + y[yid]*y[yid]/(items[14]*items[14]*items[16]*items[16]) 
                     + z[zid]*z[zid]/(items[14]*items[14]*items[17]*items[17])))
                     * cos(phi[gid]);
    wfc[gid].y = -exp(-(x[xid]*x[xid]/(items[14]*items[14]*items[15]*items[15]) 
                     + y[yid]*y[yid]/(items[14]*items[14]*items[16]*items[16]) 
                     + z[zid]*z[zid]/(items[14]*items[14]*items[17]*items[17])))
                     * sin(phi[gid]);

}

__global__ void ktorus_wfc(double *x, double *y, double *z, double *items,
                           double winding, double *phi, double2 *wfc){

    int gid = getGid3d3d();
    int xid = blockDim.x*blockIdx.x + threadIdx.x;
    int yid = blockDim.y*blockIdx.y + threadIdx.y;
    int zid = blockDim.z*blockIdx.z + threadIdx.z;

    double rad = sqrt((x[xid] - items[6]) * (x[xid] - items[6])
                      + (y[yid] - items[7]) * (y[yid] - items[7])) 
                      - 0.5*items[0]*items[12];

    wfc[gid].x = exp(-( pow((rad)/(items[14]*items[15]*0.5),2) +
                   pow((z[zid])/(items[14]*items[17]*0.5),2) ) );
    wfc[gid].y = 0.0;
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

    GV[gid].x = exp( -V[gid]*(gdt/(2*HBAR)));
    GK[gid].x = exp( -K[gid]*(gdt/HBAR));
    GV[gid].y = 0.0;
    GK[gid].y = 0.0;

    // Ax and Ay will be calculated here but are used only for
    // debugging. They may be needed later for magnetic field calc

    pAy[gid] = Ax[gid] * px[xid];
    pAx[gid] = Ay[gid] * py[yid];
    pAz[gid] = Az[gid] * pz[zid];

    GpAx[gid].x = exp(-pAx[gid]*gdt);
    GpAx[gid].y = 0;
    GpAy[gid].x = exp(-pAy[gid]*gdt);
    GpAy[gid].y = 0;
    GpAz[gid].x = exp(-pAz[gid]*gdt);
    GpAz[gid].y = 0;

    EV[gid].x=cos( -V[gid]*(dt/(2*HBAR)));
    EV[gid].y=sin( -V[gid]*(dt/(2*HBAR)));
    EK[gid].x=cos( -K[gid]*(dt/HBAR));
    EK[gid].y=sin( -K[gid]*(dt/HBAR));

    EpAz[gid].x=cos(-pAz[gid]*dt);
    EpAz[gid].y=sin(-pAz[gid]*dt);
    EpAy[gid].x=cos(-pAy[gid]*dt);
    EpAy[gid].y=sin(-pAy[gid]*dt);
    EpAx[gid].x=cos(-pAx[gid]*dt);
    EpAx[gid].y=sin(-pAx[gid]*dt);
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
