/*-------------vortex_3d.cu---------------------------------------------------//
*
* Purpose: This file intends to perform a convolution of 3d data for vortex
*          recognition in the GPUE code
*
*   Notes: We will be using the window method for convolutions because it has
*          a slightly better complxity case for a separable filter
*          (which the Sobel filter definitely is!)
*
*-----------------------------------------------------------------------------*/

#include <algorithm>
#include "../include/ds.h"
#include "../include/vortex_3d.h"
#include "../include/split_op.h"
#include "../include/kernels.h"

//We will need a few functions to deal with vortex skeletons
pos **find_skeletons(double* edges);

// Function to find the sobel operators and transfer to GPU
void find_sobel(Grid &par){

    std::string conv_type = par.sval("conv_type");
    int xDim, yDim, zDim;

    // There will be two cases to take into account here, one for fft 
    // convolution and another for window

    int index = 0;

    if (conv_type == "FFT"){
        double2 *sobel_x, *sobel_y, *sobel_z;
        xDim = par.ival("xDim");
        yDim = par.ival("yDim");
        zDim = par.ival("zDim");

        sobel_x = (double2 *) malloc(sizeof(double2) *xDim*yDim*zDim);
        sobel_y = (double2 *) malloc(sizeof(double2) *xDim*yDim*zDim);
        sobel_z = (double2 *) malloc(sizeof(double2) *xDim*yDim*zDim);

        // Now let's go ahead and pad these guys with 0's
        for (int i = 0; i < xDim; ++i){
            for (int j = 0; j < yDim; ++j){
                for (int k = 0; k < zDim; ++k){
                    index = k + j * xDim + i * yDim * zDim;
                    sobel_x[index].x = 0;
                    sobel_y[index].x = 0;
                    sobel_z[index].x = 0;
                }
            }
        }

        // There is clearly a better way to do this with matrix multiplication
        // the sobel operator is separable, so we just need to mix
        // Gradient and triangle filters, check:
        //     https://en.wikipedia.org/wiki/Sobel_operator
    
        // Starting with S_z
        int factor;
        for (int i = 0; i < 3; ++i){
            for (int j = 0; j < 3; ++j){
                for (int k = 0; k < 3; ++k){
                    index = k + zDim*j + zDim*yDim*i;
                    if (k == 0){
                        factor = -1;
                    }
                    if (k == 1){
                        factor = 0;
                    }
                    if (k == 2){
                        factor = 1;
                    }
                    if (i == 0 && j == 0){
                        sobel_z[index].x = factor * 1;
                    }
                    if (i == 0 && j == 1){
                        sobel_z[index].x = factor * 2;
                    }
                    if (i == 0 && j == 2){
                        sobel_z[index].x = factor * 1;
                    }
                    if (i == 1 && j == 0){
                        sobel_z[index].x = factor * 2;
                    }
                    if (i == 1 && j == 1){
                        sobel_z[index].x = factor * 4;
                    }
                    if (i == 1 && j == 2){
                        sobel_z[index].x = factor * 2;
                    }
                    if (i == 2 && j == 0){
                        sobel_z[index].x = factor * 1;
                    }
                    if (i == 2 && j == 1){
                        sobel_z[index].x = factor * 2;
                    }
                    if (i == 2 && j == 2){
                        sobel_z[index].x = factor * 1;
                    }
                }
            }
        }
    
        // S_y
        for (int i = 0; i < 3; ++i){
            for (int j = 0; j < 3; ++j){
                for (int k = 0; k < 3; ++k){
                    index = k + zDim*j + zDim*yDim*i;
                    if (j == 0){
                        factor = -1;
                    }
                    if (j == 1){
                        factor = 0;
                    }
                    if (j == 2){
                        factor = 1;
                    }
                    if (i == 0 && k == 0){
                        sobel_y[index].x = factor * 1;
                    }
                    if (i == 0 && k == 1){
                        sobel_y[index].x = factor * 2;
                    }
                    if (i == 0 && k == 2){
                        sobel_y[index].x = factor * 1;
                    }
                    if (i == 1 && k == 0){
                        sobel_y[index].x = factor * 2;
                    }
                    if (i == 1 && k == 1){
                        sobel_y[index].x = factor * 4;
                    }
                    if (i == 1 && k == 2){
                        sobel_y[index].x = factor * 2;
                    }
                    if (i == 2 && k == 0){
                        sobel_y[index].x = factor * 1;
                    }
                    if (i == 2 && k == 1){
                        sobel_y[index].x = factor * 2;
                    }
                    if (i == 2 && k == 2){
                        sobel_y[index].x = factor * 1;
                    }
                }
            }
        }
    
        // Now for S_x
        for (int i = 0; i < 3; ++i){
            for (int j = 0; j < 3; ++j){
                for (int k = 0; k < 3; ++k){
                    index = k + zDim*j + zDim*yDim*i;
                    if (i == 0){
                        factor = -1;
                    }
                    if (i == 1){
                        factor = 0;
                    }
                    if (i == 2){
                        factor = 1;
                    }
                    if (k == 0 && j == 0){
                        sobel_x[index].x = factor * 1;
                    }
                    if (k == 0 && j == 1){
                        sobel_x[index].x = factor * 2;
                    }
                    if (k == 0 && j == 2){
                        sobel_x[index].x = factor * 1;
                    }
                    if (k == 1 && j == 0){
                        sobel_x[index].x = factor * 2;
                    }
                    if (k == 1 && j == 1){
                        sobel_x[index].x = factor * 4;
                    }
                    if (k == 1 && j == 2){
                        sobel_x[index].x = factor * 2;
                    }
                    if (k == 2 && j == 0){
                        sobel_x[index].x = factor * 1;
                    }
                    if (k == 2 && j == 1){
                        sobel_x[index].x = factor * 2;
                    }
                    if (k == 2 && j == 2){
                        sobel_x[index].x = factor * 1;
                    }
                }
            }
        }


        par.store("sobel_x", sobel_x);
        par.store("sobel_y", sobel_y);
        par.store("sobel_z", sobel_z);

    }
    else{
        double *sobel_x, *sobel_y, *sobel_z;
        xDim = 3;
        yDim = 3;
        zDim = 3;

        sobel_x = (double *) malloc(sizeof(double) *27);
        sobel_y = (double *) malloc(sizeof(double) *27);
        sobel_z = (double *) malloc(sizeof(double) *27);

        // There is clearly a better way to do this with matrix multiplication
        // the sobel operator is separable, so we just need to mix 
        // Gradient and triangle filters, check:
        //     https://en.wikipedia.org/wiki/Sobel_operator
    
        // Starting with S_z
        int factor;
        for (int i = 0; i < 3; ++i){
            for (int j = 0; j < 3; ++j){
                for (int k = 0; k < 3; ++k){
                    index = k + 3*j + 9*i;
                    if (k == 0){
                        factor = -1;
                    }
                    if (k == 1){
                        factor = 0;
                    }
                    if (k == 2){
                        factor = 1;
                    }
                    if (i == 0 && j == 0){
                        sobel_z[index] = factor * 1;
                    }
                    if (i == 0 && j == 1){
                        sobel_z[index] = factor * 2;
                    }
                    if (i == 0 && j == 2){
                        sobel_z[index] = factor * 1;
                    }
                    if (i == 1 && j == 0){
                        sobel_z[index] = factor * 2;
                    }
                    if (i == 1 && j == 1){
                        sobel_z[index] = factor * 4;
                    }
                    if (i == 1 && j == 2){
                        sobel_z[index] = factor * 2;
                    }
                    if (i == 2 && j == 0){
                        sobel_z[index] = factor * 1;
                    }
                    if (i == 2 && j == 1){
                        sobel_z[index] = factor * 2;
                    }
                    if (i == 2 && j == 2){
                        sobel_z[index] = factor * 1;
                    }
                }
            }
        }
    
        // S_y
        for (int i = 0; i < 3; ++i){
            for (int j = 0; j < 3; ++j){
                for (int k = 0; k < 3; ++k){
                    index = k + 3*j + 9*i;
                    if (j == 0){
                        factor = -1;
                    }
                    if (j == 1){
                        factor = 0;
                    }
                    if (j == 2){
                        factor = 1;
                    }
                    if (i == 0 && k == 0){
                        sobel_y[index] = factor * 1;
                    }
                    if (i == 0 && k == 1){
                        sobel_y[index] = factor * 2;
                    }
                    if (i == 0 && k == 2){
                        sobel_y[index] = factor * 1;
                    }
                    if (i == 1 && k == 0){
                        sobel_y[index] = factor * 2;
                    }
                    if (i == 1 && k == 1){
                        sobel_y[index] = factor * 4;
                    }
                    if (i == 1 && k == 2){
                        sobel_y[index] = factor * 2;
                    }
                    if (i == 2 && k == 0){
                        sobel_y[index] = factor * 1;
                    }
                    if (i == 2 && k == 1){
                        sobel_y[index] = factor * 2;
                    }
                    if (i == 2 && k == 2){
                        sobel_y[index] = factor * 1;
                    }
                }
            }
        }
    
        // Now for S_x
        for (int i = 0; i < 3; ++i){
            for (int j = 0; j < 3; ++j){
                for (int k = 0; k < 3; ++k){
                    index = k + 3*j + 9*i;
                    if (i == 0){
                        factor = -1;
                    }
                    if (i == 1){
                        factor = 0;
                    }
                    if (i == 2){
                        factor = 1;
                    }
                    if (k == 0 && j == 0){
                        sobel_x[index] = factor * 1;
                    }
                    if (k == 0 && j == 1){
                        sobel_x[index] = factor * 2;
                    }
                    if (k == 0 && j == 2){
                        sobel_x[index] = factor * 1;
                    }
                    if (k == 1 && j == 0){
                        sobel_x[index] = factor * 2;
                    }
                    if (k == 1 && j == 1){
                        sobel_x[index] = factor * 4;
                    }
                    if (k == 1 && j == 2){
                        sobel_x[index] = factor * 2;
                    }
                    if (k == 2 && j == 0){
                        sobel_x[index] = factor * 1;
                    }
                    if (k == 2 && j == 1){
                        sobel_x[index] = factor * 2;
                    }
                    if (k == 2 && j == 2){
                        sobel_x[index] = factor * 1;
                    }
                }
            }
        }



        par.store("sobel_x", sobel_x);
        par.store("sobel_y", sobel_y);
        par.store("sobel_z", sobel_z);
    }

    transfer_sobel(par);

    par.store("found_sobel", true);
}

void find_sobel_2d(Grid &par){

    std::string conv_type = par.sval("conv_type");
    int xDim, yDim;

    // There will be two cases to take into account here, one for fft 
    // convolution and another for window

    int index = 0;

    if (conv_type == "FFT"){
        double2 *sobel_x, *sobel_y;
        xDim = par.ival("xDim");
        yDim = par.ival("yDim");

        sobel_x = (double2 *) malloc(sizeof(double2) *xDim*yDim);
        sobel_y = (double2 *) malloc(sizeof(double2) *xDim*yDim);

        // Now let's go ahead and pad these guys with 0's
        for (int i = 0; i < xDim; ++i){
            for (int j = 0; j < yDim; ++j){
                index = j + i * yDim;
                sobel_x[index].x = 0;
                sobel_y[index].x = 0;
            }
        }

        // There is clearly a better way to do this with matrix multiplication
        // the sobel operator is separable, so we just need to mix
        // Gradient and triangle filters, check:
        //     https://en.wikipedia.org/wiki/Sobel_operator
    
        // S_y
        int factor;
        for (int i = 0; i < 3; ++i){
            for (int j = 0; j < 3; ++j){
                index = j + yDim*i;
                if (j == 0){
                    factor = -1;
                }
                if (j == 1){
                    factor = 0;
                }
                if (j == 2){
                    factor = 1;
                }
                if (j == 0 && i == 0){
                    sobel_y[index].x = factor * 1;
                }
                if (j == 0 && i == 1){
                    sobel_y[index].x = factor * 2;
                }
                if (j == 0 && i == 2){
                    sobel_y[index].x = factor * 1;
                }
                if (j == 2 && i == 0){
                    sobel_y[index].x = factor * 1;
                }
                if (j == 2 && i == 1){
                    sobel_y[index].x = factor * 2;
                }
                if (j == 2 && i == 2){
                    sobel_y[index].x = factor * 1;
                }
            }
        }
    
        // Now for S_x
        for (int i = 0; i < 3; ++i){
            for (int j = 0; j < 3; ++j){
                index = j + yDim*i;
                if (i == 0){
                    factor = -1;
                }
                if (i == 1){
                    factor = 0;
                }
                if (i == 2){
                    factor = 1;
                }
                if (i == 0 && j == 0){
                    sobel_x[index].x = factor * 1;
                }
                if (i == 0 && j == 1){
                    sobel_x[index].x = factor * 2;
                }
                if (i == 0 && j == 2){
                    sobel_x[index].x = factor * 1;
                }
                if (i == 2 && j == 0){
                    sobel_x[index].x = factor * 1;
                }
                if (i == 2 && j == 1){
                    sobel_x[index].x = factor * 2;
                }
                if (i == 2 && j == 2){
                    sobel_x[index].x = factor * 1;
                }
            }
        }


        par.store("sobel_x", sobel_x);
        par.store("sobel_y", sobel_y);
    }
    else{
        double *sobel_x, *sobel_y;
        xDim = 3;
        yDim = 3;

        sobel_x = (double *) malloc(sizeof(double) *9);
        sobel_y = (double *) malloc(sizeof(double) *9);

        // There is clearly a better way to do this with matrix multiplication
        // the sobel operator is separable, so we just need to mix 
        // Gradient and triangle filters, check:
        //     https://en.wikipedia.org/wiki/Sobel_operator
    
        int factor;
        // S_y
        for (int i = 0; i < 3; ++i){
            for (int j = 0; j < 3; ++j){
                index = j + 3*i;
                if (j == 0){
                    factor = -1;
                }
                if (j == 1){
                    factor = 0;
                }
                if (j == 2){
                    factor = 1;
                }
                if (j == 0 && i == 0){
                    sobel_y[index] = factor * 1;
                }
                if (j == 0 && i == 1){
                    sobel_y[index] = factor * 2;
                }
                if (j == 0 && i == 2){
                    sobel_y[index] = factor * 1;
                }
                if (j == 2 && i == 0){
                    sobel_y[index] = factor * 1;
                }
                if (j == 2 && i == 1){
                    sobel_y[index] = factor * 2;
                }
                if (j == 2 && i == 2){
                    sobel_y[index] = factor * 1;
                }
            }
        }
    
        // Now for S_x
        for (int i = 0; i < 3; ++i){
            for (int j = 0; j < 3; ++j){
                index = j + 3*i;
                if (i == 0){
                    factor = -1;
                }
                if (i == 1){
                    factor = 0;
                }
                if (i == 2){
                    factor = 1;
                }
                if (i == 0 && j == 0){
                    sobel_x[index] = factor * 1;
                }
                if (i == 0 && j == 1){
                    sobel_x[index] = factor * 2;
                }
                if (i == 0 && j == 2){
                    sobel_x[index] = factor * 1;
                }
                if (i == 2 && j == 0){
                    sobel_x[index] = factor * 1;
                }
                if (i == 2 && j == 1){
                    sobel_x[index] = factor * 2;
                }
                if (i == 2 && j == 2){
                    sobel_x[index] = factor * 1;
                }
            }
        }



        par.store("sobel_x", sobel_x);
        par.store("sobel_y", sobel_y);
    }

    transfer_sobel(par);

}

// Function to transfer 3d sobel operators for non-fft convolution
void transfer_sobel(Grid &par){

    // Grabbing necessary parameters
    int dimnum = par.ival("dimnum");
    int xDim = par.ival("xDim");
    int yDim = par.ival("yDim");
    int zDim = par.ival("zDim");
    int gSize = xDim * yDim * zDim;
    std::string conv_type = par.sval("conv_type");

    double2 *sobel_x_gpu, *sobel_y_gpu, *sobel_z_gpu;

    cufftHandle plan_nd;
    if (dimnum == 2){
        plan_nd = par.ival("plan_2d");
    }
    else if (dimnum == 3){
        plan_nd = par.ival("plan_3d");
    }

    // creating space on device for 2 separate cases
    // Note that in the case of the FFT, the Sobel operators will be double2's
    //     while in the case of the window transform, they will be double
    if (conv_type == "FFT"){
        double2 *sobel_x = par.cufftDoubleComplexval("sobel_x");
        double2 *sobel_y = par.cufftDoubleComplexval("sobel_y");
        double2 *sobel_z = par.cufftDoubleComplexval("sobel_z");
        cudaHandleError( cudaMalloc((void**) &sobel_x_gpu, sizeof(double2) *gSize) );
        cudaHandleError( cudaMalloc((void**) &sobel_y_gpu, sizeof(double2) *gSize) );
        if (dimnum == 3){
            cudaHandleError( cudaMalloc((void**) &sobel_z_gpu, sizeof(double2) *gSize) );
        }

        // Transferring to device
    
        // Sobel_x
        cudaHandleError( cudaMemcpy(sobel_x_gpu, sobel_x, sizeof(double2)*gSize,
                                    cudaMemcpyHostToDevice));
    
        // Sobel_y
        cudaHandleError( cudaMemcpy(sobel_y_gpu, sobel_y, sizeof(double2)*gSize,
                                    cudaMemcpyHostToDevice) );
    
        if (dimnum == 3){
        // Sobel_z
        cudaHandleError( cudaMemcpy(sobel_z_gpu, sobel_z, sizeof(double2)*gSize,
                                    cudaMemcpyHostToDevice) );
        }
    
        // We only need the FFT's of the sobel operators. Let's generate those
        cufftHandleError( cufftExecZ2Z(plan_nd, sobel_x_gpu, sobel_x_gpu, CUFFT_FORWARD) );
        cufftHandleError( cufftExecZ2Z(plan_nd, sobel_y_gpu, sobel_y_gpu, CUFFT_FORWARD) );
        if (dimnum == 3){
            cufftHandleError( cufftExecZ2Z(plan_nd, sobel_z_gpu, sobel_z_gpu, CUFFT_FORWARD) );
        }
    
        // Storing in set of parameters
        par.store("sobel_x_gpu", sobel_x_gpu);
        par.store("sobel_y_gpu", sobel_y_gpu);
        par.store("sobel_z_gpu", sobel_z_gpu);
    }
    else{
        double2 *sobel_x = par.cufftDoubleComplexval("sobel_x");
        double2 *sobel_y = par.cufftDoubleComplexval("sobel_y");
        double2 *sobel_z = par.cufftDoubleComplexval("sobel_z");
        int size;
        if (dimnum == 2){
            size = 9;
        }
        else if (dimnum == 3){
            size = 27;
        }
        cudaHandleError( cudaMalloc((void**) &sobel_x_gpu, sizeof(double) *size) );
        cudaHandleError( cudaMalloc((void**) &sobel_y_gpu, sizeof(double) *size) );
        cudaHandleError( cudaMalloc((void**) &sobel_z_gpu, sizeof(double) *size) );
        // Transferring to device

        // Sobel_x
        cudaHandleError( cudaMemcpy(sobel_x_gpu, sobel_x, sizeof(double)*gSize,
                                    cudaMemcpyHostToDevice) );

        // Sobel_y
        cudaHandleError( cudaMemcpy(sobel_y_gpu, sobel_y, sizeof(double)*gSize,
                                    cudaMemcpyHostToDevice));

        if (dimnum == 3){
            // Sobel_z
            cudaHandleError( cudaMemcpy(sobel_z_gpu, sobel_z, sizeof(double)*gSize,
                                        cudaMemcpyHostToDevice) );
        }

        // Storing in set of parameters
        par.store("sobel_x_gpu", sobel_x_gpu);
        par.store("sobel_y_gpu", sobel_y_gpu);
        par.store("sobel_z_gpu", sobel_z_gpu);

    }

    par.store("found_sobel",true);

}

// function to transform a wavefunction to a field of edges
void find_edges(Grid &par,
                double2* wfc, double* edges){

    // for this, we simply need to take our sobel 3d sobel filter,
    // FFT forward, multiply, FFT back.

    int dimnum = par.ival("dimnum");
    dim3 grid = par.grid;
    dim3 threads = par.threads;

    double2 *wfc_gpu = par.cufftDoubleComplexval("wfc_gpu");
    int xDim = par.ival("xDim");
    int yDim = par.ival("yDim");
    int zDim = par.ival("zDim");
    int gSize = xDim * yDim * zDim;

    double *density_d, *edges_gpu;
    double2 *density_d2, *gradient_x_fft, *gradient_y_fft, *gradient_z_fft;

    // Now we need to grab the Sobel operators
    if (par.bval("found_sobel") == false){
        std::cout << "Finding sobel filters" << '\n';
        if (dimnum == 2){
            find_sobel_2d(par);
        }
        if (dimnum == 3){
            find_sobel(par);
        }
        cudaHandleError( cudaMalloc((void**) &density_d, sizeof(double) * gSize) );
        cudaHandleError( cudaMalloc((void**) &gradient_x_fft, sizeof(double2) * gSize) );
        cudaHandleError( cudaMalloc((void**) &gradient_y_fft, sizeof(double2) * gSize) );
        if (dimnum == 3){
            cudaHandleError( cudaMalloc((void**) &gradient_z_fft, sizeof(double2) * gSize) );
        }
        cudaHandleError( cudaMalloc((void**) &density_d2, sizeof(double2) * gSize) );
        cudaHandleError( cudaMalloc((void**) &edges_gpu, sizeof(double) * gSize) );
    }
    else{
        density_d = par.dsval("density_d");
        edges_gpu = par.dsval("edges_gpu");
        gradient_x_fft = par.cufftDoubleComplexval("gradient_x_fft");
        gradient_y_fft = par.cufftDoubleComplexval("gradient_y_fft");
        if (dimnum == 3){
            gradient_z_fft = par.cufftDoubleComplexval("gradient_z_fft");
        }
        density_d2 = par.cufftDoubleComplexval("density_d2");
    }

    cufftHandle plan_3d = par.ival("plan_3d");

    // now to perform the complexMagnitudeSquared operation
    complexMagnitudeSquared<<<grid,threads>>>(wfc_gpu, density_d);
    cudaCheckError();

    // Pulling operators from find_sobel(par)
    double2 *sobel_x_gpu = par.cufftDoubleComplexval("sobel_x_gpu");
    double2 *sobel_y_gpu = par.cufftDoubleComplexval("sobel_y_gpu");
    double2 *sobel_z_gpu;
    if (dimnum == 3){
        sobel_z_gpu = par.cufftDoubleComplexval("sobel_z_gpu");
    }

    // This should work in principle, but the D2Z transform plays tricks
    // Generating plan for d2z in 3d
    //cufftHandle plan_3d2z;
    //cufftPlan3d(&plan_3d2z, xDim, yDim, zDim, CUFFT_D2Z);

    make_cufftDoubleComplex<<<grid, threads>>>(density_d, density_d2);
    cudaCheckError();

    // Now fft forward, multiply, fft back
    cufftHandleError( cufftExecZ2Z(plan_3d, density_d2, gradient_x_fft, CUFFT_FORWARD) );
    cufftHandleError( cufftExecZ2Z(plan_3d, density_d2, gradient_y_fft, CUFFT_FORWARD) );
    if (dimnum == 3){
        cufftHandleError( cufftExecZ2Z(plan_3d, density_d2, gradient_z_fft, CUFFT_FORWARD) );
    }

    // Now to perform the multiplication
    cMult<<<grid, threads>>>(gradient_x_fft, sobel_x_gpu, gradient_x_fft);
    cudaCheckError();
    cMult<<<grid, threads>>>(gradient_y_fft, sobel_y_gpu, gradient_y_fft);
    cudaCheckError();
    if (dimnum == 3){
        cMult<<<grid, threads>>>(gradient_z_fft, sobel_z_gpu, gradient_z_fft);
        cudaCheckError();
    }
    
    // FFT back
    cufftHandleError( cufftExecZ2Z(plan_3d, gradient_x_fft, gradient_x_fft, CUFFT_INVERSE) );
    cufftHandleError( cufftExecZ2Z(plan_3d, gradient_y_fft, gradient_y_fft, CUFFT_INVERSE) );
    if (dimnum == 3){
        cufftHandleError( cufftExecZ2Z(plan_3d, gradient_z_fft, gradient_z_fft, CUFFT_INVERSE) );
    }

    if (dimnum == 2){
        l2_norm<<<grid, threads>>>(gradient_x_fft, gradient_y_fft, edges_gpu);
        cudaCheckError();
    }
    else if (dimnum == 3){
        l2_norm<<<grid, threads>>>(gradient_x_fft, gradient_y_fft, 
                                   gradient_z_fft, edges_gpu);
        cudaCheckError();
    }

    // Copying edges back
    cudaHandleError( cudaMemcpy(edges, edges_gpu, sizeof(double) * gSize, 
                                cudaMemcpyDeviceToHost) );

    // Method to find edges based on window approach -- more efficient, 
    // but difficult to implement

    par.store("density_d", density_d);
    par.store("density_d2", density_d2);
    par.store("edges_gpu", edges_gpu);
    par.store("gradient_x_fft", gradient_x_fft);
    par.store("gradient_y_fft", gradient_y_fft);
    par.store("gradient_z_fft", gradient_z_fft);
}
