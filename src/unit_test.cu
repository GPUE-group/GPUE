#include "../include/ds.h"
#include "../include/unit_test.h"
#include "../include/parser.h"
#include "../include/evolution.h"
#include "../include/init.h"
#include "../include/dynamic.h"
#include "../include/vortex_3d.h"
#include <string.h>
#include <assert.h>
#include <cufft.h>
#include <vector>
#include <fstream>
#include <sstream>

// Adding tests for mathematical operator kernels
void math_operator_test();
__global__ void add_test(double2 *a, double2 *b, double2 *c);
__global__ void subtract_test(double2 *a, double2 *b, double2 *c);
__global__ void mult_test(double2 *a, double2 *b, double2 *c);
__global__ void mult_test(double2 *a, double b, double2 *c);
__global__ void pow_test(double2 *a, int b, double2 *c);

// Tests for cufftDoubleComplex functions
void cufftDoubleComplex_functions_test();
__global__ void complexMag_test(double2 *in, double *out);
__global__ void complexMag2_test(double2 *in, double *out);
void realCompMult_test();
void cMult_test();

// Tests for quantum operations
__global__ void make_complex_kernel(double *in, int *evolution_type, 
                                    double2 *out);
void make_complex_test();
void cMultPhi_test();

// Tests for complex mathematical operations
void vecMult_test();
void scalarDiv_test();
void vecConj_test();

// AST tests
void ast_mult_test();
void ast_cmult_test();
void ast_op_mult_test();

// Other
void energyCalc_test();
void braKetMult_test();

// Test for the Grid structure with parameters in it 
void parameter_test();

// Test for the parsing function
void parser_test();

// Testing the evolve_2d function in evolution.cu
void evolve_test();

// Testing the parSum function
void parSum_test();

// Simple test of grid / cuda stuff
void grid_test2d();
void grid_test3d();

// Test of 1D fft's along all 3d grids
void fft_test();

// Test to check the equation parser for dynamic fields
void dynamic_test();

// Test to make sure the kernel for the polynomial approx. of Bessel fxns works
void bessel_test();

// Test for the vortex tracking functions in vortex_3d
void vortex3d_test();

// Test for available amount of GPU memory
void check_memory_test();

// Kernel testing will be added later
__device__ bool close(double a, double b, double threshold){
    return (abs(a-b) < threshold);
}


/*----------------------------------------------------------------------------//
* MAIN
*-----------------------------------------------------------------------------*/

void test_all(){
    std::cout << "Starting unit tests..." << '\n';
    parameter_test();

    std::cout 
        << "Beginning testing of standard mathematical operation kernels...\n";
    math_operator_test();

    std::cout << "Beginning testing of cufftDoubleComplex kernels...\n";
    cufftDoubleComplex_functions_test();

    // Do not uncomment these 2
    //parser_test();

    grid_test2d();
    grid_test3d();
    parSum_test();
    fft_test();
    dynamic_test();
    bessel_test();
    //vortex3d_test();
    make_complex_test();
    cMultPhi_test();
    evolve_test();

    check_memory_test();

    std::cout << "All tests completed. GPUE passed." << '\n';
}

void math_operator_test(){

    // First, we need to create a set of grids and threads to read into the 
    // kernels for testing
    dim3 grid = {1,1,1};
    dim3 threads = {1,1,1};

    double2 *ha, *hb, *hc;
    double2 *da, *db, *dc;
    
    // Allocating single-element arrays to test kernels with. 
    ha = (double2*)malloc(sizeof(double2));
    hb = (double2*)malloc(sizeof(double2));
    hc = (double2*)malloc(sizeof(double2));

    ha[0].x = 0.01;
    ha[0].y = 0.1;
    hb[0].x = 0.02;
    hb[0].y = 0.2;

    cudaMalloc((void**) &da, sizeof(double2));
    cudaMalloc((void**) &db, sizeof(double2));
    cudaMalloc((void**) &dc, sizeof(double2));

    cudaMemcpy(da, ha, sizeof(double2), cudaMemcpyHostToDevice);
    cudaMemcpy(db, hb, sizeof(double2), cudaMemcpyHostToDevice);

    add_test<<<grid, threads>>>(da, db, dc);
    cudaMemcpy(hc, dc, sizeof(double2), cudaMemcpyDeviceToHost);

    if (abs(hc[0].x - 0.03) > 1e-16 || abs(hc[0].y - 0.3) > 1e-16){
        std::cout << "Complex addition test failed!\n";
        exit(1);
    }

    subtract_test<<<grid, threads>>>(da, db, dc);
    cudaMemcpy(hc, dc, sizeof(double2), cudaMemcpyDeviceToHost);

    if (hc[0].x != -0.01 || hc[0].y != -0.1){
        std::cout << "Complex subtraction test failed!\n";
        exit(1);
    }

    pow_test<<<grid, threads>>>(da, 3, dc);
    cudaMemcpy(hc, dc, sizeof(double2), cudaMemcpyDeviceToHost);

    if (abs(hc[0].x + 0.000299) > 1e-16 || abs(hc[0].y + 0.00097) > 1e-16){
        std::cout << "Complex power test failed!\n";
        exit(1);
    }

    mult_test<<<grid, threads>>>(da, db, dc);
    cudaMemcpy(hc, dc, sizeof(double2), cudaMemcpyDeviceToHost);

    if (abs(hc[0].x + 0.0198) > 1e-16 || abs(hc[0].y - 0.004) > 1e-16){
        std::cout << "Complex multiplication test failed!\n";
        exit(1);
    }

    mult_test<<<grid, threads>>>(da, 3.0, dc);
    cudaMemcpy(hc, dc, sizeof(double2), cudaMemcpyDeviceToHost);

    if (abs(hc[0].x - 0.03) > 1e-16 || abs(hc[0].y - 0.3) > 1e-16){
        std::cout << "Complex multiplication test with real number failed!\n";
        exit(1);
    }

    std::cout << "Complex addition, subtraction, multiplication, and powers have been tested\n";

    cudaFree(da);
    cudaFree(db);
    cudaFree(dc);


    std::cout << "Now testing the derive() kernels...\n";
    // Now testing the derive function
    int dim = 4;
    double2 *darray, *darray_gpu, *darray_out;
    darray = (double2 *)malloc(sizeof(double2)*dim*dim*dim);
    cudaMalloc((void**) &darray_gpu, sizeof(double2)*dim*dim*dim);
    cudaMalloc((void**) &darray_out, sizeof(double2)*dim*dim*dim);

    for (int i = 0; i < dim; ++i){
        for (int j = 0; j < dim; ++j){
            for (int k = 0; k < dim; ++k){
                int index = k + j*dim + i*dim*dim;
                darray[index].x = i + j + k;
                darray[index].y = i + j + k;
            }
        }
    }

    cudaMemcpy(darray_gpu, darray, sizeof(double2)*dim*dim*dim,
               cudaMemcpyHostToDevice);

    grid = {1, dim, dim};
    threads = {dim, 1, 1};

    derive<<<grid, threads>>>(darray_gpu, darray_out, 1, dim*dim*dim,1);
    cudaMemcpy(darray, darray_out, sizeof(double2)*dim*dim*dim,
               cudaMemcpyDeviceToHost);
    for (int i = 0; i < dim-1; ++i){
        for (int j = 0; j < dim-1; ++j){
            for (int k = 0; k < dim-1; ++k){
                int index = k + j*dim + i*dim*dim;
                assert(darray[index].x == 1);
                assert(darray[index].y == 1);
            }
        }
    }

    derive<<<grid, threads>>>(darray_gpu, darray_out, dim, dim*dim*dim,1);
    cudaMemcpy(darray, darray_out, sizeof(double2)*dim*dim*dim,
               cudaMemcpyDeviceToHost);
    for (int i = 0; i < dim-1; ++i){
        for (int j = 0; j < dim-1; ++j){
            for (int k = 0; k < dim-1; ++k){
                int index = k + j*dim + i*dim*dim;
                assert(darray[index].x == 1);
                assert(darray[index].y == 1);
            }
        }
    }

    derive<<<grid, threads>>>(darray_gpu, darray_out, dim*dim, dim*dim*dim,1);
    cudaMemcpy(darray, darray_out, sizeof(double2)*dim*dim*dim,
               cudaMemcpyDeviceToHost);
    for (int i = 0; i < dim-1; ++i){
        for (int j = 0; j < dim-1; ++j){
            for (int k = 0; k < dim-1; ++k){
                int index = k + j*dim + i*dim*dim;
                assert(darray[index].x == 1);
                assert(darray[index].y == 1);
            }
        }
    }


    std::cout << "derive functions passed!\n";
}

__global__ void add_test(double2 *a, double2 *b, double2 *c){
    c[0] = add(a[0],b[0]);
}

__global__ void subtract_test(double2 *a, double2 *b, double2 *c){
    c[0] = subtract(a[0],b[0]);
}

__global__ void pow_test(double2 *a, int b, double2 *c){
    c[0] = pow(a[0],b);
}

__global__ void mult_test(double2 *a, double2 *b, double2 *c){
    c[0] = mult(a[0],b[0]);
}

__global__ void mult_test(double2 *a, double b, double2 *c){
    c[0] = mult(a[0],b);
}

void cufftDoubleComplex_functions_test(){

    // first creating the grid and threads
    dim3 grid = {1,1,1};
    dim3 threads = {1,1,1};

    double *hval_double, *dval_double, *dout, *hout;
    double2 *hval_double2, *dval_double2;

    double2 *hin, *din;

    hval_double = (double*)malloc(sizeof(double));
    hval_double2 = (double2*)malloc(sizeof(double2));
    hout = (double*)malloc(sizeof(double));
    hin = (double2*)malloc(sizeof(double2));

    hval_double[0] = 3.0;
    hval_double2[0].x = 0.3;
    hval_double2[0].y = 0.4;

    hin[0].x = 3.0;
    hin[0].y = 4.0;

    cudaMalloc((void**)&dval_double, sizeof(double));
    cudaMalloc((void**)&dval_double2, sizeof(double2));
    cudaMalloc((void**)&dout, sizeof(double));
    cudaMalloc((void**)&din, sizeof(double2));


    // Testing make_cufftDoubleComplex function
    cudaMemcpy(dval_double, hval_double, sizeof(double), 
               cudaMemcpyHostToDevice);

    make_cufftDoubleComplex<<<grid, threads>>>(dval_double, dval_double2);

    cudaMemcpy(hval_double2, dval_double2, sizeof(double2),
               cudaMemcpyDeviceToHost);

    if (hval_double2[0].x != 3.0 || hval_double2[0].y != 0){
        std::cout << "Test of make_cufftDoubleComplex failed!\n";
        exit(1);
    }

    // testing device complexMagnitude function
    cudaMemcpy(din, hin, sizeof(double2), cudaMemcpyHostToDevice);
    complexMag_test<<<grid, threads>>>(din, dout);

    cudaMemcpy(hout, dout, sizeof(double), cudaMemcpyDeviceToHost);

    if (hout[0] != 5.0){
        std::cout << hout[0] << '\n';
        std::cout << "Test of device complexMagnitude failed!\n";
        exit(1);
    }

    // Testing global complexMagnitude function
    complexMagnitude<<<grid, threads>>>(din, dout);
    cudaMemcpy(hout, dout, sizeof(double), cudaMemcpyDeviceToHost);

    if (hout[0] != 5.0){
        std::cout << hout[0] << '\n';
        std::cout << "Test of global complexMagnitude failed!\n";
        exit(1);
    }

    complexMag2_test<<<grid, threads>>>(din, dout);

    cudaMemcpy(hout, dout, sizeof(double), cudaMemcpyDeviceToHost);

    if (hout[0] != 25.0){
        std::cout << hout[0] << '\n';
        std::cout << "Test of device complexMagnitudeSquared failed!\n";
        exit(1);
    }

    // Testing global complexMagnitude function
    complexMagnitudeSquared<<<grid, threads>>>(din, dout);
    cudaMemcpy(hout, dout, sizeof(double), cudaMemcpyDeviceToHost);

    if (hout[0] != 25.0){
        std::cout << hout[0] << '\n';
        std::cout << "Test of global complexMagnitudeSquared failed!\n";
        exit(1);
    }


    std::cout << "make_cufftDoubleComplex, and complexMagnitude[Squared] have been tested\n";

}

__global__ void complexMag_test(double2 *in, double *out){
    out[0] = complexMagnitude(in[0]);
}

__global__ void complexMag2_test(double2 *in, double *out){
    out[0] = complexMagnitudeSquared(in[0]);
}

// Test to check the equation parser for dynamic fields
// For this test, we will need a general set of parameters to read in and a
// standard equation string to look at. 
void dynamic_test(){

    std::cout << "Beginning test of dynamic functions..." <<'\n';
    std::string eqn_string = "(((3*x)+7)+(5-7)+cos(0)*1)+pow(120,2)";

    Grid par;
    par.store("x",5.0);
    par.store("y",5.0);
    par.store("z",5.0);
    par.store("omegaX",5.0);
    par.store("omegaY",5.0);
    par.store("omegaZ",5.0);
    par.store("xMax",5.0);
    par.store("yMax",5.0);
    par.store("zMax",5.0);
    par.store("fudge",5.0);
    par.store("mass",5.0);
    std::string val_string = "check_var";

    EqnNode eqn_tree = parse_eqn(par, eqn_string, val_string);

    std::cout << "finding the number of elements in abstract syntax tree...\n";

    int num = 0;
    find_element_num(eqn_tree, num);
    int element_num = num;

    std::cout << "Total number of elements is: " << num << '\n';

    std::cout << "Now to copy the tree to the GPU..." << '\n';

    EqnNode_gpu *eqn_gpu, *eqn_cpu;
    eqn_cpu = (EqnNode_gpu *)malloc(sizeof(EqnNode_gpu)*element_num);
    num = 0;
    tree_to_array(eqn_tree, eqn_cpu, num);

/*
    for (int i = 0; i < num; ++i){
        std::cout << eqn_cpu[i].val << '\n';
        std::cout << eqn_cpu[i].left << '\n';
        std::cout << eqn_cpu[i].right << '\n' << '\n';
    }
*/

    cudaMalloc((void**)&eqn_gpu, sizeof(EqnNode_gpu)*element_num);
    cudaMemcpy(eqn_gpu, eqn_cpu, sizeof(EqnNode_gpu)*element_num,
               cudaMemcpyHostToDevice);

    // Now to check some simple evaluation
    std::cout << "Now to check simple GPU evaluation..." << '\n';
    int n = 64;
    double *array, *array_gpu;
    array = (double *)malloc(sizeof(double)*n);
    cudaMalloc(&array_gpu, sizeof(double)*n);

    int threads = 64;
    int grid = (int)ceil((float)n/threads);

    //zeros<<<grid, threads>>>(array_gpu, n);
    find_field<<<grid, threads>>>(array_gpu, 1, 0.0, 0.0, 0.0, 1,1,1,eqn_gpu);

    cudaMemcpy(array, array_gpu, sizeof(double)*n, cudaMemcpyDeviceToHost);

    for (int i = 0; i < n; ++i){
        double eqn_val = (((3*i)+7)+(5-7)+cos(0)*1)+pow(120,2);
        if (array[i] != eqn_val){
            std::cout << "GPU evaluation failed in dynamic test!\n";
            assert(array[i] == eqn_val);
        }
    }

    // Now testing simple parsing of example "example.cfg"
#ifdef _CFG_FILE_PATH
    std::stringstream ss;
    ss << _CFG_FILE_PATH;
    std::string cfg_path;
    ss >> cfg_path;
#else
    std::string cfg_path = "src/example.cfg";
#endif

    std::cout << "Testing simple parameter parsing." << '\n';
    par.store("param_file", cfg_path);
    parse_param_file(par);
    EqnNode_gpu *eqn = par.astval("V");
    find_field<<<grid, threads>>>(array_gpu, 0.1, 0.1, 0.1, 1, 1, 1, 0, eqn);
    cudaDeviceSynchronize();
    
    cudaMemcpy(array, array_gpu, sizeof(double)*n, cudaMemcpyDeviceToHost);

/*
    for (int i = 0; i < n; ++i){
        std::cout << array[i] << '\n';
    }
*/

    std::cout << "Dynamic tests passed" <<'\n';
}

__global__ void bessel_test_kernel(double *j, double *j_poly, bool *val){
    int xid = blockDim.x*blockIdx.x + threadIdx.x;
    j[xid] = j0(xid * 2.0 / 128);
    j_poly[xid] = poly_j(0,xid * 2.0 / 128, 40);

    if (!close(j[xid],j_poly[xid], 0.0001)){
        val[0] = false;
        printf("Error at element %u in Bessel test!\tValues: %f, %f\n", 
               xid, j[xid], j_poly[xid]);
    }
}

// Test for bessel functions
void bessel_test(){

    std::cout << "Testing Bessel Functions..." << '\n';

    double *j_gpu, *j_poly_gpu;
    bool *val, *val_gpu;
    int n = 128;
    cudaMalloc((void **)&j_gpu, sizeof(double)*n);
    cudaMalloc((void **)&j_poly_gpu, sizeof(double)*n);

    cudaMalloc((void **)&val_gpu, sizeof(bool));
    val = (bool *)malloc(sizeof(bool));
    val[0] = true;
    cudaMemcpy(val_gpu, val, sizeof(bool), cudaMemcpyHostToDevice);

    bessel_test_kernel<<<64,2>>>(j_gpu, j_poly_gpu, val_gpu);
    cudaMemcpy(val, val_gpu, sizeof(bool), cudaMemcpyDeviceToHost);

    if(val[0]){
        std::cout << "Bessel Test Passed!" << '\n';
    }
    else{
        std::cout << "Bessel Test Failed!" << '\n';
        exit(1);
    }

}

// Test of 1D fft's along all 3d grids
// In particular, we need to test the generate_plan_other3d function
// These will be checked against 1d 
void fft_test(){

    std::cout << "Beginning cufft test.\n";

    // For these tests, we are assuming that the x, y and z dimensions are 
    // All the same (2x2x2)
    // Note that yDim needs to be singled out differently, but z/x need no loops

    // now we need to create the necessary parameters and store everything
    int xDim = 2;
    int yDim = 2;
    int zDim = 2;
    int gsize = xDim * yDim * zDim;

    Grid par;
    par.store("xDim", xDim);
    par.store("yDim", yDim);
    par.store("zDim", zDim);

    cufftHandle plan_x, plan_y, plan_z;
    // Now creating the plans
    generate_plan_other3d(&plan_x, par, 0);
    generate_plan_other3d(&plan_y, par, 1);
    generate_plan_other3d(&plan_z, par, 2);

    // And the result / error
    cudaError_t err;
    cufftResult result;

    // Creating the initial array for the x dimension fft
    double2 *array, *gpu_array;
    array = (double2 *) malloc(sizeof(double2)*gsize);
    cudaMalloc((void**) &gpu_array, sizeof(double2)*gsize);
    for (int i = 0; i < gsize; i++){
        array[i].x = 1;
        array[i].y = 0;
    }

    // transferring to gpu
    err = cudaMemcpy(gpu_array, array, sizeof(double2)*gsize,
                     cudaMemcpyHostToDevice);
    if (err != cudaSuccess){
        std::cout << "Could not coppy array to device!" << '\n';
        std::cout << "error code: " << err << '\n';
        exit(1);
    }

    // Performing the x transformation
    for (int i = 0; i < yDim; i++){
        result = cufftExecZ2Z(plan_y, &gpu_array[i*xDim*yDim], 
                                      &gpu_array[i*xDim*yDim], CUFFT_FORWARD);
    }
    //result = cufftExecZ2Z(plan_z, gpu_array, gpu_array, CUFFT_FORWARD);

    // transferring back to host to check output
    err = cudaMemcpy(array, gpu_array, sizeof(double2)*gsize, 
                     cudaMemcpyDeviceToHost);
    if (err != cudaSuccess){
        std::cout << "Could not coppy gpu_array to host!" << '\n';
        std::cout << "error code: " << err << '\n';
        exit(1);
    }

/*
    for (int i = 0; i < gsize; i++){
        std::cout << array[i].x << '\t' << array[i].y << '\n';
    }
*/

    // Now to try the inverse direction

    for (int i = 0; i < yDim; i++){
        result = cufftExecZ2Z(plan_y, &gpu_array[i*xDim*yDim], 
                                      &gpu_array[i*xDim*yDim], CUFFT_INVERSE);
    }
    //result = cufftExecZ2Z(plan_z, gpu_array, gpu_array, CUFFT_INVERSE);

    // copying back
    err = cudaMemcpy(array, gpu_array, sizeof(double2)*gsize, 
                     cudaMemcpyDeviceToHost);
    if (err != cudaSuccess){
        std::cout << "Could not coppy gpu_array to host!" << '\n';
        std::cout << "error code: " << err << '\n';
        exit(1);
    }

/*
    for (int i = 0; i < gsize; i++){
        std::cout << array[i].x << '\t' << array[i].y << '\n';
    }
*/

    std::cout << "cufft test passed!\n";



}

// Simple test of CUDA grid stuff
void grid_test2d(){

    std::cout << "testing grid / threads and stuff" << '\n';

    int max_threads = 256;

    int xDim = 1024;
    int yDim = 1024;
    int zDim = 1;

    int xD = 1, yD = 1, zD = 1;

    int gsize = xDim * yDim;

    // Now to set up the CUDA grid / threads
    dim3 block;
    dim3 grid;

    Grid par;
    par.store("xDim",xDim);
    par.store("yDim",yDim);
    par.store("zDim",zDim);
    par.store("dimnum", 2);

    generate_grid(par);

    if (xDim <= max_threads){
        block.x = xDim;
        block.y = 1;
        block.z = 1;

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

        block.x = dim_tmp;
        block.y = 1;
        block.z = 1;
        xD = pow(2,count);
        yD = yDim;
        zD = 1;
    }

    std::cout << "threads in x are: " << block.x << '\n';
    std::cout << "dimensions are: " << xD << '\t' << yD << '\t' << zD << '\n';

    grid.x=xD; 
    grid.y=yD; 
    grid.z=zD; 

    if (grid.x != par.grid.x || block.x != par.threads.x ||
        grid.y != par.grid.y || block.y != par.threads.y ||
        grid.z != par.grid.z || block.z != par.threads.z){
        std::cout << "Gridding test 2D failed! Improper generation of threads!"
                  << '\n';
        assert(grid.x == par.grid.x);
        assert(grid.y == par.grid.y);
        assert(grid.z == par.grid.z);
        assert(block.x == par.threads.x);
        assert(block.y == par.threads.y);
        assert(block.z == par.threads.z);
    }

    int total_threads = block.x * block.y * block.z;

    // Now we need to initialize our double * and send it to the gpu
    double *host_array, *device_array;
    host_array = (double *) malloc(sizeof(double)*gsize);
    cudaMalloc((void**) &device_array, sizeof(double)*gsize);

    // initializing 2d array
    for (int i = 0; i < gsize; i++){
        host_array[i] = -1;
    }

    // Now to copy to device
    cudaMemcpy(device_array, host_array,
               sizeof(double)*gsize,
               cudaMemcpyHostToDevice);

    // Test
    thread_test<<<grid,block>>>(device_array,device_array);

    // Now to copy back and print
    cudaMemcpy(host_array, device_array,
               sizeof(double)*gsize,
               cudaMemcpyDeviceToHost);
    
    
    for (int i = 0; i < xDim; ++i){
        for (int j = 0; j < yDim; ++j){
            int index = i*yDim + j;
            if (host_array[index] != index){
                std::cout << "Threadding values improperly set!\n";
                assert(host_array[index] == index);
            }
        }
    }

    std::cout << "2d grid tests completed. Now for 3d cases" << '\n';

}

// Simple test of CUDA grid stuff
void grid_test3d(){

    std::cout << "testing grid / threads and stuff for 3d" << '\n';

    int max_threads = 256;

    int xDim = 256;
    int yDim = 256;
    int zDim = 256;

    int xD = 1, yD = 1, zD = 1;

    int gsize = xDim * yDim * zDim;

    // Now to set up the CUDA grid / threads
    dim3 block;
    dim3 grid;

    Grid par;
    par.store("xDim",xDim);
    par.store("yDim",yDim);
    par.store("zDim",zDim);
    par.store("dimnum", 3);

    generate_grid(par);


    if (xDim <= max_threads){
        block.x = xDim;
        block.y = 1;
        block.z = 1;

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

        block.x = dim_tmp;
        block.y = 1;
        block.z = 1;
        xD = pow(2,count);
        yD = yDim;
        zD = zDim;
    }

    std::cout << "threads in x are: " << block.x << '\n';
    std::cout << "dimensions are: " << xD << '\t' << yD << '\t' << zD << '\n';

    grid.x=xD; 
    grid.y=yD; 
    grid.z=zD; 

    if (grid.x != par.grid.x || block.x != par.threads.x ||
        grid.y != par.grid.y || block.y != par.threads.y ||
        grid.z != par.grid.z || block.z != par.threads.z){
        std::cout << "Gridding test 3D failed! Improper generation of threads!"
                  << '\n';
        assert(grid.x == par.grid.x);
        assert(grid.y == par.grid.y);
        assert(grid.z == par.grid.z);
        assert(block.x == par.threads.x);
        assert(block.y == par.threads.y);
        assert(block.z == par.threads.z);
    }


    int total_threads = block.x * block.y * block.z;

    // Now we need to initialize our double * and send it to the gpu
    double *host_array, *device_array;
    host_array = (double *) malloc(sizeof(double)*gsize);
    cudaMalloc((void**) &device_array, sizeof(double)*gsize);

    // initializing 2d array
    for (int i = 0; i < gsize; i++){
        host_array[i] = -1;
    }

    // Now to copy to device
    cudaMemcpy(device_array, host_array,
               sizeof(double)*gsize,
               cudaMemcpyHostToDevice);

    // Test
    thread_test<<<grid,block>>>(device_array,device_array);

    // Now to copy back and print
    cudaMemcpy(host_array, device_array,
               sizeof(double)*gsize,
               cudaMemcpyDeviceToHost);
    
    
    for (int i = 0; i < xDim; ++i){
        for (int j = 0; j < yDim; ++j){
            for (int k = 0; k < zDim; ++k){
                int index = i*yDim*zDim + j*yDim + k;
                if (host_array[index] != index){
                    std::cout << "Threadding values improperly set!\n";
                    assert(host_array[index] == index);
                }
            }
        }
    }

    std::cout << "3d grid tests completed. now for 3d cases" << '\n';

}

// Test of the parSum function in 3d
void parSum_test(){

    std::cout << "Beginning test of parallel summation.\n";

    // Setting error
    cudaError_t err;

    // first, we need to initialize the Grid and Cuda classes
    Grid par;

    // 2D test first

    // For now, we will assume an 64x64 array for summing
    dim3 threads(16, 1, 1);
    int total_threads = threads.x*threads.y*threads.z;

    double dx = 0.1;
    double dy = 0.1;
    double dz = 0.1;

    par.store("dimnum", 2);
    par.store("xDim", 64);
    par.store("yDim", 64);
    par.store("zDim", 1);
    par.store("dx",dx);
    par.store("dy",dy);
    par.store("dz",dz);
    par.threads = threads;

    // Now we need to initialize the grid for the getGid3d3d kernel
    int gsize = 64*64;
    dim3 grid;
    grid.x = 4;
    grid.y = 64;

    par.grid = grid;

    // now we need to initialize the wfc to all 1's;
    double2 *wfc, *host_sum;
    wfc = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * gsize);
    host_sum = (cufftDoubleComplex *) 
               malloc(sizeof(cufftDoubleComplex) * gsize / total_threads);

    // init wfc
    for (int i = 0; i < gsize; i++){
        wfc[i].x = 2;
        wfc[i].y = 2;
    }

    double2 *gpu_wfc;
    cudaMalloc((void**) &gpu_wfc, sizeof(cufftDoubleComplex)*gsize);

    // copying wfc to device
    err = cudaMemcpy(gpu_wfc, wfc, sizeof(cufftDoubleComplex)*gsize,
                     cudaMemcpyHostToDevice);

    if (err!=cudaSuccess){
        std::cout << "ERROR: Could not copy wfc to device!" << '\n';
    }

    // Creating parsum on device

    parSum(gpu_wfc, par);

    // copying parsum back
    err = cudaMemcpy(wfc, gpu_wfc, 
                     sizeof(cufftDoubleComplex)*gsize, 
                     cudaMemcpyDeviceToHost);
    if (err!=cudaSuccess){
        std::cout << err << '\n';
        std::cout << "ERROR: Could not copy wfc to the host!" << '\n';
        exit(1);
    }

    for (int i = 0; i < gsize; ++i){
        if (wfc[i].x != 2/sqrt(32768.0*dx*dy) ||
            wfc[i].y != 2/sqrt(32768.0*dx*dy)){
            std::cout << "Wavefunction not normalized!" << '\n';
            std::cout << i << '\t' << wfc[i].x << '\t' << wfc[i].y << '\n';
            assert(wfc[i].x == 2/sqrt(32768.0*dx*dy));
            assert(wfc[i].y == 2/sqrt(32768.0*dx*dy));
        }
    }

    // Now for the 3d case
    // For now, we will assume a 16x16x16 array for summing
    par.store("dimnum", 3);
    par.store("xDim", 16);
    par.store("yDim", 16);
    par.store("zDim", 16);
    par.store("dx",dx);
    par.store("dy",dy);
    par.store("dz",dz);

    // Now we need to initialize the grid for the getGid3d3d kernel
    grid.x = 16;
    grid.y = 16;
    grid.z = 16;

    par.grid = grid;

    // copying host wfc back to device
    err = cudaMemcpy(gpu_wfc, wfc, sizeof(cufftDoubleComplex)*gsize,
                     cudaMemcpyHostToDevice);

    parSum(gpu_wfc, par);

    // copying parsum back
    err = cudaMemcpy(wfc, gpu_wfc, 
                     sizeof(cufftDoubleComplex)*gsize, 
                     cudaMemcpyDeviceToHost);
    if (err!=cudaSuccess){
        std::cout << "ERROR: Could not copy wfc to the host!" << '\n';
        exit(1);
    }

    for (int i = 0; i < gsize; ++i){
        if (wfc[i].x != 2/sqrt(32768.0*dx*dy*dz) ||
            wfc[i].y != 2/sqrt(32768.0*dx*dy*dz)){
            std::cout << "Wavefunction not normalized!" << '\n';
            std::cout << wfc[i].x << '\t' << wfc[i].y << '\n';
            assert(wfc[i].x == 2/sqrt(32768.0*dx*dy*dz));
            assert(wfc[i].y == 2/sqrt(32768.0*dx*dy*dz));
        }
    }

    std::cout << "Parallel summation test passed in 2 and 3D!\n";

}

// Test for the Grid structure with paramters in it
// Initialize all necessary variables and read them back out
void parameter_test(){
    // For this test, we simply need to read in and out stuff from each 
    // class and structure in ds.h / ds.cc
    
    // Certain variables will be used multiple times. 
    double *dstar_var;
    dstar_var = (double *)malloc(sizeof(double) * 5);
    cufftDoubleComplex *cdc_var;
    cdc_var = (cufftDoubleComplex *)malloc(sizeof(cufftDoubleComplex) * 5);
    for (int i = 0; i < 5; ++i){
        dstar_var[i] = (double)i * 0.5;
        cdc_var[i].x = (double)i * 0.5;
        cdc_var[i].y = (double)i * 0.5;
    }

    double dvar = 1.05;
    int ivar = 5;
    bool bvar = true;

    // Now testing the Grid class
    Grid grid_test;
    grid_test.store("dstar_var",dstar_var);
    grid_test.store("dvar", dvar);
    grid_test.store("ivar", ivar);
    grid_test.store("bvar", bvar);

    assert(dstar_var == grid_test.dsval("dstar_var"));
    assert(dvar == grid_test.dval("dvar"));
    assert(ivar == grid_test.ival("ivar"));
    assert(bvar == grid_test.bval("bvar"));

    std::cout << "Grid class checked, now checking the Cuda class..." << '\n';

    std::cout << "All data structures checked" << '\n';

}

// Test for the parsing function
void parser_test(){

    // Testing the command-line parser with defaults and with modifications
    std::cout << "Testing command-line parser with no arguments..." << '\n';

    // First testing default values in and out of the parser function
    char *fake_noargv[] = {NULL};
    Grid noarg_grid;
    noarg_grid = parseArgs(0,fake_noargv);

    // Checking contents of noarg_grid:
    assert(noarg_grid.ival("xDim") == 256);
    assert(noarg_grid.ival("yDim") == 256);
    assert(noarg_grid.ival("zDim") == 256);
    assert(noarg_grid.dval("omega") == 0);
    assert(noarg_grid.dval("gammaY") == 1.0);
    assert(noarg_grid.ival("gsteps") == 1);
    assert(noarg_grid.ival("esteps") == 1);
    assert(noarg_grid.dval("gdt") == 1e-4);
    assert(noarg_grid.dval("dt") == 1e-4);
    assert(noarg_grid.ival("device") == 0);
    assert(noarg_grid.ival("atoms") == 1);
    assert(noarg_grid.bval("read_wfc") == false);
    assert(noarg_grid.ival("printSteps") == 100);
    assert(noarg_grid.dval("winding") == 0);
    assert(noarg_grid.bval("corotating") == false);
    assert(noarg_grid.bval("gpe") == false);
    assert(noarg_grid.dval("omegaZ") == 6.283);
    assert(noarg_grid.dval("interaction") == 1);
    assert(noarg_grid.dval("laser_power") == 0);
    assert(noarg_grid.dval("angle_sweep") == 0);
    assert(noarg_grid.ival("kick_it") == 0);
    assert(noarg_grid.bval("write_it") == false);
    assert(noarg_grid.dval("x0_shift") == 0);
    assert(noarg_grid.dval("y0_shift") == 0);
    assert(noarg_grid.dval("z0_shift") == 0);
    assert(noarg_grid.dval("sepMinEpsilon") == 0);
    assert(noarg_grid.bval("graph") == false);
    assert(noarg_grid.bval("unit_test") == false);
    assert(noarg_grid.dval("omegaX") == 6.283);
    assert(noarg_grid.dval("omegaY") == 6.283);
    assert(noarg_grid.sval("data_dir") == "data/");
    assert(noarg_grid.bval("ramp") == false);
    assert(noarg_grid.ival("ramp_type") == 1);
    assert(noarg_grid.ival("dimnum") == 2);
    assert(noarg_grid.bval("write_file") == true);
    assert(noarg_grid.dval("fudge") == 1.0);
    assert(noarg_grid.ival("kill_idx") == -1);
    assert(noarg_grid.dval("mask_2d") == 1.5e-4);
    assert(noarg_grid.dval("box_size") == 2.5e-5);
    assert(noarg_grid.bval("found_sobel") == false);
    assert(noarg_grid.Afn == "rotation");
    assert(noarg_grid.Kfn == "rotation_K");
    assert(noarg_grid.Vfn == "2d");
    assert(noarg_grid.Wfcfn == "2d");
    assert(noarg_grid.sval("conv_type") == "FFT");
    assert(noarg_grid.ival("charge") == 0);
    assert(noarg_grid.bval("flip") == false);

    // Now testing all values specified by command-line arguments
    std::cout << "Testing command-line parser with all arguments..." << '\n';
    std::vector<std::string> argarray(10);

    // I apologize for the mess... If you have a better way of creating the 
    // char ** for this without running into memory issues, let me know!
    char *fake_fullargv[] = {strdup("./gpue"), 
                             strdup("-A"), strdup("rotation"), 
                             strdup("-a"),
                             strdup("-b"), strdup("2.5e-5"), 
                             strdup("-C"), strdup("0"), 
                             strdup("-c"), strdup("3"), 
                             strdup("-D"), strdup("data"), 
                             strdup("-E"), 
                             strdup("-e"), strdup("1"), 
                             strdup("-f"), 
                             strdup("-G"), strdup("1"),
                             strdup("-g"), strdup("1"), 
                             strdup("-i"), strdup("1"), 
                             strdup("-K"), strdup("0"), 
                             strdup("-k"), strdup("0"),
                             strdup("-L"), strdup("0"), 
                             strdup("-l"), 
                             strdup("-n"), strdup("1"), 
                             strdup("-O"), strdup("0"),
                             strdup("-P"), strdup("0"), 
                             strdup("-p"), strdup("100"),
                             strdup("-Q"), strdup("0"), 
                             strdup("-q"), strdup("0"), 
                             strdup("-R"), strdup("1"), 
                             //strdup("-r"),
                             strdup("-S"), strdup("0"), 
                             strdup("-s"),
                             strdup("-T"), strdup("1e-4"), 
                             strdup("-t"), strdup("1e-4"), 
                             strdup("-U"), strdup("0"), 
                             strdup("-V"), strdup("0"), 
                             strdup("-W"), 
                             strdup("-w"), strdup("0"), 
                             strdup("-X"), strdup("1.0"),
                             strdup("-x"), strdup("256"), 
                             strdup("-Y"), strdup("1.0"), 
                             strdup("-y"), strdup("256"),
                             strdup("-Z"), strdup("6.283"), 
                             strdup("-z"), strdup("256"), 
                             NULL};
    int fake_argc = sizeof(fake_fullargv) / sizeof(char *) - 1;

    // Now to read into gpue and see what happens
    Grid fullarg_grid;
    fullarg_grid = parseArgs(fake_argc, fake_fullargv);

    // Checking contents of fullarg_grid:
    assert(fullarg_grid.ival("xDim") == 256);
    assert(fullarg_grid.ival("yDim") == 256);
    assert(fullarg_grid.ival("zDim") == 256);
    assert(fullarg_grid.dval("omega") == 0);
    assert(fullarg_grid.dval("gammaY") == 1.0);
    assert(fullarg_grid.ival("gsteps") == 1);
    assert(fullarg_grid.ival("esteps") == 1);
    assert(fullarg_grid.dval("gdt") == 1e-4);
    assert(fullarg_grid.dval("dt") == 1e-4);
    assert(fullarg_grid.ival("device") == 0);
    assert(fullarg_grid.ival("atoms") == 1);
    assert(fullarg_grid.bval("read_wfc") == false);
    assert(fullarg_grid.ival("printSteps") == 100);
    assert(fullarg_grid.dval("winding") == 0);
    assert(fullarg_grid.bval("corotating") == true);
    assert(fullarg_grid.bval("gpe") == true);
    assert(fullarg_grid.dval("omegaZ") == 6.283);
    assert(fullarg_grid.dval("interaction") == 1);
    assert(fullarg_grid.dval("laser_power") == 0);
    assert(fullarg_grid.dval("angle_sweep") == 0);
    assert(fullarg_grid.ival("kick_it") == 0);
    assert(fullarg_grid.bval("write_it") == true);
    assert(fullarg_grid.dval("x0_shift") == 0);
    assert(fullarg_grid.dval("y0_shift") == 0);
    assert(fullarg_grid.dval("z0_shift") == 0);
    assert(fullarg_grid.dval("sepMinEpsilon") == 0);
    assert(fullarg_grid.bval("graph") == true);
    assert(fullarg_grid.bval("unit_test") == false);
    assert(fullarg_grid.dval("omegaX") == 1.0);
    assert(fullarg_grid.dval("omegaY") == 1.0);
    assert(fullarg_grid.sval("data_dir") == "data/");
    assert(fullarg_grid.bval("ramp") == true);
    assert(fullarg_grid.ival("ramp_type") == 1);
    assert(fullarg_grid.ival("dimnum") == 3);
    assert(fullarg_grid.bval("write_file") == false);
    assert(fullarg_grid.dval("fudge") == 1.0);
    assert(fullarg_grid.ival("kill_idx") == 0);
    assert(fullarg_grid.dval("mask_2d") == 1.5e-4);
    assert(fullarg_grid.dval("box_size") == 2.5e-5);
    assert(fullarg_grid.bval("found_sobel") == false);
    assert(fullarg_grid.Afn == "rotation");
    assert(fullarg_grid.Kfn == "rotation_K3d");
    assert(fullarg_grid.Vfn == "3d");
    assert(fullarg_grid.Wfcfn == "3d");
    assert(fullarg_grid.sval("conv_type") == "FFT");
    assert(fullarg_grid.ival("charge") == 0);
    assert(fullarg_grid.bval("flip") == true);

}

// Testing the evolve function in evolution.cu
// This test will also test the energy function
// Run the simulation in imaginary time for a simple harmonic oscillator and 
//     check energies in 1, 2, and 3D.
// Run the simulation in real time to make sure the energy remains the same
void evolve_test(){

    std::cout << "Starting test of evolution function in nD...\n";

    // Setting default values
    Grid par;

    int res = 32;
    par.store("omega", 0.0);
    par.store("gammaY", 1.0);
    par.store("device", 0);
    par.store("read_wfc", false);
    par.store("winding", 0.0);
    par.store("corotating", false);
    par.store("gpe", false);
    par.store("interaction",1.0);
    par.store("laser_power",0.0);
    par.store("angle_sweep",0.0);
    par.store("kick_it", 0);
    par.store("x0_shift",0.0);
    par.store("y0_shift",0.0);
    par.store("z0_shift",0.0);
    par.store("sepMinEpsilon",0.0);
    par.store("graph", false);
    par.store("unit_test",false);
    par.store("ramp", false);
    par.store("ramp_type", 1);
    par.store("dimnum", 2);
    par.store("fudge", 0.0);
    par.store("kill_idx", -1);
    par.store("mask_2d", 1.5e-4);
    par.store("found_sobel", false);
    par.store("use_param_file", false);
    par.store("param_file","param.cfg");
    par.store("data_dir", (std::string)"data/");
    par.store("cyl_coord",false);
    par.Afn = "rotation";
    par.Kfn = "rotation_K";
    par.Vfn = "2d";
    par.Wfcfn = "2d";
    par.store("conv_type", (std::string)"FFT");
    par.store("charge", 0);
    par.store("flip", false);
    par.store("thresh_const", 1.0);


    double thresh = 0.01;
    std::string buffer;
    int gsteps = 30001;
    int esteps = 30001;

    par.store("gdt", 1e-4);
    par.store("dt", 1e-4);
    par.store("atoms", 1);
    par.store("omegaZ", 1.0);
    par.store("omegaX", 1.0);
    par.store("omegaY", 1.0);
    par.store("esteps", esteps);
    par.store("gsteps", gsteps);
    par.store("printSteps", 30000);
    par.store("write_file", false);
    par.store("write_it", false);
    par.store("energy_calc", true);
    par.store("corotating", true);
    par.store("omega",0.0);
    par.store("box_size", 0.00007);
    par.store("xDim", res);
    par.store("yDim", 1);
    par.store("zDim", 1);



    // Running through all the dimensions to check the energy
    for (int i = 1; i <= 3; ++i){
        if (i == 2){
            par.store("yDim", res);
        }
        if (i == 3){
            par.store("zDim", res);
        }
        par.store("dimnum",i);
        init(par);

        if (par.bval("write_file")){
            FileIO::writeOutParam(buffer, par, "data/Params.dat");
        }

        double omegaX = par.dval("omegaX");
        set_variables(par, 0);

        // Evolve and find the energy
        evolve(par, gsteps, 0, buffer);

        // Check that the energy is correct
        double energy = par.dval("energy");
        double energy_check = 0;
        energy_check = (double)i * 0.5 * HBAR * omegaX;

        if (abs(energy - energy_check) > thresh*energy_check){
            std::cout << "Energy is not correct in imaginary-time for " 
                      << i << "D!\n";
            assert(energy == energy_check);
        }

        // Run in real time to make sure that the energy is constant
        set_variables(par, 1);
        evolve(par, esteps, 1, buffer);
        double energy_ev = par.dval("energy");

        if (abs(energy - energy_ev) > thresh*energy_check){
            std::cout << "Energy is not constant in real-time for " 
                      << i << "D!\n";
            assert(energy == energy_ev);
        }
    }
    
}

__global__ void make_complex_kernel(double *in, int *evolution_type, 
                                    double2 *out){

    //int id = threadIdx.x + blockIdx.x*blockDim.x;
    //out[id] = make_complex(in[id], evolution_type[id]);
    for (int i = 0; i < 3; ++i){
        out[i] = make_complex(in[i], evolution_type[i]);
    }
}

void make_complex_test(){

    // Creating a simple array to hold the 3 possible make_complex options
    double *input_array, *dinput_array;
    double2 *output_array, *doutput_array;
    int *evolution_type, *devolution_type;

    input_array = (double *)malloc(sizeof(double)*3);
    output_array = (double2 *)malloc(sizeof(double2)*3);
    evolution_type = (int *)malloc(sizeof(int)*3);

    input_array[0] = 10;
    input_array[1] = 10;
    input_array[2] = 10;

    evolution_type[0] = 0;
    evolution_type[1] = 1;
    evolution_type[2] = 2;

    cudaMalloc((void **)&dinput_array, sizeof(double)*3);
    cudaMalloc((void **)&doutput_array, sizeof(double2)*3);
    cudaMalloc((void **)&devolution_type, sizeof(int)*3);

    cudaMemcpy(dinput_array, input_array, sizeof(double)*3,
               cudaMemcpyHostToDevice);
    cudaMemcpy(devolution_type, evolution_type, sizeof(int)*3,
               cudaMemcpyHostToDevice);

    dim3 threads = {1,1,1};
    dim3 grid = {1,1,1};

    make_complex_kernel<<<1,1>>>(dinput_array, devolution_type,
                                           doutput_array);
    cudaDeviceSynchronize();

    cudaMemcpy(output_array, doutput_array, sizeof(double2)*3, 
               cudaMemcpyDeviceToHost);

    bool pass = true;
    double thresh = 0.000001;

    if (abs(output_array[0].x - input_array[0]) > thresh || 
        (output_array[0].y) > thresh){
        std::cout << "failed 1\n";
        pass = false;
    }
    if (abs(output_array[1].x - exp(-input_array[1])) > thresh || 
        abs(output_array[1].y) > thresh){
        std::cout << "failed 2\n";
        pass = false;
    }
    if (abs(output_array[2].x - cos(-input_array[2])) > thresh || 
        abs(output_array[2].y - sin(-input_array[2])) > thresh){
        std::cout << "failed 3\n";
        pass = false;
    }

    if(pass){
        std::cout << "make_complex test passed!\n";
    }
    else{
        std::cout << "make_complex test failed!\n";
        exit(1);
    }
    
}

void cMultPhi_test(){
    // first, we are creating a double2 array to work with
    int n = 32;
    double2 *in1, *out;
    double *in2;
    double2 *din1, *dout;
    double *din2;

    in1 = (double2 *)malloc(sizeof(double2)*n);
    in2 = (double *)malloc(sizeof(double)*n);
    out = (double2 *)malloc(sizeof(double2)*n);

    cudaMalloc((void **)&din1, sizeof(double2)*n);
    cudaMalloc((void **)&din2, sizeof(double)*n);
    cudaMalloc((void **)&dout, sizeof(double2)*n);

    for (int i = 0; i < n; ++i){
        in1[i].x = i;
        in1[i].y = n-i;
        in2[i] = n-i;
    }

    cudaMemcpy(din1, in1, sizeof(double2)*n, cudaMemcpyHostToDevice);
    cudaMemcpy(din2, in2, sizeof(double)*n, cudaMemcpyHostToDevice);

    cMultPhi<<<1,n>>>(din1, din2, dout);
    cudaDeviceSynchronize();

    cudaMemcpy(out, dout, sizeof(double2)*n, cudaMemcpyDeviceToHost);

    double thresh = 0.000001;
    bool result = true;
    for (int i = 0; i < n; ++i){
        if (abs(out[i].x-cos(in2[i])*in1[i].x-in1[i].y*sin(in2[i])) < thresh ||
            abs(out[i].y-in1[i].x*sin(in2[i])+in1[i].y*cos(in2[i])) < thresh){
            result = false;
        }
    }

    if (result){
        std::cout << "cMultPhi test passed!\n";
    }
    else{
        std::cout << "cMultPhi test failed!\n";
        exit(1);
    }

}

// Test for available amount of GPU memory
void check_memory_test(){
    Grid par;
    par.store("xDim",10);
    par.store("yDim",10);
    par.store("zDim",10);

    par.store("energy_calc",true);

    check_memory(par);

    std::cout << "CUDA memory check passed!\n";
}
