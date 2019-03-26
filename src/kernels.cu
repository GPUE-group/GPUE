#include "../include/constants.h"
#include "../include/dynamic.h"
#include "../include/ds.h"
#include <stdio.h>

__device__ double2 subtract(double2 a, double2 b){
    return {a.x-b.x, a.y-b.y};
}
__device__ double2 add(double2 a, double2 b){
    return {a.x+b.x, a.y+b.y};
}
__device__ double2 pow(double2 a, int b){
    double r = sqrt(a.x*a.x + a.y*a.y);
    double theta = atan(a.y / a.x);
    return{pow(r,b)*cos(b*theta),pow(r,b)*sin(b*theta)};
}

__device__ double2 mult(double2 a, double b){
    return {a.x*b, a.y*b};
}
__device__ double2 mult(double2 a, double2 b){
    return {a.x*b.x - a.y*b.y, a.x*b.y + a.y*b.x};
}

//Evaluted in MATLAB: N*4*HBAR*HBAR*PI*(4.67e-9/mass)*sqrt(mass*(omegaZ)/(2*PI*HBAR))
//__constant__ double gDenConst = 6.6741e-40;

__device__ unsigned int getGid3d3d(){
    int blockId = blockIdx.x + blockIdx.y * gridDim.x
                  + gridDim.x * gridDim.y * blockIdx.z;
    int threadId = blockId * (blockDim.x * blockDim.y * blockDim.z)
                   + (threadIdx.y * blockDim.x)
                   + (threadIdx.z * (blockDim.x * blockDim.y)) + threadIdx.x;
    return threadId;
}

// global kernel to perform derivative
// For xDim derivative, stride = 1
// For yDim derivative, stride = xDim
// For zDim derivative, stride = xDim*yDim
__global__ void derive(double *data, double *out, int stride, int gsize,
                       double dx){
    int gid = getGid3d3d();
    if (gid + stride < gsize){
        out[gid] = (data[gid+stride] - data[gid])/dx;
    }
    else{
        out[gid] = data[gid]/dx;
    }
}

// global kernel to perform derivative
// For xDim derivative, stride = 1
// For yDim derivative, stride = xDim
// For zDim derivative, stride = xDim*yDim
__global__ void derive(double2 *data, double2 *out, int stride, int gsize,
                       double dx){
    int gid = getGid3d3d();
    if (gid + stride < gsize){
        out[gid].x = (data[gid+stride].x - data[gid].x)/dx;
        out[gid].y = (data[gid+stride].y - data[gid].y)/dx;
    }
    else{
        out[gid].x = data[gid].x/dx;
        out[gid].y = data[gid].y/dx;
    }
}

__global__ void is_eq(bool *a, bool *b, bool *ans){
    int gid = getGid3d3d();
    ans[0] = true;
    if (a[gid] != b[gid]){
        ans[0] = false;
    }
}


// Function to convert a double* to double2*
__global__ void make_cufftDoubleComplex(double *in, double2 *out){
    int gid = getGid3d3d();
    out[gid].x = in[gid];
    out[gid].y = 0;
}

// Function to copy double2* values
__global__ void copy(double2 *in, double2 *out){
    int gid = getGid3d3d();
    out[gid] = in[gid];
}


// function to perform a transposition (2d) or permutation (3d)
// Note: The 3 ints represent the final placement of that data direction
//       after transposition
inline __device__ unsigned int permuteGid(int d1, int d2, int d3){

    // I cannot seem to think of any way to write this in a general case...

    unsigned int x, y, z;

    // If the three axes are in the original directions.
    if (d1 == 0 && d2 == 1 && d3 == 2){
        return getGid3d3d();
    } 

    else if (d1 == 1 && d2 == 2 && d3 == 0){
        x = blockIdx.x * blockDim.x + threadIdx.x;
        z = blockDim.z * (x + blockIdx.z) + threadIdx.z;
        y = blockDim.y * (z + blockIdx.y) + threadIdx.y;
        return y;
    }

    else if (d1 == 2 && d2 == 0 && d3 == 1){
        y = blockIdx.y * blockDim.y + threadIdx.y;
        x = blockDim.x * (y + blockIdx.x) + threadIdx.x;
        z = blockDim.z * (x + blockIdx.z) + threadIdx.z;
        return z;
    }

    else if (d1 == 0 && d2 == 2 && d3 == 1){
        y = blockIdx.y * blockDim.y + threadIdx.y;
        z = blockDim.z * (y + blockIdx.z) + threadIdx.z;
        x = blockDim.x * (z + blockIdx.x) + threadIdx.x;
        return x;
    }

    else if (d1 == 1 && d2 == 0 && d3 == 2){
        //z = blockIdx.z * blockDim.z + threadIdx.z;
        //y = blockDim.y * (z + blockIdx.x) + threadIdx.y;
        //x = blockDim.x * (y + blockIdx.y) + threadIdx.x;
        //x = blockDim.x * (z + blockIdx.x) + threadIdx.x;
        //y = blockDim.y * (x + blockIdx.y) + threadIdx.y;
        //return y;
        //return x;
        z = blockIdx.z*blockDim.z + threadIdx.z;
        x = blockIdx.x*blockDim.x + threadIdx.x;
        y = blockIdx.y*blockDim.y + threadIdx.y;
        return x + blockDim.x*y + blockDim.y*blockDim.x*z;
    }

    else if (d1 == 2 && d2 == 1 && d3 == 0){
        x = blockIdx.x * blockDim.x + threadIdx.x;
        y = blockDim.y * (x + blockIdx.y) + threadIdx.y;
        z = blockDim.z * (y + blockIdx.z) + threadIdx.z;
        return z;
    }
    else{
        return 0;
    }

}

__device__ unsigned int getBid3d3d(){
    return blockIdx.x + gridDim.x*(blockIdx.y + gridDim.y * blockIdx.z);
}

__device__ unsigned int getTid3d3d(){
    return blockDim.x * ( blockDim.y * ( blockDim.z + ( threadIdx.z * blockDim.y ) )  + threadIdx.y )  + threadIdx.x;
}

__device__ double2 make_complex(double in, int evolution_type){
    double2 result;

    switch(evolution_type){
        // No change
        case 0:
            result.x = in;
            result.y = 0;
            break;
        // Im. Time evolution
        case 1:
            result.x = exp(-in);
            result.y = 0;
            break;
        // Real Time evolution
        case 2:
            result.x = cos(-in);
            result.y = sin(-in);
            break;
    }

    return result;
}

__device__ double2 conjugate(double2 in){
    double2 result = in;
    result.y = -result.y;
    return result;
}

__device__ double2 realCompMult(double scalar, double2 comp){
    double2 result;
    result.x = scalar * comp.x;
    result.y = scalar * comp.y;
    return result;
}

__device__ double complexMagnitude(double2 in){
    return sqrt(in.x*in.x + in.y*in.y);
}

__global__ void energy_sum(double2 *in1, double2 *in2, double *out){
    int gid = getGid3d3d();
    out[gid] = in1[gid].x + in2[gid].x;
}

__global__ void energy_lsum(double *in1, double2 *in2, double *out){
    int gid = getGid3d3d();
    out[gid] = in1[gid] + in2[gid].x;
}

__global__ void sum(double2 *in1, double2 *in2, double2 *out){
    int gid = getGid3d3d();
    out[gid].x = in1[gid].x + in2[gid].x;
    out[gid].y = in1[gid].y + in2[gid].y;
}

__global__ void complexAbsSum(double2 *in1, double2 *in2, double *out){
    int gid = getGid3d3d();
    double2 temp;
    temp.x = in1[gid].x + in2[gid].x;
    temp.y = in1[gid].y + in2[gid].y;
    out[gid] = sqrt(temp.x*temp.x + temp.y*temp.y);
}

__global__ void complexAbsSum(double2 *in1, double2 *in2, double2 *in3,
                              double *out){
    int gid = getGid3d3d();
    double2 temp;
    temp.x = in1[gid].x + in2[gid].x + in3[gid].x;
    temp.y = in1[gid].y + in2[gid].y + in3[gid].y;
    out[gid] = sqrt(temp.x*temp.x + temp.y*temp.y);
}

__global__ void complexMagnitude(double2 *in, double *out){
    int gid = getGid3d3d();
    out[gid] = sqrt(in[gid].x*in[gid].x + in[gid].y*in[gid].y);
}

__host__ __device__ double complexMagnitudeSquared(double2 in){
    return in.x*in.x + in.y*in.y;
}

__global__ void complexMagnitudeSquared(double2 *in, double *out){
    int gid = getGid3d3d();
    out[gid] = in[gid].x*in[gid].x + in[gid].y*in[gid].y;
}

__global__ void complexMagnitudeSquared(double2 *in, double2 *out){
    int gid = getGid3d3d();
    out[gid].x = in[gid].x*in[gid].x + in[gid].y*in[gid].y;
    out[gid].y = 0;
}

__host__ __device__ double2 complexMultiply(double2 in1, double2 in2){
    double2 result;
    result.x = (in1.x*in2.x - in1.y*in2.y);
    result.y = (in1.x*in2.y + in1.y*in2.x);
    return result;
}

__global__ void complexMultiply(double2 *in1, double2 *in2, double2 *out){
    int gid = getGid3d3d();
    out[gid] = complexMultiply(in1[gid], in2[gid]);
}


/*
* Used to perform conj(in1)*in2; == < in1 | in2 >
*/
inline __device__ double2 braKetMult(double2 in1, double2 in2){
    return complexMultiply(conjugate(in1),in2);
}

/**
 * Performs complex multiplication of in1 and in2, giving result as out. 
 */
__global__ void cMult(double2* in1, double2* in2, double2* out){
    unsigned int gid = getGid3d3d();
    double2 result;
    double2 tin1 = in1[gid];
    double2 tin2 = in2[gid];
    result.x = (tin1.x*tin2.x - tin1.y*tin2.y);
    result.y = (tin1.x*tin2.y + tin1.y*tin2.x);
    out[gid] = result;
}

__global__ void cMultPhi(double2* in1, double* in2, double2* out){
    double2 result;
    unsigned int gid = getGid3d3d();
    result.x = cos(in2[gid])*in1[gid].x - in1[gid].y*sin(in2[gid]);
    result.y = in1[gid].x*sin(in2[gid]) + in1[gid].y*cos(in2[gid]);
    out[gid] = result;
}

/**
 * Performs multiplication of double* with double2*
 */
__global__ void vecMult(double2 *in, double *factor, double2 *out){
    double2 result;
    unsigned int gid = getGid3d3d();
    result.x = in[gid].x * factor[gid];
    result.y = in[gid].y * factor[gid];
    out[gid] = result;
}

__global__ void vecMult(double *in, double *factor, double *out){
    double result;
    unsigned int gid = getGid3d3d();
    result = in[gid] * factor[gid];
    out[gid] = result;
}


__global__ void vecSum(double2 *in, double *factor, double2 *out){
    double2 result;
    unsigned int gid = getGid3d3d();
    result.x = in[gid].x + factor[gid];
    result.y = in[gid].y + factor[gid];
    out[gid] = result;
}

__global__ void vecSum(double *in, double *factor, double *out){
    double result;
    unsigned int gid = getGid3d3d();
    result = in[gid] + factor[gid];
    out[gid] = result;
}


__global__ void l2_norm(double *in1, double *in2, double *in3, double *out){

    int gid = getGid3d3d();
    out[gid] = sqrt(in1[gid]*in1[gid] + in2[gid]*in2[gid] + in3[gid]*in3[gid]);
}

__global__ void l2_norm(double2 *in1, double2 *in2, double2 *in3, double *out){

    int gid = getGid3d3d();
    out[gid] = sqrt(in1[gid].x*in1[gid].x + in1[gid].y*in1[gid].y
                    + in2[gid].x*in2[gid].x + in2[gid].y*in2[gid].y
                    + in3[gid].x*in3[gid].x + in3[gid].y*in3[gid].y);
}

__global__ void l2_norm(double *in1, double *in2, double *out){

    int gid = getGid3d3d();
    out[gid] = sqrt(in1[gid]*in1[gid] + in2[gid]*in2[gid]);
}

__global__ void l2_norm(double2 *in1, double2 *in2, double *out){

    int gid = getGid3d3d();
    out[gid] = sqrt(in1[gid].x*in1[gid].x + in1[gid].y*in1[gid].y
                    + in2[gid].x*in2[gid].x + in2[gid].y*in2[gid].y);
}

/**
 * Performs the non-linear evolution term of Gross--Pitaevskii equation.
 */
__global__ void cMultDensity(double2* in1, double2* in2, double2* out, double dt, int gstate, double gDenConst){
    double2 result;
    double gDensity;

    int gid = getGid3d3d();
    double2 tin1 = in1[gid];
    double2 tin2 = in2[gid];
    gDensity = gDenConst*complexMagnitudeSquared(in2[gid])*(dt/HBAR);

    if(gstate == 0){
        double tmp = in1[gid].x*exp(-gDensity);
        result.x = (tmp)*tin2.x - (tin1.y)*tin2.y;
        result.y = (tmp)*tin2.y + (tin1.y)*tin2.x;
    }
    else{
        double2 tmp;
        tmp.x = tin1.x*cos(-gDensity) - tin1.y*sin(-gDensity);
        tmp.y = tin1.y*cos(-gDensity) + tin1.x*sin(-gDensity);
        
        result.x = (tmp.x)*tin2.x - (tmp.y)*tin2.y;
        result.y = (tmp.x)*tin2.y + (tmp.y)*tin2.x;
    }
    out[gid] = result;
}

//cMultDensity for ast V
__global__ void cMultDensity_ast(EqnNode_gpu *eqn, double2* in, double2* out, 
                                 double dx, double dy, double dz, double time,
                                 int e_num, double dt, int gstate, 
                                 double gDenConst){
    double2 result;
    double gDensity;

    int gid = getGid3d3d();
    int xid = blockIdx.x*blockDim.x + threadIdx.x;
    int yid = blockIdx.y*blockDim.y + threadIdx.y;
    int zid = blockIdx.z*blockDim.z + threadIdx.z;

    double2 tin = in[gid];
    gDensity = gDenConst*complexMagnitudeSquared(in[gid])*(dt/HBAR);
    double2 val = make_complex(evaluate_eqn_gpu(eqn, xid*dx, yid*dy, zid*dz,
                                                time, e_num), gstate+1);

    if(gstate == 0){
        double tmp = val.x*exp(-gDensity);
        result.x = (tmp)*tin.x - (val.y)*tin.y;
        result.y = (tmp)*tin.y + (val.y)*tin.x;
    }
    else{
        double2 tmp;
        tmp.x = val.x*cos(-gDensity) - val.y*sin(-gDensity);
        tmp.y = val.y*cos(-gDensity) + val.x*sin(-gDensity);

        result.x = (tmp.x)*tin.x - (tmp.y)*tin.y;
        result.y = (tmp.x)*tin.y + (tmp.y)*tin.x;
    }
    out[gid] = result;
}

/**
 * Divides both components of vector type "in", by the value "factor".
 * Results given with "out".
 */
__global__ void scalarDiv(double2* in, double factor, double2* out){
    double2 result;
    unsigned int gid = getGid3d3d();
    result.x = (in[gid].x / factor);
    result.y = (in[gid].y / factor);
    out[gid] = result;
}

/**
 * Divides both components of vector type "in", by the value "factor".
 * Results given with "out".
 */
__global__ void scalarDiv(double* in, double factor, double* out){
    double result;
    unsigned int gid = getGid3d3d();
    result = (in[gid] / factor);
    out[gid] = result;
}


/**
 * Multiplies both components of vector type "in", by the value "factor".
 * Results given with "out". 
 */
__global__ void scalarMult(double2* in, double factor, double2* out){
    double2 result;
    //extern __shared__ double2 tmp_in[];
    unsigned int gid = getGid3d3d();
    result.x = (in[gid].x * factor);
    result.y = (in[gid].y * factor);
    out[gid] = result;
}

__global__ void scalarMult(double* in, double factor, double* out){
    double result;
    unsigned int gid = getGid3d3d();
    result = (in[gid] * factor);
    out[gid] = result;
}


/**
 * As above, but normalises for wfc
 */
__global__ void scalarDiv_wfcNorm(double2* in, double dr, double* pSum, double2* out){
    unsigned int gid = getGid3d3d();
    double2 result;
    double norm = sqrt((pSum[0])*dr);
    result.x = (in[gid].x/norm);
    result.y = (in[gid].y/norm);
    out[gid] = result;
}

/**
 * Raises in to the power of param
 */
__global__ void scalarPow(double2* in, double param, double2* out){
    unsigned int gid = getGid3d3d();
    double2 result;
    result.x = pow(in[gid].x, param);
    result.y = pow(in[gid].y, param);
    out[gid] = result;
}

/**
 * Finds conjugate for double2*
 */
__global__ void vecConjugate(double2 *in, double2 *out){
    double2 result;
    unsigned int gid = getGid3d3d(); 
    result.x = in[gid].x;
    result.y = -in[gid].y;
    out[gid] = result;
}

__global__ void angularOp(double omega, double dt, double2* wfc, double* xpyypx, double2* out){
    unsigned int gid = getGid3d3d();
    double2 result;
    double op;
    op = exp( -omega*xpyypx[gid]*dt);
    result.x=wfc[gid].x*op;
    result.y=wfc[gid].y*op;
    out[gid]=result;
}

/**
 * Kernel for a quick test of the threads and such for GPU computing
 */
__global__ void thread_test(double *in, double *out){

    unsigned int Gid = getGid3d3d();

    // Now we set each element in the 
    out[Gid] = Gid;
    //in[Gid] = Gid;
}

/**
 * Routine for parallel summation. Can be looped over from host.
 */
__global__ void multipass(double2* input, double2* output, int pass){
    unsigned int tid = threadIdx.x + threadIdx.y*blockDim.x 
                       + threadIdx.z * blockDim.x * blockDim.y;
    unsigned int bid = blockIdx.x + blockIdx.y * gridDim.x
                       + gridDim.x * gridDim.y * blockIdx.z;

    //unsigned int tid = getTid3d3d();
    //unsigned int bid = getBid3d3d();
    // printf("bid0=%d\n",bid);

    unsigned int gid = getGid3d3d();
    extern __shared__ double2 sdata[];
    sdata[tid] = input[gid];
    __syncthreads();
    for(int i = blockDim.x>>1; i > 0; i>>=1){
        if(tid < i){
            sdata[tid].x += sdata[tid + i].x;
            sdata[tid].y += sdata[tid + i].y;
        }
        __syncthreads();
    }
    if(tid==0){
        output[bid] = sdata[0];
    }
}

/**
 * Routine for parallel summation. Can be looped over from host.
 */
__global__ void multipass(double* input, double* output){
    unsigned int tid = threadIdx.x + threadIdx.y*blockDim.x
                       + threadIdx.z * blockDim.x * blockDim.y;
    unsigned int bid = blockIdx.x + blockIdx.y * gridDim.x
                       + gridDim.x * gridDim.y * blockIdx.z;

    //unsigned int tid = getTid3d3d();
    //unsigned int bid = getBid3d3d();
    // printf("bid0=%d\n",bid);

    unsigned int gid = getGid3d3d();
    extern __shared__ double sdatad[];
    sdatad[tid] = input[gid];
    __syncthreads();

    for(int i = blockDim.x>>1; i > 0; i>>=1){
        if(tid < i){
            sdatad[tid] += sdatad[tid + i];
        }
        __syncthreads();
    }
    if(tid==0){
        output[bid] = sdatad[0];
    }
}

/*
* Calculates all of the energy of the current state. sqrt_omegaz_mass = sqrt(omegaZ/mass), part of the nonlin interaction term
*/
__global__ void energyCalc(double2 *wfc, double2 *op, double dt, double2 *energy, int gnd_state, int op_space, double sqrt_omegaz_mass, double gDenConst){
    unsigned int gid = getGid3d3d();
    //double hbar_dt = HBAR/dt;
    double2 result;
    if(op_space == 0){
        double g_local = gDenConst*sqrt_omegaz_mass*complexMagnitudeSquared(wfc[gid]);
        op[gid].x += g_local;
    }

    if (op_space < 2){
        result = braKetMult(wfc[gid], energy[gid]);
        energy[gid].x += result.x;
        energy[gid].y += result.y;
    }
    else{
        result = complexMultiply(op[gid],wfc[gid]);
    }
    result = braKetMult(wfc[gid], complexMultiply(op[gid],wfc[gid]));
    energy[gid].x += result.x;
    energy[gid].y += result.y;

}

// Function to multiply a double* with an astval
__global__ void ast_mult(double *array, double *array_out, EqnNode_gpu *eqn,
                         double dx, double dy, double dz, double time,
                         int element_num){
    int gid = getGid3d3d();
    int xid = blockIdx.x*blockDim.x + threadIdx.x;
    int yid = blockIdx.y*blockDim.y + threadIdx.y;
    int zid = blockIdx.z*blockDim.z + threadIdx.z;

    double val = evaluate_eqn_gpu(eqn, xid*dx, yid*dy, zid*dz, 
                                  time, element_num);

    array_out[gid] = array[gid] * val;
}

// Function to multiply a double* with an astval
__global__ void ast_cmult(double2 *array, double2 *array_out, EqnNode_gpu *eqn,
                          double dx, double dy, double dz, double time,
                          int element_num){
    int gid = getGid3d3d();
    int xid = blockIdx.x*blockDim.x + threadIdx.x;
    int yid = blockIdx.y*blockDim.y + threadIdx.y;
    int zid = blockIdx.z*blockDim.z + threadIdx.z;

    double val = evaluate_eqn_gpu(eqn, xid*dx, yid*dy, zid*dz, 
                                  time, element_num);

    array_out[gid].x = array[gid].x * val;
    array_out[gid].y = array[gid].y * val;
}

// Function to multiply an AST V in real or imaginary time evolution
__global__ void ast_op_mult(double2 *array, double2 *array_out, 
                            EqnNode_gpu *eqn,
                            double dx, double dy, double dz, double time,
                            int element_num, int evolution_type, double dt){
    int gid = getGid3d3d();
    int xid = blockIdx.x*blockDim.x + threadIdx.x;
    int yid = blockIdx.y*blockDim.y + threadIdx.y;
    int zid = blockIdx.z*blockDim.z + threadIdx.z;

    double val = evaluate_eqn_gpu(eqn, xid*dx, yid*dy, zid*dz, 
                                  time, element_num);
    double2 complex_val = make_complex(val*dt, evolution_type);

    array_out[gid].x = array[gid].x * complex_val.x 
                       - array[gid].y * complex_val.y;
    array_out[gid].y = array[gid].x * complex_val.y 
                       + array[gid].y * complex_val.x;
}


// Function to find the ast in real-time dynamics
__device__ double2 real_ast(double val, double dt){

    return {cos(-val*dt), sin(-val*dt)};
}

// Function to find the ast in real-time dynamics
__device__ double2 im_ast(double val, double dt){

    return {exp(-val*dt), 0};
}

__global__ void zeros(bool *in, bool *out){
    int gid = getGid3d3d();
    out[gid] = 0;
}

__global__ void zeros(double *in, double *out){
    int gid = getGid3d3d();
    out[gid] = 0;
}

__global__ void zeros(double2 *in, double2 *out){
    int gid = getGid3d3d();
    out[gid].x = 0;
    out[gid].y = 0;
}

__global__ void set_eq(double *in1, double *in2){
    int gid = getGid3d3d();
    in2[gid] = in1[gid];
}

//##############################################################################
//##############################################################################

/**
 * Routine for parallel summation. Can be looped over from host.
 */
template<typename T> __global__ void pSumT(T* in1, T* output, int pass){
        unsigned int tid = threadIdx.x;
        unsigned int bid = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x;// printf("bid0=%d\n",bid);
        unsigned int gid = getGid3d3d();
        extern __shared__ T sdata[];
        for(int i = blockDim.x>>1; i > 0; i>>=1){
                if(tid < blockDim.x>>1){
                        sdata[tid] += sdata[tid + i];
                }
                __syncthreads();
        }
        if(tid==0){
                output[bid] = sdata[0];
        }
}

/**
 * Routine for parallel summation. Can be looped over from host. BETA
 */
__global__ void pSum(double* in1, double* output, int pass){
        unsigned int tid = threadIdx.x;
        unsigned int bid = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x;// printf("bid0=%d\n",bid);
        unsigned int gid = getGid3d3d();
        extern __shared__ double sdata2[];
        for(int i = blockDim.x>>1; i > 0; i>>=1){
                if(tid < blockDim.x>>1){
                        sdata2[tid] += sdata2[tid + i];
                }
                __syncthreads();
        }
        if(tid==0){
                output[bid] = sdata2[0];
        }
}

//##############################################################################
//##############################################################################
