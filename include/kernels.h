///@endcond
//##############################################################################
/**
 *  @file    kernels.h
 *  @author  Lee J. O'Riordan (mlxd)
 *  @date    12/11/2015
 *  @version 0.1
 *
 *  @brief GPU kernel definitions
 *
 *  @section DESCRIPTION
 *  Kernel definitions for all CUDA-enabled routines for solving GPE.
 */
//##############################################################################

#ifndef KERNELS_H
#define KERNELS_H
#include<stdio.h>

/**
* @brief	subtraction operation for 2 double2 values
* @ingroup	gpu
*/
__device__ double2 subtract(double2 a, double2 b);

/**
* @brief	addition operation for 2 double2 values
* @ingroup	gpu
*/
__device__ double2 add(double2 a, double2 b);

/**
* @brief	power operation for a double2
* @param        base number
* @param        power
* @ingroup	gpu
*/
__device__ double2 pow(double2 a, int b);

/**
* @brief	multiplication operation for double2 and double values
* @ingroup	gpu
*/
__device__ double2 mult(double2 a, double b);

/**
* @brief	multiplication operation for 2 double2 values
* @ingroup	gpu
*/
__device__ double2 mult(double2 a, double2 b);

/**
* @brief	transforms an array of doubles into double2's
* @ingroup	gpu
*/
__global__ void make_cufftDoubleComplex(double *in, double2 *out);

/**
* @brief	Indexing of threads on grid
* @ingroup	gpu
*/
__device__ unsigned int getGid3d3d();

/**
* @brief	Indexing of blocks on device
* @ingroup	gpu
*/
__device__ unsigned int getBid3d3d();

/**
* @brief	Indexing of threads in a block on device
* @ingroup	gpu
*/
__device__ unsigned int getTid3d3d();

/**
* @brief	checks to arrays to see if they are equal
* @ingroup	gpu
*/
__global__ void is_eq(bool *a, bool *b, bool *ans);

//##############################################################################
/**
* Helper functions for complex numbers
*/
//##############################################################################

/**
* @brief	Calculates magnitude of complex number. $|a + ib|$
* @ingroup	gpu
* @param	in Complex number
* @return	Magnitude of complex number
*/
__device__ double complexMagnitude(double2 in);

/**
* @brief	Complex multiplication of two input arrays
* @ingroup	gpu
* @param	Input 1
* @param	Input 2
* @param	Output
*/
__global__ void complexMultiply(double2 *in1, double2 *in2, double2 *out);

/**
* @brief	Complex multiplication of two input values
* @ingroup	gpu
* @param	Input 1
* @param	Input 2
* @return	Output
*/
__host__ __device__ double2 complexMultiply(double2 in1, double2 in2);

/**
* @brief        Transforms field value into operator
* @ingroup      gpu
* @param        Input value
* @param        Evolution type (0 for imaginary, 1 for real)
* @return       complex output
*/
__device__ double2 make_complex(double in, int evolution_type);

/**
* @brief        Sums the absolute value of two complex arrays
* @ingroup      gpu
* @param        Array 1
* @param        Array 2
* @param        Output
*/
__global__ void complexAbsSum(double2 *in1, double2 *in2, double *out);

/**
* @brief        Sums the absolute value of two complex arrays
* @ingroup      gpu
* @param        Array 1
* @param        Array 2
* @param        Array 3
* @param        Output
*/
__global__ void complexAbsSum(double2 *in1, double2 *in2, double2 *in3,
                              double *out);

/**
* @brief        Complex magnitude of a double2 array
* @ingroup      gpu
*/
__global__ void complexMagnitude(double2 *in, double *out);

/**
* @brief	Return the squared magnitude of a complex number. $|(a+\textrm{i}b)*(a-\textrm{i}b)|$
* @ingroup	gpu
* @param	in Complex number
* @return	Absolute-squared complex number
*/
__device__ double complexMagnitudeSquared(double2 in);

/**
* @brief        Complex magnitude of a double2 array
* @ingroup      gpu
*/
__global__ void complexMagnitudeSquared(double2 *in, double *out);

/**
* @brief        Complex magnitude of a double2 array
* @ingroup      gpu
*/
__global__ void complexMagnitudeSquared(double2 *in, double2 *out);

/**
* @brief	Returns conjugate of the a complex number
* @ingroup	gpu
* @param	in Number to be conjugated
* @return	Conjugated complex number
*/
__device__ double2 conjugate(double2 in);

/**
* @brief	Multiply real scalar by a complex number
* @ingroup	gpu
* @param	scalar Scalar multiplier
* @param	comp Complex multiplicand
* @return	Result of scalar * comp
*/
__device__ double2 realCompMult(double scalar, double2 comp);

//##############################################################################
/**
 * Multiplication for linear, non-linear and phase-imprinting of the condensate.
 */
//##############################################################################

/**
* @brief	Kernel for complex multiplication
* @ingroup	gpu
* @param	in1 Wavefunction input
* @param	in2 Evolution operator input
* @param	out Pass by reference output for multiplcation result
*/
__global__ void cMult(double2* in1, double2* in2, double2* out);

/**
* @brief	Kernel for multiplcation with real array and complex array
* @ingroup	gpu
* @param	in1 Wavefunction input
* @param	in2 Evolution operator input
* @param	out Pass by reference output for multiplcation result
*/
__global__ void cMultPhi(double2* in1, double* in2, double2* out);

/**
* @brief	Kernel for complex multiplication with nonlinear density term
* @ingroup	gpu
* @param	in1 Wavefunction input
* @param	in2 Evolution operator input
* @param	out Pass by reference output for multiplication result
* @param	dt Timestep for evolution
* @param	gState If performing real (1) or imaginary (0) time evolution
* @param	gDenConst a constant for evolution
*/
__global__ void cMultDensity(double2* in1, double2* in2, double2* out, double dt, int gstate, double gDenConst);

/**
* @brief        Kernel for complex multiplication with nonlinear density term
* @ingroup      gpu
* @param        GPU AST
* @param        in Wavefunction input
* @param        out Wavefunction output
* @param        dx
* @param        dy
* @param        dz
* @param        time
* @param        element number in AST
* @param        dt Timestep for evolution
* @param        gState If performing real (1) or imaginary (0) time evolution
* @param        gDenConst a constant for evolution
*/

__global__ void cMultDensity_ast(EqnNode_gpu *eqn, double2* in, double2* out,
                                 double dx, double dy, double dz, double time,
                                 int e_num, double dt, int gstate,
                                 double gDenConst);


/**
* @brief	Hold vortex at specified position. Not implemented. cMultPhi should implement required functionality.
* @ingroup	gpu
* @param	in1 Wavefunction input
* @param	in2 Evolution operator input
* @param	out Pass by reference output for multiplcation result
*/
__global__ void pinVortex(double2* in1, double2* in2, double2* out);


/**
* @brief        Complex field scaling and renormalisation. Used mainly post-FFT.
* @ingroup      gpu
* @param        in Complex field to be scaled (divided, not multiplied)
* @param        factor Scaling vector to be used
* @param        out Pass by reference output for result
*/
__global__ void vecMult(double2 *in, double *factor, double2 *out);

/**
* @brief        Complex field Summation
* @ingroup      gpu
* @param        in Complex field to be scaled (divided, not multiplied)
* @param        factor Scaling vector to be used
* @param        out Pass by reference output for result
*/
__global__ void vecSum(double2 *in, double *factor, double2 *out);

/**
* @brief        field scaling
* @ingroup      gpu
* @param        in field to be scaled (divided, not multiplied)
* @param        factor Scaling vector to be used
* @param        out Pass by reference output for result
*/
__global__ void vecMult(double *in, double *factor, double *out);

/**
* @brief        field Summation
* @ingroup      gpu
* @param        in field to be scaled (divided, not multiplied)
* @param        factor Scaling vector to be used
* @param        out Pass by reference output for result
*/
__global__ void vecSum(double *in, double *factor, double *out);

/**
* @brief        performs the l2 normalization of the provided terms
* @ingroup      gpu
*/
__global__ void l2_norm(double *in1, double *in2, double *in3, double *out);

/**
* @brief        performs the l2 normalization of the provided terms
* @ingroup      gpu
*/
__global__ void l2_norm(double2 *in1, double2 *in2, double2 *in3, double *out);

/**
* @brief        performs the l2 normalization of the provided terms
* @ingroup      gpu
*/
__global__ void l2_norm(double *in1, double *in2, double *out);

/**
* @brief        performs the l2 normalization of the provided terms
* @ingroup      gpu
*/
__global__ void l2_norm(double2 *in1, double2 *in2, double *out);

/**
* @brief	Complex field scaling and renormalisation. Used mainly post-FFT.
* @ingroup	gpu
* @param	in Complex field to be scaled (divided, not multiplied)
* @param	factor Scaling factor to be used
* @param	out Pass by reference output for result
*/
__global__ void scalarDiv(double2* in, double factor, double2* out);

/**
* @brief        Real field scaling and renormalisation. Used mainly post-FFT.
* @ingroup      gpu
* @param        in Real field to be scaled (divided, not multiplied)
* @param        factor Scaling factor to be used
* @param        out Pass by reference output for result
*/
__global__ void scalarDiv(double* in, double factor, double* out);

/**
* @brief        Complex field scaling and renormalisation. Used mainly post-FFT.
* @ingroup      gpu
* @param        in Complex field to be scaled (multiplied, not divided)
* @param        factor Scaling factor to be used
* @param        out Pass by reference output for result
*/
__global__ void scalarMult(double2* in, double factor, double2* out);

/**
* @brief        field scaling and renormalisation. Used mainly post-FFT.
* @ingroup      gpu
* @param        in field to be scaled (multiplied, not divided)
* @param        factor Scaling factor to be used
* @param        out Pass by reference output for result
*/
__global__ void scalarMult(double* in, double factor, double* out);

/**
* @brief        Complex field raised to a power
* @ingroup      gpu
* @param        in Complex field to be scaled (multiplied, not divided)
* @param        power parameter
* @param        out Pass by reference output for result
*/
__global__ void scalarPow(double2* in, double param, double2* out);

/**
* @brief        Conjugate of double2*.
* @ingroup      gpu
* @param        in Complex field to be conjugated
* @param        out Pass by reference output for result
*/
__global__ void vecConjugate(double2 *in, double2 *out);

/**
* @brief	Complex field scaling and renormalisation. Not implemented. Use scalarDiv
* @ingroup	gpu
*/
__global__ void scalarDiv1D(double2*, double2*);

/**
* @brief	Complex field scaling and renormalisation. Not implemented. Use scalarDiv
* @ingroup	gpu
*/
__global__ void scalarDiv2D(double2*, double2*);

/**
* @brief	Used as part of multipass to renormalise the wavefucntion
* @ingroup	gpu
* @param	in Complex field to be renormalised
* @param	dr Smallest area element of grid (dx*dy)
* @param	pSum GPU array used to store intermediate results during parallel summation
*/
__global__ void scalarDiv_wfcNorm(double2* in, double dr, double* pSum, double2* out);

//##############################################################################

/**
* @brief	Not implemented
* @ingroup	gpu
* @param	in Input field
* @param	out Output values
*/
__global__ void reduce(double2* in, double* out);

/**
* @brief        Performs wavefunction renormalisation using parallel summation and applying scalarDiv_wfcNorm
* @ingroup      gpu
* @param        input Wavefunction to be renormalised
* @param        output Pass by reference return of renormalised wavefunction
* @param        pass Number of passes performed by routine
*/
__global__ void thread_test(double* input, double* output);

/**
* @brief	Performs wavefunction renormalisation using parallel summation and applying scalarDiv_wfcNorm
* @ingroup	gpu
* @param	input Wavefunction to be renormalised
* @param	output Pass by reference return of renormalised wavefunction
* @param	pass Number of passes performed by routine
*/
__global__ void multipass(double2* input, double2* output, int pass);

/**
* @brief        Performs parallel summation of double arrays
* @ingroup      gpu
* @param        input double array
* @param        input double array after summation
*/
__global__ void multipass(double* input, double* output);

//##############################################################################

/**
* @brief	Calculates angular momentum. Not fully implemented. Handled in post-processing instead.
* @ingroup	gpu
* @param	omega Harmonic trap rotation frequency
* @param	dt Time-step for evolution
* @param	wfc Wavefunction
* @param	xpyypx L_z operator
* @param	out Output of calculation
*/
__global__ void angularOp(double omega, double dt, double2* wfc, double* xpyypx, double2* out);

/**
* @brief        Multiplication of array with AST
* @ingroup      gpu
*/
__global__ void ast_mult(double *array, double *array_out, EqnNode_gpu *eqn,
                         double dx, double dy, double dz, double time,
                         int element_num);
/**
* @brief        Complex multiplication of array with AST
* @ingroup      gpu
*/
__global__ void ast_cmult(double2 *array, double2 *array_out, EqnNode_gpu *eqn,
                          double dx, double dy, double dz, double time,
                          int element_num);
/**
* @brief        Multiplication of array with AST Operator
* @ingroup      gpu
*/
__global__ void ast_op_mult(double2 *array, double2 *array_out,
                            EqnNode_gpu *eqn,
                            double dx, double dy, double dz, double time,
                            int element_num, int evolution_type, double dt);

/**
* @brief        Function to find AST operator in real-time
* @ingroup      gpu
*/
__device__ double2 real_ast(double val, double dt);

/**
* @brief        Function to find AST operator in imaginary-time
* @ingroup      gpu
*/
__device__ double2 im_ast(double val, double dt);

/**
* @brief        Sets boolean array to 0
* @ingroup      gpu
*/
__global__ void zeros(bool *in, bool *out);

/**
* @brief        Sets in2 to be equal to in1
* @ingroup      gpu
*/
__global__ void set_eq(double *in1, double *in2);

/**
* @brief	Calculates energy of the current state during evolution. Not implemented.
* @ingroup	gpu
* @param	wfc Wavefunction
* @param	op Operator to calculate energy for.
* @param	dt Time-step for evolution
* @param	energy Energy result output
* @param	gnd_state Wavefunction
* @param	op_space Check if position space with non-linear term or not.
* @param	sqrt_omegaz_mass sqrt(omegaZ/mass), part of the nonlin interaction term.
*/
__global__ void energyCalc(double2 *wfc, double2 *op, double dt, double2 *energy, int gnd_state, int op_space, double sqrt_omegaz_mass, double gDenConst);

/**
* @brief	Performs bra-ket state multiplication. Not fully implemented.
* @ingroup	gpu
* @param	in1 Bra
* @param	in2 Ket
* @return	<Bra|Ket>
*/
inline __device__ double2 braKetMult(double2 in1, double2 in2);


/**
* @brief	Performs parallel sum. Not verified. I use multipass instead.
* @ingroup	gpu
* @param	in1 That which must be summed
* @param	output That which has been summed
* @param	pass Number of passes
*/
__global__ void pSum(double* in1, double* output, int pass);


#endif
