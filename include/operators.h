//##############################################################################
/**
 *  @file    operators.h
 *  @author  James R. Schloss (Leios)
 *  @date    1/1/2017
 *  @version 0.1
 *
 *  @brief File to hold all operators for finding fields on the GPU
 *
 *  @section DESCRIPTION
 *      This file holds all operators for finding fields on the GPU
 */
 //#############################################################################

#ifndef OPERATORS_H
#define OPERATORS_H

#include "../include/ds.h"
#include "../include/constants.h"
#include <sys/stat.h>
#include <unordered_map>
//#include <boost/math/special_functions.hpp>

// function to find Bz, the curl of the gauge field
 /**
 * @brief       Finds Bz, the curl of the gauge field
 * @ingroup     data
 * @param       Grid simulation data
 * @param       gauge field Ax
 * @param       gauge field Ay
 * @return      Bz, the curl of A
 */
double *curl2d(Grid &par, double *Ax, double *Ay);

 /**
 * @brief       Finds Bx, the curl of the gauge field
 * @ingroup     data
 * @param       Grid simulation data
 * @param       gauge field Ax
 * @param       gauge field Ay
 * @param       gauge field Az
 * @return      Bx, the curl of A
 */
double *curl3d_x(Grid &par, double *Ax, double *Ay, double *Az);

 /**
 * @brief       Finds By, the curl of the gauge field
 * @ingroup     data
 * @param       Grid simulation data
 * @param       gauge field Ax
 * @param       gauge field Ay
 * @param       gauge field Az
 * @return      By, the curl of A
 */
double *curl3d_y(Grid &par, double *Ax, double *Ay, double *Az);

 /**
 * @brief       Finds Bz, the curl of the gauge field
 * @ingroup     data
 * @param       Grid simulation data
 * @param       gauge field Ax
 * @param       gauge field Ay
 * @param       gauge field Az
 * @return      Bz, the curl of A
 */
double *curl3d_z(Grid &par, double *Ax, double *Ay, double *Az);

 /**
 * @brief       Finds Br, the curl of the gauge field
 * @ingroup     data
 * @param       Grid simulation data
 * @param       gauge field Ax
 * @param       gauge field Ay
 * @param       gauge field Az
 * @return      Br, the curl of A
 */
double *curl3d_r(Grid &par, double *Bx, double *By);

 /**
 * @brief       Finds Bphi, the curl of the gauge field
 * @ingroup     data
 * @param       Grid simulation data
 * @param       gauge field Ax
 * @param       gauge field Ay
 * @param       gauge field Az
 * @return      Bphi, the curl of A
 */
double *curl3d_phi(Grid &par, double *Bx, double *By);


 /**
 * @brief       Determines if file exists, requests new file if it does not
 * @ingroup     data
 */

std::string filecheck(std::string filename);

 /**
 * @brief       Reads A from file
 * @param       filename
 * @param       A field array
 * @param       omega multiplicative constant
 * @ingroup     data
 */
void file_A(std::string filename, double *A, double omega);

/*----------------------------------------------------------------------------//
* GPU KERNELS
*-----------------------------------------------------------------------------*/

 /**
 * @brief      Function to generate momentum grids
 */
void generate_p_space(Grid &par);

 /**
 * @brief       This function calls the appropriate K kernel
 */
void generate_K(Grid &par);

 /**
 * @brief       Simple kernel for generating K
 */
__global__ void simple_K(double *xp, double *yp, double *zp, double mass,
                         double *K);

 /**
 * @brief       Function to generate game fields
 */
void generate_gauge(Grid &par);

 /**
 * @brief       constant Kernel A
 */
__global__ void kconstant_A(double *x, double *y, double *z,
                            double xMax, double yMax, double zMax,
                            double omegaX, double omegaY, double omegaZ,
                            double omega, double fudge, double *A);

 /**
 * @brief       Kernel for simple rotational case, Ax
 */
__global__ void krotation_Ax(double *x, double *y, double *z,
                             double xMax, double yMax, double zMax,
                             double omegaX, double omegaY, double omegaZ,
                             double omega, double fudge, double *A);

 /**
 * @brief       Kernel for simple rotational case, Ay
 */
__global__ void krotation_Ay(double *x, double *y, double *z,
                             double xMax, double yMax, double zMax,
                             double omegaX, double omegaY, double omegaZ,
                             double omega, double fudge, double *A);

 /**
 * @brief       Kernel for simple triangular lattice of rings, Ax
 */
__global__ void kring_rotation_Ax(double *x, double *y, double *z,
                                  double xMax, double yMax, double zMax,
                                  double omegaX, double omegaY, double omegaZ,
                                  double omega, double fudge, double *A);

 /**
 * @brief       Kernel for simple triangular lattice of rings, Ay
 */
__global__ void kring_rotation_Ay(double *x, double *y, double *z,
                                  double xMax, double yMax, double zMax,
                                  double omegaX, double omegaY, double omegaZ,
                                  double omega, double fudge, double *A);

 /**
 * @brief       Kernel for simple triangular lattice of rings, Az
 */
__global__ void kring_rotation_Az(double *x, double *y, double *z,
                                  double xMax, double yMax, double zMax,
                                  double omegaX, double omegaY, double omegaZ,
                                  double omega, double fudge, double *A);


 /**
 * @brief       Kernel for testing Ay
 */
__global__ void ktest_Ay(double *x, double *y, double *z,
                         double xMax, double yMax, double zMax,
                         double omegaX, double omegaY, double omegaZ,
                         double omega, double fudge, double *A);

 /**
 * @brief       Kernel for testing Ax
 */
__global__ void ktest_Ax(double *x, double *y, double *z,
                         double xMax, double yMax, double zMax,
                         double omegaX, double omegaY, double omegaZ,
                         double omega, double fudge, double *A);

 /**
 * @brief       Kernel for testing Az
 */
__global__ void kring_Az(double *x, double *y, double *z,
                         double xMax, double yMax, double zMax,
                         double omegaX, double omegaY, double omegaZ,
                         double omega, double fudge, double *A);


 /**
 * @brief       Function to generate V
 */
void generate_fields(Grid &par);

 /**
 * @brief       Kernel to generate harmonic V
 */
__global__ void kharmonic_V(double *x, double *y, double *z, double *items,
                            double *Ax, double *Ay, double *Az, double *V);

 /**
 * @brief       Kernel to generate toroidal V (3d)
 */
__global__ void ktorus_V(double *x, double *y, double *z, double *items,
                         double *Ax, double *Ay, double *Az, double *V);

 /**
 * @brief       Kernel to generate simple gaussian wavefunction
 */
__global__ void kstd_wfc(double *x, double *y, double *z, double *items,
                         double winding, double *phi, double2 *wfc);

 /**
 * @brief       Kernel to generate toroidal wavefunction
 */
__global__ void ktorus_wfc(double *x, double *y, double *z, double *items,
                           double winding, double *phi, double2 *wfc);
 
 /**
 * @brief       Kernel to generate all auxiliary fields
 */

__global__ void aux_fields(double *V, double *K, double gdt, double dt,
                           double* Ax, double *Ay, double* Az,
                           double *px, double *py, double *pz,
                           double* pAx, double* pAy, double* pAz,
                           double2* GV, double2* EV, double2* GK, double2* EK,
                           double2* GpAx, double2* GpAy, double2* GpAz,
                           double2* EpAx, double2* EpAy, double2* EpAz);
// Function to generate grid and treads
void generate_grid(Grid &par);

#endif
