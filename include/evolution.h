///@endcond
//##############################################################################
/**
 *  @file    evolution.h
 *  @author  James Ryan Schloss (leios)
 *  @date    5/31/2016
 *  @version 0.1
 *
 *  @brief function for evolution.
 *
 *  @section DESCRIPTION
 *  These functions and variables are necessary for carrying out the GPUE
 *	simulations. This file will be re-written in an improved form in some
 *	future release.
 */
//##############################################################################

#ifndef EVOLUTION_H
#define EVOLUTION_H

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <ctime>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include <ctype.h>
#include <getopt.h>
#include "tracker.h"
#include "ds.h"
#include "split_op.h"
#include "kernels.h"
#include "constants.h"
#include "fileIO.h"
#include "lattice.h"
#include "manip.h"
#include "unit_test.h"
#include "vortex_3d.h"



// UPDATE LIST LATER
 /**
 * @brief       performs real or imaginary time evolution
 * @ingroup     data
 * @param       Parameter set
 * @param       Total number of steps
 * @param       Real (1) or imaginary (1) time evolution
 * @param       String buffer for writing files
 * @return      0 for success. See CUDA failure codes in cuda.h for other values
 */
void evolve(Grid &par,
            int numSteps,
            unsigned int gstate,
            std::string buffer);

void apply_gauge(Grid &par, double2 *wfc, double2 *Ax, double2 *Ay,
                 double2 *Az, double renorm_factor_x,
                 double renorm_factor_y, double renorm_factor_z, bool flip,
                 cufftHandle plan_1d, cufftHandle plan_dim2,
                 cufftHandle plan_dim3, double dx, double dy, double dz,
                 double time, int yDim, int size);

void apply_gauge(Grid &par, double2 *wfc, double2 *Ax, double2 *Ay,
                 double renorm_factor_x, double renorm_factor_y, bool flip,
                 cufftHandle plan_1d, cufftHandle plan_dim2, double dx,
                 double dy, double dz, double time);

#endif
