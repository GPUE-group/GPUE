///@endcond
//##############################################################################
/**
 *  @file    manip.cu
 *  @author  Lee J. O'Riordan (mlxd)
 *  @date    11/08/2015
 *  @version 0.1
 */
 //##############################################################################
#include "manip.h"

void WFC::phaseWinding(double *phi, int winding, double *x, double *y, double dx, double dy, double posx, double posy, int dim){
    for(int ii=0; ii < dim; ++ii){
        for(int jj=0; jj < dim; ++jj){
            phi[ii*dim +jj] = fmod(winding*atan2(y[jj]-(posy-dim/2+1)*dx, x[ii]-(posx-dim/2+1)*dy),2*PI);
        }
    }
}

void WFC::phaseWinding(double *phi, int winding, double *x, double *y, double dx, double dy, double* posx, double* posy, int sites, int dim){
    memset(phi,0,dim*dim);
    for(int a = 0; a< sites; ++a){
        for(int ii=0; ii < dim; ++ii){
            for(int jj=0; jj < dim; ++jj){
                phi[ii*dim +jj] = fmod(phi[ii*dim +jj] 
                    + (double)winding*(atan2(y[jj]-(posy[a]-dim/2+1)*dx, 
                                     x[ii]-(posx[a]-dim/2+1)*dy) ),2*PI);
            }
        }
    }
}

void WFC::applyPhase(double *phi, double2 *wfc_array, int dim){
    for(int ii=dim-1; ii >= 0; --ii){
        for(int jj=dim-1; jj >= 0; --jj){
            wfc_array[ii*dim + jj].x *= cos(phi[ii*dim + jj]);
            wfc_array[ii*dim + jj].y *= -sin(phi[ii*dim + jj]);
        }
    }
}
