
#include "../include/tracker.h"
#include "../include/fileIO.h"
#include "../include/minions.h"
#include "../include/constants.h"
#include "../include/vort.h"

/**
 *  Contains all the glorious info you need to track vortices and see what they are up to.
 **/
namespace Tracker {
    char bufferT[1024];

    /**
     * Determines the vortex separation at the centre of the lattice.
     */
    double vortSepAvg(const std::vector<std::shared_ptr<Vtx::Vortex> > &vArray, const std::shared_ptr<Vtx::Vortex> centre){
        double min = 0.0;
        double min_tmp = 0.0;
        min = sqrt( pow(centre->getCoordsD().x - vArray[0]->getCoordsD().x,2) + pow(centre->getCoordsD().y - vArray[0]->getCoordsD().y,2));
        for (int j=1; j<vArray.size(); ++j){
            min_tmp    = sqrt( pow(centre->getCoordsD().x - vArray[j]->getCoordsD().x,2) + pow(centre->getCoordsD().y - vArray[j]->getCoordsD().y,2));
            if(min > min_tmp && min_tmp > 1e-7){//100nm length
                min = min_tmp;
            }
        }
        return min;
    }

    /**
     * Finds the maxima of the optical lattice. Deprecated.
     */
    [[deprecated]]
    int findOLMaxima(int *marker, double *Vopt, double radius, int xDim, double* x){
        double gridValues[9];
        int i,j,found;
        found=0;
        for (i=1; i<xDim-1; ++i ){
            for(j=1; j<xDim-1;++j){
                if(sqrt(x[i]*x[i] + x[j]*x[j]) < radius){
                    gridValues[0] = Vopt[(i-1)*xDim + (j-1)];
                    gridValues[1] = Vopt[(i-1)*xDim + j];
                    gridValues[2] = Vopt[(i-1)*xDim + (j+1)];
                    gridValues[3] = Vopt[i*xDim + (j-1)];
                    gridValues[4] = Vopt[i*xDim + j];
                    gridValues[5] = Vopt[i*xDim + (j+1)];
                    gridValues[6] = Vopt[(i+1)*xDim + (j-1)];
                    gridValues[7] = Vopt[(i+1)*xDim + j];
                    gridValues[8] = Vopt[(i+1)*xDim + (j+1)];
                    if(fabs((gridValues[4]-Minions::maxValue(gridValues,9))/gridValues[4]) <= 1e-7){
                        //printf ("%d,%d\n",i,j);
                        (marker)[i*xDim + j] = 1;
                        ++found;
                    }
                }
            }
        }
        return found;
    }

    #ifdef VORT_MIN
    [[deprecated]]
    int findVortex(int *marker, const double2* wfc, double radius, int xDim,
                   const double* x, int timestep){

        double gridValues[9];
        int2 vIndex[1024];
        int2 index;
        int i,j,found;
        found=0;
    //    #pragma omp parallel for private(j)
        for (i=1; i<xDim-1; ++i ){
            for(j=1; j<xDim-1;++j){
                if(sqrt(x[i]*x[i] + x[j]*x[j]) < radius){
                    gridValues[0] = Minions::psi2(wfc[(i-1)*xDim + (j-1)]);
                    gridValues[1] = Minions::psi2(wfc[(i-1)*xDim + j]);
                    gridValues[2] = Minions::psi2(wfc[(i-1)*xDim + (j+1)]);
                    gridValues[3] = Minions::psi2(wfc[(i)*xDim + (j-1)]);
                    gridValues[4] = Minions::psi2(wfc[(i)*xDim + j]);
                    gridValues[5] = Minions::psi2(wfc[(i)*xDim + (j+1)]);
                    gridValues[6] = Minions::psi2(wfc[(i+1)*xDim + (j-1)]);
                    gridValues[7] = Minions::psi2(wfc[(i+1)*xDim + j]);
                    gridValues[8] = Minions::psi2(wfc[(i+1)*xDim + (j+1)]);
                    if(fabs((gridValues[4]-Minions::minValue(gridValues,9)) /
                             gridValues[4]) < 1e-7){
                        //printf ("%d,%d\n",i,j);
                        (marker)[i*xDim + j] = 1;
                        index.x=i;
                        index.y=j;
                        vIndex[found] = index;
                        found++;
                    }
                }
            }
        }
        return found;
    }
    #else
    /**
     * Phase winding method to determine vortex positions. Calculates the phase around a loop and checks if ~ +/-2Pi.
     */
    int findVortex(int *marker, const double2* wfc, double radius, int xDim, const double *x, int timestep){
            double2 *g = (double2*) malloc(sizeof(double2)*4);
            double *phiDelta = (double*) malloc(sizeof(double)*4);
        int i,j,found;
        int cond_x, cond_y;
        cond_x = 0; cond_y = 0;
        found = 0;
        long rnd_value = 0;
        double sum = 0.0;
        for ( i=0; i < xDim-1; ++i ){
            for( j=0; j < xDim-1; ++j ){
                if(sqrt(x[i]*x[i] + x[j]*x[j]) < radius){
                    g[0] = Minions::complexScale(
                        Minions::complexDiv( wfc[i*xDim + j],
                        wfc[(i+1)*xDim + j] ),
                        (Minions::complexMag(wfc[(i+1)*xDim+j])
                        / Minions::complexMag(wfc[i*xDim+j])));
                    g[1] = Minions::complexScale(
                        Minions::complexDiv(wfc[(i+1)*xDim+j],
                        wfc[(i+1)*xDim + (j+1)] ),
                        (Minions::complexMag(wfc[(i+1)*xDim+(j+1)])
                        / Minions::complexMag( wfc[(i+1)*xDim + j] )));
                    g[2] = Minions::complexScale(
                        Minions::complexDiv( wfc[(i+1)*xDim + (j+1)],
                        wfc[i*xDim + (j+1)] ),
                        ( Minions::complexMag( wfc[i*xDim + (j+1)])
                        / Minions::complexMag( wfc[(i+1)*xDim + (j+1)] )));
                    g[3] = Minions::complexScale(
                        Minions::complexDiv( wfc[i*xDim + (j+1)],
                        wfc[i*xDim + j] ),
                       ( Minions::complexMag( wfc[i*xDim + j])
                       / Minions::complexMag( wfc[i*xDim + (j+1)] )));
                    for (int k=0; k<4; ++k){
                        phiDelta[k] = atan2( g[k].y, g[k].x );
                        if(phiDelta[k] <= -PI){
                            phiDelta[k] += 2*PI;
                        }
                    }
                    sum = phiDelta[0] + phiDelta[1] + phiDelta[2] + phiDelta[3];
                    rnd_value = lround(sum/(2*PI));
                    if( sum >= 1.9*PI && cond_x <= 0 && cond_y <= 0){
                        marker[i*xDim + j] = rnd_value;
                        ++found;
                        sum = 0.0;
                        cond_x = 2; cond_y = 2;
                    }
                    else if( sum <= -1.9*PI && cond_x <= 0 && cond_y <= 0 )  {
                        marker[i*xDim + j] = -rnd_value;
                        ++found;
                        sum = 0.0;
                        cond_x = 2; cond_y = 2;
                    }
                --cond_x;
                --cond_y;
                }
            }
        }
        free(g); free(phiDelta);
        return found;
    }
    #endif

    /**
     * Accepts matrix of vortex locations as argument, returns array of x,y coordinates of locations and first encountered vortex angle
     */
     [[deprecated]]
    void olPos(int *marker, int2 *olLocation, int xDim){
        int i,j;
        unsigned int counter=0;
        for(i=0; i<xDim; ++i){
            for(j=0; j<xDim; ++j){
                if((marker)[i*xDim + j] == 1){
                    (olLocation)[ counter ].x=i;
                    (olLocation)[ counter ].y=j;
                    ++counter;
                }
            }
        }
    }

    /**
     * Tests the phase winding of the wavefunction, looking for vortices
     */
    int phaseTest(int2 vLoc, const double2* wfc, int xDim){
        int result = 0;
        double2 gridValues[4];
        double phiDelta[4];
        double sum=0.0;
        int i=vLoc.x, j=vLoc.y;
        gridValues[0] = Minions::complexScale( Minions::complexDiv(wfc[i*xDim + j],wfc[(i+1)*xDim + j]),             (Minions::complexMag(wfc[(i+1)*xDim + j])     / Minions::complexMag(wfc[i*xDim + j])));
        gridValues[1] = Minions::complexScale( Minions::complexDiv(wfc[(i+1)*xDim + j],wfc[(i+1)*xDim + (j+1)]),     (Minions::complexMag(wfc[(i+1)*xDim + (j+1)])/ Minions::complexMag(wfc[(i+1)*xDim + j])));
        gridValues[2] = Minions::complexScale( Minions::complexDiv(wfc[(i+1)*xDim + (j+1)],wfc[i*xDim + (j+1)]),     (Minions::complexMag(wfc[i*xDim + (j+1)])     / Minions::complexMag(wfc[(i+1)*xDim + (j+1)])));
        gridValues[3] = Minions::complexScale( Minions::complexDiv(wfc[i*xDim + (j+1)],wfc[i*xDim + j]),             (Minions::complexMag(wfc[i*xDim + j])         / Minions::complexMag(wfc[i*xDim + (j+1)])));

        for (int k=0; k<4; ++k){
            phiDelta[k] = atan2(gridValues[k].y,gridValues[k].x);
                    if(phiDelta[k] <= -PI){
                        phiDelta[k] += 2*PI;
            }
        }
        sum = phiDelta[0] + phiDelta[1] + phiDelta[2] + phiDelta[3];
        if(sum >=1.8*PI){
            result = 1;
        }
        return result;
    }

    /**
     * Accepts matrix of vortex locations as argument, returns array of x,y coordinates of locations and first encountered vortex angle
     */
    void vortPos(const int *marker, std::vector<std::shared_ptr<Vtx::Vortex> > &vLocation, int xDim, const double2 *wfc){
        int i,j;
        // unsigned int counter=0;

        int2 coords; double2 coordsD;
        coords.x=0; coords.y = 0;
        coordsD.x=0.; coordsD.y = 0.;

        for( i = 0; i < xDim; ++i){
            for( j = 0; j < xDim; ++j){
                if( abs((marker)[i*xDim + j]) >= 1){
                    coords.x = i; coords.y = j;
                    vLocation.push_back(std::make_shared<Vtx::Vortex>(coords, coordsD, marker[i*xDim + j], false, 0));
                    // ++counter;
                }
            }
        }
    }

    /**
     * Ensures the vortices are tracked and arranged in the right order based on minimum distance between previous and current positions
     */
    void vortArrange(std::vector<std::shared_ptr<Vtx::Vortex> > &vCoordsC, const std::vector<std::shared_ptr<Vtx::Vortex> > &vCoordsP){
        int dist, dist_t;
        int i, j, index;
        for ( i = 0; i < vCoordsC.size(); ++i ){
            dist = 0x7FFFFFFF; //arbitrary big value fo initial distance value
            index = i;
            for ( j = i; j < vCoordsC.size() ; ++j){//Changed to C and P from num_vort[0] size for both. May be an issue here with inconsistent sizing
                dist_t = ( (vCoordsP[i]->getCoordsD().x 
                            - vCoordsC[j]->getCoordsD().x)
                         * (vCoordsP[i]->getCoordsD().x 
                            - vCoordsC[j]->getCoordsD().x) 
                         + (vCoordsP[i]->getCoordsD().y 
                            - vCoordsC[j]->getCoordsD().y)
                         * (vCoordsP[i]->getCoordsD().y 
                            - vCoordsC[j]->getCoordsD().y) );
                if(dist > dist_t ){
                    dist = dist_t;
                    index = j;
                }
            }
            std::swap(vCoordsC[index], vCoordsC[i]); // Swap the elements at the given positions. Remove call to Minions::coordSwap(vCoordsC,index,i);
        }
    }

    /**
     * Determines the coords of the vortex closest to the central position. Useful for centering the optical lattice over v. lattice*
    */
    std::shared_ptr<Vtx::Vortex> vortCentre(const std::vector<std::shared_ptr<Vtx::Vortex> > &cArray, int xDim){
        int i, counter=0;
        int valX, valY;
        double valueTest, value = 0.0;
        valX = (cArray)[0]->getCoordsD().x - ((xDim/2)-1);
        valY = (cArray)[0]->getCoordsD().y - ((xDim/2)-1);
        value = sqrt( valX*valX + valY*valY );//Calcs the sqrt(x^2+y^2) from central position. try to minimise this value
        for ( i=1; i<cArray.size(); ++i ){
            valX = (cArray)[i]->getCoordsD().x - ((xDim/2)-1);
            valY = (cArray)[i]->getCoordsD().y - ((xDim/2)-1);
            valueTest = sqrt(valX*valX + valY*valY);
            if(value > valueTest){
                value = valueTest;
                counter = i;
            }
        }
        return (cArray)[counter];
    }

    /**
     * Determines the angle of the vortex lattice relative to the x-axis
     */
    double vortAngle(const std::vector<std::shared_ptr<Vtx::Vortex>> &vortCoords, const std::shared_ptr<Vtx::Vortex> central){
        int location = 0;
        double minVal=1e300;//(pow(central.x - vortCoords[0].x,2) + pow(central.y - vortCoords[0].y,2));
        for (int i=0; i < vortCoords.size(); ++i){//Assuming healing length on the order of 2 um
            if (minVal > (pow(central->getCoordsD().x - vortCoords[i]->getCoordsD().x,2) + pow(central->getCoordsD().y - vortCoords[i]->getCoordsD().y,2)) && std::abs(central->getCoordsD().x - vortCoords[i]->getCoordsD().x) > 2e-6 && std::abs(central->getCoordsD().y - vortCoords[i]->getCoordsD().y) > 2e-6){
                minVal = (pow(central->getCoordsD().x - vortCoords[i]->getCoordsD().x,2) + pow(central->getCoordsD().y - vortCoords[i]->getCoordsD().y,2));
                location = i;
            }
        }
        double ang=(fmod(atan2( (vortCoords[location]->getCoordsD().y - central->getCoordsD().y), (vortCoords[location]->getCoordsD().x - central->getCoordsD().x) ),PI/3));
        printf("Angle=%e\n",ang);
        return PI/3 - ang;

        //return PI/2 + fmod(atan2(vortCoords[location].y-central.y, vortCoords[location].x - central.x), PI/3);
        //return PI/2 - sign*acos( ( (central.x - vortCoords[location].x)*(central.x - vortCoords[location].x) ) / ( minVal*(central.x - vortCoords[location].x) ) );
    }

    /**
     * Sigma of vortex lattice and optical lattice. Deprecated
     */
    [[deprecated]]
    double sigVOL(const std::vector<std::shared_ptr<Vtx::Vortex> > &vArr, const int2 *opLatt, const double *x){
        double sigma = 0.0;
        double dx = std::abs(x[1]-x[0]);
        for (int i=0; i<vArr.size(); ++i){
            sigma += pow( std::abs( sqrt( (vArr[i]->getCoordsD().x - opLatt[i].x)*(vArr[i]->getCoordsD().x - opLatt[i].x) + (vArr[i]->getCoordsD().y - opLatt[i].y)*(vArr[i]->getCoordsD().y - opLatt[i].y) )*dx),2);
        }
        sigma /= vArr.size();
        return sigma;
    }

    /**
     * Performs least squares fitting to get exact vortex core position.
     */
    void lsFit(std::vector<std::shared_ptr<Vtx::Vortex> > &vortCoords, const double2 *wfc, int xDim){
        double2 *wfc_grid = (double2*) malloc(sizeof(double2)*4);
        double2 *res = (double2*) malloc(sizeof(double2)*3);
        double2 R;
        double2 coordsAdjusted;
        double det=0.0;
        for(int ii=0; ii < vortCoords.size(); ++ii){
            //vortCoords[ii]->getCoordsD().x = 0.0; vortCoords[ii]->getCoordsD().y = 0.0;
            coordsAdjusted.x=0.; coordsAdjusted.y=0.;

            wfc_grid[0] = wfc[vortCoords[ii]->getCoords().x*xDim + vortCoords[ii]->getCoords().y];
            wfc_grid[1] = wfc[(vortCoords[ii]->getCoords().x + 1)*xDim + vortCoords[ii]->getCoords().y];
            wfc_grid[2] = wfc[vortCoords[ii]->getCoords().x*xDim + (vortCoords[ii]->getCoords().y + 1)];
            wfc_grid[3] = wfc[(vortCoords[ii]->getCoords().x + 1)*xDim + (vortCoords[ii]->getCoords().y + 1)];

            for(int jj=0; jj<3; ++jj) {
                res[jj].x = lsq[jj][0]*wfc_grid[0].x + lsq[jj][1]*wfc_grid[1].x + lsq[jj][2]*wfc_grid[2].x + lsq[jj][3]*wfc_grid[3].x;
                res[jj].y = lsq[jj][0]*wfc_grid[0].y + lsq[jj][1]*wfc_grid[1].y + lsq[jj][2]*wfc_grid[2].y + lsq[jj][3]*wfc_grid[3].y;
            }

            //Solve Ax=b here. A = res[0,1], b = - res[2]. Solution -> X
            det = 1.0/(res[0].x*res[1].y - res[0].y*res[1].x);
            R.x = det*(res[1].y*res[2].x - res[0].y*res[2].y);
            R.y = det*(-res[1].x*res[2].x + res[0].x*res[2].y);
            coordsAdjusted.x = vortCoords[ii]->getCoords().x - R.x;
            coordsAdjusted.y = vortCoords[ii]->getCoords().y - R.y;
            vortCoords[ii]->updateCoordsD(coordsAdjusted);
        }
    }

    void updateVortices(std::shared_ptr<Vtx::VtxList> vLCurrent, std::shared_ptr<Vtx::VtxList> &vLPrev){
        //Iterate through the previous vortices, and compare the distances
        for (auto vtxPrev : vLPrev->getVortices()) {
            auto vtxMin = vLCurrent->minDistPair(vtxPrev, 1);
            //If the vortex has a min-pairing, and has a UID < 0 then we can assign it the same UID as the paired vortex
            if (vtxMin.second != nullptr && vtxMin.second->getUID() < 0 ) {
                vtxMin.second->updateUID(vtxPrev->getUID());
                vtxMin.second->updateIsOn(true);
            }
                //If no pairing found, then the vortex has disappeared or been killed. Switch it off, and add it to the current list with the given UID
            else{
                vtxPrev->updateIsOn(false);
                vLCurrent->addVtx(vtxPrev);//Will this cause trouble? Maybe rethink the UID determination
            }
        }
        //Find new vortices, assign them UIDs and switch them on
        for (auto v: vLCurrent->getVortices()) {
            if (v->getUID() < 0){
                v->updateUID(vLCurrent->getMax_Uid()++);
                v->updateIsOn(true);
            }
        }
        //Sort the list based on vortex UIDS. This may not be necessary, but helps for now with debugging things
        std::sort(
                vLCurrent->getVortices().begin(),
                vLCurrent->getVortices().end(),
                []( std::shared_ptr<Vtx::Vortex> a,
                    std::shared_ptr<Vtx::Vortex> b) {
                        return b->getUID() < a->getUID();
        });
        //Overwrite previous list with current list.
        vLPrev->getVortices().swap(vLCurrent->getVortices());
    }
}
