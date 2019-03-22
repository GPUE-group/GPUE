#include "../include/evolution.h"
#include "../include/vortex_3d.h"

// 3D
void apply_gauge(Grid &par, double2 *in, double2 *out, double2 *Ax, double2 *Ay,
                 double2 *Az, double renorm_factor_x,
                 double renorm_factor_y, double renorm_factor_z, bool flip,
                 cufftHandle plan_1d, cufftHandle plan_dim2,
                 cufftHandle plan_dim3, double dx, double dy, double dz,
                 double time, int yDim, int size){

    cufftResult result;
    dim3 grid = par.grid;
    dim3 threads = par.threads;

    if (flip){

        // 1d forward / mult by Az
        result = cufftExecZ2Z(plan_1d, in, out, CUFFT_FORWARD);
        scalarMult<<<grid,threads>>>(out, renorm_factor_z, out);
        if(par.bval("Az_time")){
            EqnNode_gpu* Az_eqn = par.astval("Az");
            int e_num = par.ival("Az_num");
            ast_cmult<<<grid,threads>>>(out, out, Az_eqn, dx, dy, dz,
                                        time, e_num);
        }
        else{
            cMult<<<grid,threads>>>(out, (cufftDoubleComplex*) Az, out);
        }
        result = cufftExecZ2Z(plan_1d, out, out, CUFFT_INVERSE);
        scalarMult<<<grid,threads>>>(out, renorm_factor_z, out);

        // loop to multiply by Ay
        for (int i = 0; i < yDim; i++){
            result = cufftExecZ2Z(plan_dim2,  &out[i*size],
                     &out[i*size], CUFFT_FORWARD);
        }

        scalarMult<<<grid,threads>>>(out, renorm_factor_y, out);
        if(par.bval("Ay_time")){
            EqnNode_gpu* Ay_eqn = par.astval("Ay");
            int e_num = par.ival("Ay_num");
            ast_cmult<<<grid,threads>>>(out, out, Ay_eqn, dx, dy, dz,
                                        time, e_num);
        }
        else{
            cMult<<<grid,threads>>>(out, (cufftDoubleComplex*) Ay, out);
        }

        for (int i = 0; i < yDim; i++){
            //size = xDim * zDim;
            result = cufftExecZ2Z(plan_dim2, &out[i*size],
                     &out[i*size], CUFFT_INVERSE);
        }
        scalarMult<<<grid,threads>>>(out, renorm_factor_y, out);

        // 1D FFT to Ax
        result = cufftExecZ2Z(plan_dim3, out, out, CUFFT_FORWARD);
        scalarMult<<<grid,threads>>>(out, renorm_factor_x, out);

        if(par.bval("Ax_time")){
            EqnNode_gpu* Ax_eqn = par.astval("Ax");
            int e_num = par.ival("Ax_num");
            ast_cmult<<<grid,threads>>>(out, out, Ax_eqn, dx, dy, dz,
                                        time, e_num);
        }
        else{
            cMult<<<grid,threads>>>(out, (cufftDoubleComplex*) Ax, out);
        }

        result = cufftExecZ2Z(plan_dim3, out, out, CUFFT_INVERSE);
        scalarMult<<<grid,threads>>>(out, renorm_factor_x, out);

    }
    else{

        // 1D FFT to Ax
        result = cufftExecZ2Z(plan_dim3, out, out, CUFFT_FORWARD);
        scalarMult<<<grid,threads>>>(out, renorm_factor_x, out);

        if(par.bval("Ax_time")){
            EqnNode_gpu* Ax_eqn = par.astval("Ax");
            int e_num = par.ival("Ax_num");
            ast_cmult<<<grid,threads>>>(out, out, Ax_eqn, dx, dy, dz,
                                        time, e_num);
        }
        else{
            cMult<<<grid,threads>>>(out, (cufftDoubleComplex*) Ax, out);
        }

        result = cufftExecZ2Z(plan_dim3, out, out, CUFFT_INVERSE);
        scalarMult<<<grid,threads>>>(out, renorm_factor_x, out);


        // loop to multiply by Ay
        for (int i = 0; i < yDim; i++){
            result = cufftExecZ2Z(plan_dim2,  &out[i*size],
                     &out[i*size], CUFFT_FORWARD);
        }

        scalarMult<<<grid,threads>>>(out, renorm_factor_y, out);
        if(par.bval("Ay_time")){
            EqnNode_gpu* Ay_eqn = par.astval("Ay");
            int e_num = par.ival("Ay_num");
            ast_cmult<<<grid,threads>>>(out, out, Ay_eqn, dx, dy, dz,
                                        time, e_num);
        }
        else{
            cMult<<<grid,threads>>>(out, (cufftDoubleComplex*) Ay, out);
        }

        for (int i = 0; i < yDim; i++){
            //size = xDim * zDim;
            result = cufftExecZ2Z(plan_dim2, &out[i*size],
                     &out[i*size], CUFFT_INVERSE);
        }
        scalarMult<<<grid,threads>>>(out, renorm_factor_y, out);

        // 1d forward / mult by Az
        result = cufftExecZ2Z(plan_1d, out, out, CUFFT_FORWARD);
        scalarMult<<<grid,threads>>>(out, renorm_factor_z, out);
        if(par.bval("Az_time")){
            EqnNode_gpu* Az_eqn = par.astval("Az");
            int e_num = par.ival("Az_num");
            ast_cmult<<<grid,threads>>>(out, out, Az_eqn, dx, dy, dz,
                                        time, e_num);
        }
        else{
            cMult<<<grid,threads>>>(out, (cufftDoubleComplex*) Az, out);
        }
        result = cufftExecZ2Z(plan_1d, out, out, CUFFT_INVERSE);
        scalarMult<<<grid,threads>>>(out, renorm_factor_z, out);

    }

}

// 2D
void apply_gauge(Grid &par, double2 *in, double2 *out, double2 *Ax, double2 *Ay,
                 double renorm_factor_x, double renorm_factor_y, bool flip,
                 cufftHandle plan_1d, cufftHandle plan_dim2, double dx,
                 double dy, double dz, double time){

    cufftResult result;
    dim3 grid = par.grid;
    dim3 threads = par.threads;

    if (flip){

        // 1d forward / mult by Ay
        result = cufftExecZ2Z(plan_1d, in, out, CUFFT_FORWARD);
        scalarMult<<<grid,threads>>>(out, renorm_factor_y ,out);
        if(par.bval("Ay_time")){
            EqnNode_gpu* Ay_eqn = par.astval("Ay");
            int e_num = par.ival("Ay_num");
            ast_cmult<<<grid,threads>>>(out, out, Ay_eqn, dx, dy, dz,
                                        time, e_num);
        }
        else{
            cMult<<<grid,threads>>>(out, (cufftDoubleComplex*) Ay, out);
        }
        result = cufftExecZ2Z(plan_1d, out, out, CUFFT_INVERSE);
        scalarMult<<<grid,threads>>>(out, renorm_factor_y, out);


        // 1D FFT to out_pAx
        result = cufftExecZ2Z(plan_dim2, out, out, CUFFT_FORWARD);
        scalarMult<<<grid,threads>>>(out, renorm_factor_x, out);
        if(par.bval("Ax_time")){
            EqnNode_gpu* Ax_eqn = par.astval("Ax");
            int e_num = par.ival("Ax_num");
            ast_cmult<<<grid,threads>>>(out, out, Ax_eqn, dx, dy, dz,
                                        time, e_num);
        }
        else{
            cMult<<<grid,threads>>>(out, (cufftDoubleComplex*) Ax, out);
        }

        result = cufftExecZ2Z(plan_dim2, out, out, CUFFT_INVERSE);
        scalarMult<<<grid,threads>>>(out, renorm_factor_x, out);
    }
    else{

        // 1D FFT to out_pAx
        result = cufftExecZ2Z(plan_dim2, out, out, CUFFT_FORWARD);
        scalarMult<<<grid,threads>>>(out, renorm_factor_x, out);
        if(par.bval("Ax_time")){
            EqnNode_gpu* Ax_eqn = par.astval("Ax");
            int e_num = par.ival("Ax_num");
            ast_cmult<<<grid,threads>>>(out, out, Ax_eqn, dx, dy, dz,
                                        time, e_num);
        }
        else{
            cMult<<<grid,threads>>>(out, (cufftDoubleComplex*) Ax, out);
        }

        result = cufftExecZ2Z(plan_dim2, out, out, CUFFT_INVERSE);
        scalarMult<<<grid,threads>>>(out, renorm_factor_x, out);


        // 1d forward / mult by Ay
        result = cufftExecZ2Z(plan_1d, out, out, CUFFT_FORWARD);
        scalarMult<<<grid,threads>>>(out, renorm_factor_y ,out);
        if(par.bval("Ay_time")){
            EqnNode_gpu* Ay_eqn = par.astval("Ay");
            int e_num = par.ival("Ay_num");
            ast_cmult<<<grid,threads>>>(out, out, Ay_eqn, dx, dy, dz,
                                        time, e_num);
        }
        else{
            cMult<<<grid,threads>>>(out, (cufftDoubleComplex*) Ay, out);
        }
        result = cufftExecZ2Z(plan_1d, out, out, CUFFT_INVERSE);
        scalarMult<<<grid,threads>>>(out, renorm_factor_y, out);

    }

}


void evolve(Grid &par,
            int numSteps,
            unsigned int gstate,
            std::string buffer){

    // Re-establishing variables from parsed Grid class
    std::string data_dir = par.sval("data_dir");
    int dimnum = par.ival("dimnum");
    double omega = par.dval("omega");
    double angle_sweep = par.dval("angle_sweep");
    double gdt = par.dval("gdt");
    double dt = par.dval("dt");
    double omegaX = par.dval("omegaX");
    double omegaY = par.dval("omegaY");
    double mass = par.dval("mass");
    double dx = par.dval("dx");
    double dy = 1;
    double dz = 1;
    double interaction = par.dval("interaction");
    double laser_power = par.dval("laser_power");
    double gDenConst = par.dval("gDenConst");
    double thresh_const = par.dval("thresh_const");
    double *x = par.dsval("x");
    double *y;
    double *V = par.dsval("V");
    double *Phi = par.dsval("Phi");
    double2 *gpu1dpAx = par.cufftDoubleComplexval("pAx_gpu");
    double2 *gpu1dpAy;
    double2 *gpu1dpAz;
    double *Phi_gpu = par.dsval("Phi_gpu");
    bool write_it = par.bval("write_it");
    bool graph = par.bval("graph");
    int N = par.ival("atoms");
    int printSteps = par.ival("printSteps");
    bool nonlin = par.bval("gpe");
    bool lz = par.bval("corotating");
    bool ramp = par.bval("ramp");
    int ramp_type = par.ival("ramp_type");
    int xDim = par.ival("xDim");
    int yDim = 1;
    int zDim = 1;
    cufftDoubleComplex *wfc = par.cufftDoubleComplexval("wfc");
    cufftDoubleComplex *gpuWfc = par.cufftDoubleComplexval("wfc_gpu");
    cufftDoubleComplex *K_gpu =
        par.cufftDoubleComplexval("K_gpu");
    cufftDoubleComplex *V_gpu =
        par.cufftDoubleComplexval("V_gpu");

    if (dimnum > 1){
        dy = par.dval("dy");
        y = par.dsval("y");
        gpu1dpAy = par.cufftDoubleComplexval("pAy_gpu");
        yDim = par.ival("yDim");
    }
    if (dimnum > 2){
        dz = par.dval("dz");
        gpu1dpAz = par.cufftDoubleComplexval("pAz_gpu");
        zDim = par.ival("zDim");
    }

    int gridSize = xDim * yDim * zDim;

    // getting data from Cuda class
    cufftResult result;
    cufftHandle plan_1d = par.ival("plan_1d");
    cufftHandle plan_2d = par.ival("plan_2d");
    cufftHandle plan_other2d = par.ival("plan_other2d");
    cufftHandle plan_3d = par.ival("plan_3d");
    cufftHandle plan_dim2 = par.ival("plan_dim2");
    cufftHandle plan_dim3 = par.ival("plan_dim3");
    dim3 threads = par.threads;
    dim3 grid = par.grid;

    // Because no two operations are created equally.
    // Multiplication is faster than divisions.
    double renorm_factor_nd=1.0/pow(gridSize,0.5);
    double renorm_factor_x=1.0/pow(xDim,0.5);
    double renorm_factor_y=1.0/pow(yDim,0.5);
    double renorm_factor_z=1.0/pow(zDim,0.5);

    clock_t begin, end;
    double time_spent;
    double Dt;
    if(gstate==0){
        Dt = gdt;
        printf("Timestep for groundstate solver set as: %E\n",Dt);
    }
    else{
        Dt = dt;
        printf("Timestep for evolution set as: %E\n",Dt);
    }
    begin = clock();
    double omega_0=omega*omegaX;

    // ** ############################################################## ** //
    // **         HERE BE DRAGONS OF THE MOST DANGEROUS KIND!            ** //
    // ** ############################################################## ** //

    // 2D VORTEX TRACKING

    double mask_2d = par.dval("mask_2d");
    int x0_shift = par.dval("x0_shift");
    int y0_shift = par.dval("y0_shift");
    int charge = par.ival("charge");
    int kill_idx = par.ival("kill_idx");
    cufftDoubleComplex *EV_opt = par.cufftDoubleComplexval("EV_opt");
    int kick_it = par.ival("kick_it");
    double *V_opt = par.dsval("V_opt");
    // Double buffering and will attempt to thread free and calloc operations to
    // hide time penalty. Or may not bother.
    int num_vortices[2] = {0,0};

    // binary matrix of size xDim*yDim,
    // 1 for vortex at specified index, 0 otherwise
    int* vortexLocation;
    //int* olMaxLocation = (int*) calloc(xDim*yDim,sizeof(int));

    std::shared_ptr<Vtx::Vortex> central_vortex; //vortex closest to the central position
    /*
    central_vortex.coords.x = -1;
    central_vortex.coords.y = -1;
    central_vortex.coordsD.x = -1.;
    central_vortex.coordsD.y = -1.;
    central_vortex.wind = 0;
    */

    // Angle of vortex lattice. Add to optical lattice for alignment.
    double vort_angle;

    // array of vortex coordinates from vortexLocation 1's
    //struct Vtx::Vortex *vortCoords = NULL;


    std::shared_ptr<Vtx::VtxList> vortCoords = std::make_shared<Vtx::VtxList>(7);
    //std::vector<std::shared_ptr<Vtx::Vortex> vortCoords;

    //Previous array of vortex coordinates from vortexLocation 1's
    //struct Vtx::Vortex *vortCoordsP = NULL;
    //std::vector<struct Vtx::Vortex> vortCoordsP;
    std::shared_ptr<Vtx::VtxList> vortCoordsP = std::make_shared<Vtx::VtxList>(7);


    LatticeGraph::Lattice lattice; //Vortex lattice graph.
    double* adjMat;

    // Assuming triangular lattice at rotatio


    //std::cout << "numSteps is: " << numSteps << '\n';
    // Iterating through all of the steps in either g or esteps.
    for(int i=0; i < numSteps; ++i){
        double time = Dt*i;
        if (ramp){

            //Adjusts omega for the appropriate trap frequency.
            if (ramp_type == 1){
                if (i == 0){
                    omega_0 = (double)omega;
                }
                else{
                    omega_0 = (double)i / (double)(i+1);
                }
            }
            else{
                if (i == 0){
                    omega_0=(double)omega/(double)(numSteps);
                }
                else{
                    omega_0 = (double)(i+1) / (double)i;
                }
            }
        }

        // Print-out at pre-determined rate.
        // Vortex & wfc analysis performed here also.
        if(i % printSteps == 0) {
            // If the unit_test flag is on, we need a special case
            printf("Step: %d    Omega: %lf\n", i, omega_0);
            cudaMemcpy(wfc, gpuWfc, sizeof(cufftDoubleComplex)*xDim*yDim*zDim, 
                       cudaMemcpyDeviceToHost);

            // Printing out time of iteration
            end = clock();
            time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
            printf("Time spent: %lf\n", time_spent);
            std::string fileName = "";
            //printf("ramp=%d        gstate=%d    rg=%d        \n",
            //       ramp, gstate, ramp | (gstate << 1));
            switch (ramp | (gstate << 1)) {
                case 0: //Groundstate solver, constant Omega value.
                {
                    fileName = "wfc_0_const";
                    break;
                }
                case 1: //Groundstate solver, ramped Omega value.
                {
                    fileName = "wfc_0_ramp";
                    break;
                }
                case 2: //Real-time evolution, constant Omega value.
                {
                    if (dimnum == 3){
                        // Note: In the case of 3d, we need to think about
                        //       vortex tracking in a new way.
                        //       It may be as simple as splitting the problem
                        //       into 2D elements and working from there, but
                        //       let's look into it when we need it in the
                        //       future.
                        std::cout << "commencing 3d vortex tracking" << '\n';

                        // Creating the necessary double* values
                        double* edges = (double *)malloc(sizeof(double)
                                                         *gridSize);


                        find_edges(par, wfc, edges);
                        double* edges_gpu = par.dsval("edges_gpu");

                        // Now we need to output everything
                        if (write_it){
                            FileIO::writeOutDouble(buffer, data_dir+"Edges",
                                                   edges, gridSize, i);
                        }
                        free(edges);

                    }
                    else if (dimnum == 2){
                        vortexLocation = (int *) calloc(xDim * yDim, sizeof(int));
                        num_vortices[0] = Tracker::findVortex(vortexLocation, wfc,
                                                              mask_2d, xDim, x, i);
                        // If initial step, locate vortices, least-squares to find
                        // exact centre, calculate lattice angle, generate optical
                        // lattice.
                        if (i == 0) {
                             if(num_vortices[0] > 0){
                             //Reserve enough space for the vortices
                             //reserve(num_vortices[0]);
                             vortCoords = std::make_shared<Vtx::VtxList>
                                              (num_vortices[0]);
                             vortCoordsP = std::make_shared<Vtx::VtxList>
                                               (num_vortices[0]);
    
                             //Locate the vortex positions to the nearest grid, then
                             //perform a least-squares fit to determine the location
                             //to sub-grid reolution.
                             Tracker::vortPos(vortexLocation,
                                     vortCoords->getVortices(), xDim, wfc);
                             Tracker::lsFit(vortCoords->getVortices(),wfc,xDim);
    
                             //Find the centre-most vortex in the lattice
                             central_vortex = Tracker::vortCentre(vortCoords->
                                     getVortices(), xDim);
                             //Determine the Angle formed by the lattice relative to
                             //the x-axis
                             vort_angle = Tracker::vortAngle(vortCoords->
                                     getVortices(), central_vortex);
    
                            //Store the vortex angle in the parameter file
                             par.store("Vort_angle", vort_angle);
                                 
                            //Determine average lattice spacing.
                             double sepAvg = Tracker::vortSepAvg(vortCoords->
                                     getVortices(), central_vortex);
    
                             par.store("Central_vort_x",
                                      (double) central_vortex->getCoords().x);
                             par.store("Central_vort_y",
                                      (double) central_vortex->getCoords().y);
                             par.store("Central_vort_winding",
                                      (double) central_vortex->getWinding());
                             par.store("Num_vort", (double) vortCoords->
                                     getVortices().size());
    
                            //Setup the optical lattice to match the spacing and
                            // angle+angle_sweep of the vortex lattice.
                            // Amplitude matched by setting laser_power
                            // parameter switch.
                             optLatSetup(central_vortex, V, 
                                        vortCoords->getVortices(),
                                        vort_angle + PI * angle_sweep / 180.0,
                                        laser_power * HBAR * sqrt(omegaX * omegaY),
                                        V_opt, x, y, par);
    
    
			    }
                            // If kick_it param is 2, perform a single kick of the 
                            // optical lattice for the first timestep only.
                            // This is performed by loading the
                            // EV_opt exp(V + V_opt) array into GPU memory
                            // for the potential.
                            if (kick_it == 2) {
                                printf("Kicked it 1\n");
                                cudaMemcpy(V_gpu, EV_opt,
                                           sizeof(cufftDoubleComplex) * xDim * yDim,
                                           cudaMemcpyHostToDevice);
                            }
                            // Write out the newly specified potential
                            // and exp potential to files
                            if(write_it){
                                FileIO::writeOutDouble(buffer, data_dir + "V_opt_1",
                                                       V_opt, xDim * yDim, 0);
                                FileIO::writeOut(buffer, data_dir + "EV_opt_1", EV_opt,
                                                 xDim * yDim, 0);
    
                                //Store necessary parameters to Params.dat file.
                                FileIO::writeOutParam(buffer, par,
                                                      data_dir + "Params.dat");
                            }
                        }
                        //If i!=0 and the number of vortices changes
                        // if num_vortices[1] < num_vortices[0] ... Fewer vortices
                        else {
                             if (num_vortices[0] > 0){
                        	    Tracker::vortPos(vortexLocation, 
                                        vortCoords->getVortices(), xDim, wfc);
                        	    Tracker::lsFit(vortCoords->getVortices(), 
                                                   wfc, xDim);
                        	    Tracker::vortArrange(vortCoords->getVortices(),
                                                        vortCoordsP->getVortices());
                                    if(write_it){
                    		        FileIO::writeOutInt(buffer, data_dir + "vLoc_",
                                                       vortexLocation, xDim * yDim, i);
                                    }
                             }
                        }
    
                        // The following will eventually be modified and moved into
                        // a new library that works closely wy0_shiftUE. Used to 
                        // also defined for vortex elimination using graph positions
                        // and UID numbers.
                        if (graph && num_vortices[0] > 0) {
                            for (int ii = 0; ii < vortCoords->getVortices().size();
                                 ++ii) {
                                std::shared_ptr<LatticeGraph::Node>
                                    n(new LatticeGraph::Node(
                                        *vortCoords->getVortices().at(ii).get()));
                                lattice.addVortex(std::move(n));
                            }
                            unsigned int *uids = (unsigned int *) malloc(
                                    sizeof(unsigned int) *
                                    lattice.getVortices().size());
                            for (size_t a=0; a < lattice.getVortices().size(); ++a){
                                uids[a] = lattice.getVortexIdx(a)->getUid();
                            }
                            if(i==0) {
                                //Lambda for vortex annihilation/creation.
                                auto killIt=[&](int idx, int winding, 
                                                double delta_x, double delta_y) {
                                    if (abs(delta_x) > 0 || abs(delta_y) > 0){
                                        // Killing initial vortex and then 
                                        // imprinting new one
                                        WFC::phaseWinding(Phi, 1, x,y, dx,dy,
                                            lattice.getVortexUid(idx)->
                                                getData().getCoordsD().x,
                                            lattice.getVortexUid(idx)->
                                                getData().getCoordsD().y,
                                            xDim);
    
                                        cudaMemcpy(Phi_gpu, Phi, 
                                                   sizeof(double) * xDim * yDim, 
                                                   cudaMemcpyHostToDevice);
                                        cMultPhi <<<grid, threads>>>(gpuWfc,Phi_gpu,
                                                                      gpuWfc);
    
                                        // Imprinting new one
                                        int cval = -winding;
                                        WFC::phaseWinding(Phi, cval, x,y, dx,dy,
                                            lattice.getVortexUid(idx)->
                                                getData().getCoordsD().x + delta_x,
                                            lattice.getVortexUid(idx)->
                                                getData().getCoordsD().y + delta_y,
                                            xDim);
    
                                        // Sending to device for imprinting
                                        cudaMemcpy(Phi_gpu, Phi, 
                                                   sizeof(double) * xDim * yDim, 
                                                   cudaMemcpyHostToDevice);
                                        cMultPhi <<<grid, threads>>>(gpuWfc,Phi_gpu,
                                                                      gpuWfc);
                                    }
                                    else{
                                        int cval = -(winding-1);
                                        WFC::phaseWinding(Phi, cval, x,y,dx,dy,
                                            lattice.getVortexUid(idx)->
                                                getData().getCoordsD().x,
                                            lattice.getVortexUid(idx)->
                                                getData().getCoordsD().y,
                                            xDim);
                                        cudaMemcpy(Phi_gpu, Phi, 
                                                   sizeof(double) * xDim * yDim, 
                                                   cudaMemcpyHostToDevice);
                                        cMultPhi <<<grid, threads>>>(gpuWfc,Phi_gpu,
                                                                      gpuWfc);
                                    }
                                };
                                if (kill_idx > 0){
                                    killIt(kill_idx, charge, x0_shift, y0_shift);
                                }
                            }
                            lattice.createEdges(1.5 * 2e-5 / dx);
    
                            //Assumes that vortices
                            //only form edges when the distance is upto 1.5*2e-5.
                            //Replace with delaunay triangulation determined edges 
                            //for better computational scaling (and sanity)
    
                            //O(n^2) -> terrible implementation. It works for now.
                            //Generates the adjacency matrix from the graph, and
                            //outputs to a Mathematica compatible format.
                            adjMat = (double *)calloc(lattice.getVortices().size() *
                                                      lattice.getVortices().size(),
                                                       sizeof(double));
                            lattice.genAdjMat(adjMat);
                            if (write_it){
                                FileIO::writeOutAdjMat(buffer, data_dir + "graph",
                                                       adjMat, uids,
                                                       lattice.getVortices().size(), i);
                           }
    
                            //Free and clear all memory blocks
                            free(adjMat);
                            free(uids);
                            lattice.getVortices().clear();
                            lattice.getEdges().clear();
                            //exit(0);
                        }

                        //Write out the vortex locations
                        if(write_it){
                            FileIO::writeOutVortex(buffer, data_dir + "vort_arr",
                                                   vortCoords->getVortices(), i);
                        }
                        printf("Located %d vortices\n", 
                               vortCoords->getVortices().size());
    
                        //Free memory block for now.
                        free(vortexLocation);

                        //Current values become previous values.
                        num_vortices[1] = num_vortices[0];
                        vortCoords->getVortices().swap(vortCoordsP->getVortices());
		                vortCoords->getVortices().clear();
    
                    }
                    fileName = "wfc_ev";
                    break;
                }
                case 3:
                {
                    fileName = "wfc_ev_ramp";
                    break;
                }
                default:
                {
                    break;
                }
            }

            //std::cout << "writing" << '\n';
            if (write_it) {
                FileIO::writeOut(buffer, data_dir + fileName,
                                 wfc, xDim*yDim*zDim, i);
            }
            //std::cout << "written" << '\n';
            if (par.bval("energy_calc")){
                par.store("time", time);
                double energy = energy_calc(par,gpuWfc);
                // Now opening and closing file for writing.
                std::ofstream energy_out;
                std::string mode = "energyi.dat";
                if (gstate == 1){
                    mode = "energy.dat";
                }
                if (i == 0){
                    energy_out.open(data_dir + mode);
                }
                else{
                    energy_out.open(data_dir + mode, std::ios::out |
                                                     std::ios::app);
                }
                energy_out << energy << '\n';
                energy_out.close();
                printf("Energy[t@%d]=%E\n",i,energy);
                par.store("energy", energy);
            }

        }

        // No longer writing out

        // ** ########################################################## ** //
        // **                     More F'n' Dragons!                     ** //
        // ** ########################################################## ** //

        // U_r(dt/2)*wfc
        if(nonlin == 1){
            if(par.bval("V_time")){
                EqnNode_gpu* V_eqn = par.astval("V");
                int e_num = par.ival("V_num");
                cMultDensity_ast<<<grid,threads>>>(V_eqn,gpuWfc,gpuWfc,
                    dx, dy, dz, time, e_num, 0.5*Dt,
                    gstate,interaction*gDenConst);
            }
            else{
                cMultDensity<<<grid,threads>>>(V_gpu,gpuWfc,gpuWfc,0.5*Dt,
                    gstate,interaction*gDenConst);
            }
        }
        else {
            if(par.bval("V_time")){ 
                EqnNode_gpu* V_eqn = par.astval("V");
                int e_num = par.ival("V_num");
                ast_op_mult<<<grid,threads>>>(gpuWfc,gpuWfc, V_eqn,
                    dx, dy, dz, time, e_num, gstate+1, Dt);
            }
            else{
                cMult<<<grid,threads>>>(V_gpu,gpuWfc,gpuWfc);
            }
        }

        // U_p(dt)*fft2(wfc)
        result = cufftExecZ2Z(plan_3d,gpuWfc,gpuWfc,CUFFT_FORWARD);

        // Normalise
        scalarMult<<<grid,threads>>>(gpuWfc,renorm_factor_nd,gpuWfc);
        if (par.bval("K_time")){
            EqnNode_gpu* k_eqn = par.astval("k");
            int e_num = par.ival("k_num");
            ast_op_mult<<<grid,threads>>>(gpuWfc,gpuWfc, k_eqn,
                dx, dy, dz, time, e_num, gstate+1, Dt);
        }
        else{
            cMult<<<grid,threads>>>(K_gpu,gpuWfc,gpuWfc);
        }
        result = cufftExecZ2Z(plan_3d,gpuWfc,gpuWfc,CUFFT_INVERSE);

        // Normalise
        scalarMult<<<grid,threads>>>(gpuWfc,renorm_factor_nd,gpuWfc);

        // U_r(dt/2)*wfc
        if(nonlin == 1){
            if(par.bval("V_time")){
                EqnNode_gpu* V_eqn = par.astval("V");
                int e_num = par.ival("V_num");
                cMultDensity_ast<<<grid,threads>>>(V_eqn,gpuWfc,gpuWfc,
                    dx, dy, dz, time, e_num, 0.5*Dt,
                    gstate,interaction*gDenConst);
            }
            else{
                cMultDensity<<<grid,threads>>>(V_gpu,gpuWfc,gpuWfc,0.5*Dt,
                    gstate,interaction*gDenConst);
            }
        }
        else {
            if(par.bval("V_time")){  
                EqnNode_gpu* V_eqn = par.astval("V");
                int e_num = par.ival("V_num");
                ast_op_mult<<<grid,threads>>>(gpuWfc,gpuWfc, V_eqn,
                    dx, dy, dz, time, e_num, gstate+1, Dt);
            }
            else{
                cMult<<<grid,threads>>>(V_gpu,gpuWfc,gpuWfc);
            }

        }

        // Angular momentum pAy-pAx (if engaged)  //
        if(lz == true){
            // Multiplying by ramping factor if necessary
            // Note: using scalarPow to do the scaling inside of the exp
            if (ramp){
                scalarPow<<<grid,threads>>>((cufftDoubleComplex*) gpu1dpAy, 
                                            omega_0,
                                            (cufftDoubleComplex*) gpu1dpAy);
                if (dimnum > 1){
                    scalarPow<<<grid,threads>>>((cufftDoubleComplex*) gpu1dpAx, 
                                                omega_0,
                                                (cufftDoubleComplex*) gpu1dpAx);
                }
                if (dimnum > 2){
                    scalarPow<<<grid,threads>>>((cufftDoubleComplex*) gpu1dpAz, 
                                                omega_0,
                                                (cufftDoubleComplex*) gpu1dpAz);
                }
            }
            int size = xDim*zDim;
            if (dimnum == 3){
                apply_gauge(par, gpuWfc, gpuWfc, gpu1dpAx, gpu1dpAy, gpu1dpAz,
                            renorm_factor_x, renorm_factor_y, renorm_factor_z,
                            i%2, plan_1d, plan_dim2, plan_dim3,
                            dx, dy, dz, time, yDim, size);
            }
            else if (dimnum == 2){
                apply_gauge(par, gpuWfc, gpuWfc, gpu1dpAx, gpu1dpAy,
                            renorm_factor_x, renorm_factor_y, i%2, plan_1d,
                            plan_other2d, dx, dy, dz, time);
            }
            else if (dimnum == 1){
                result = cufftExecZ2Z(plan_1d,gpuWfc,gpuWfc,CUFFT_FORWARD);
                scalarMult<<<grid,threads>>>(gpuWfc,renorm_factor_x,gpuWfc);
                if(par.bval("Ax_time")){
                    EqnNode_gpu* Ax_eqn = par.astval("Ax");
                    int e_num = par.ival("Ax_num");
                    ast_cmult<<<grid,threads>>>(gpuWfc, gpuWfc,
                        Ax_eqn, dx, dy, dz, time, e_num);
                }
                else{
                    cMult<<<grid,threads>>>(gpuWfc,
                        (cufftDoubleComplex*) gpu1dpAx, gpuWfc);
                }

                result = cufftExecZ2Z(plan_1d,gpuWfc,gpuWfc, CUFFT_INVERSE);
                scalarMult<<<grid,threads>>>(gpuWfc, renorm_factor_x, gpuWfc);
            }
        }

        if(gstate==0){
            parSum(gpuWfc, par);
        }
    }

    par.store("wfc", wfc);
    par.store("wfc_gpu", gpuWfc);
}
