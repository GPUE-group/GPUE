#include "evolution.h"
#include "vortex_3d.h"

// 3D
void apply_gauge(Grid &par, double2 *wfc, double2 *Ax, double2 *Ay,
                 double2 *Az, double renorm_factor_x,
                 double renorm_factor_y, double renorm_factor_z, bool flip,
                 cufftHandle plan_1d, cufftHandle plan_dim2,
                 cufftHandle plan_dim3, double dx, double dy, double dz,
                 double time, int yDim, int size){

    dim3 grid = par.grid;
    dim3 threads = par.threads;

    if (flip){

        // 1d forward / mult by Ax
        cufftHandleError( cufftExecZ2Z(plan_1d, wfc, wfc, CUFFT_FORWARD) );
        scalarMult<<<grid,threads>>>(wfc, renorm_factor_x, wfc);
        cudaCheckError();
        if(par.bval("Ax_time")){
            EqnNode_gpu* Ax_eqn = par.astval("Ax");
            int e_num = par.ival("Ax_num");
            ast_cmult<<<grid,threads>>>(wfc, wfc, Ax_eqn, dx, dy, dz,
                                        time, e_num);
            cudaCheckError();
        }
        else{
            cMult<<<grid,threads>>>(wfc, (cufftDoubleComplex*) Ax, wfc);
            cudaCheckError();
        }
        cufftHandleError( cufftExecZ2Z(plan_1d, wfc, wfc, CUFFT_INVERSE) );
        scalarMult<<<grid,threads>>>(wfc, renorm_factor_x, wfc);
        cudaCheckError();

        // loop to multiply by Ay
        for (int i = 0; i < yDim; i++){
            cufftHandleError(cufftExecZ2Z(plan_dim2,  &wfc[i*size],
                             &wfc[i*size], CUFFT_FORWARD));
        }

        scalarMult<<<grid,threads>>>(wfc, renorm_factor_y, wfc);
        cudaCheckError();
        if(par.bval("Ay_time")){
            EqnNode_gpu* Ay_eqn = par.astval("Ay");
            int e_num = par.ival("Ay_num");
            ast_cmult<<<grid,threads>>>(wfc, wfc, Ay_eqn, dx, dy, dz,
                                        time, e_num);
            cudaCheckError();
        }
        else{
            cMult<<<grid,threads>>>(wfc, (cufftDoubleComplex*) Ay, wfc);
            cudaCheckError();
        }

        for (int i = 0; i < yDim; i++){
            //size = xDim * zDim;
            cufftHandleError(cufftExecZ2Z(plan_dim2, &wfc[i*size],
                             &wfc[i*size], CUFFT_INVERSE));
        }
        scalarMult<<<grid,threads>>>(wfc, renorm_factor_y, wfc);
        cudaCheckError();

        // 1D FFT to Az
        cufftHandleError( cufftExecZ2Z(plan_dim3, wfc, wfc, CUFFT_FORWARD) );
        scalarMult<<<grid,threads>>>(wfc, renorm_factor_z, wfc);
        cudaCheckError();

        if(par.bval("Az_time")){
            EqnNode_gpu* Az_eqn = par.astval("Az");
            int e_num = par.ival("Az_num");
            ast_cmult<<<grid,threads>>>(wfc, wfc, Az_eqn, dx, dy, dz,
                                        time, e_num);
            cudaCheckError();
        }
        else{
            cMult<<<grid,threads>>>(wfc, (cufftDoubleComplex*) Az, wfc);
            cudaCheckError();
        }

        cufftHandleError( cufftExecZ2Z(plan_dim3, wfc, wfc, CUFFT_INVERSE) );
        scalarMult<<<grid,threads>>>(wfc, renorm_factor_z, wfc);
        cudaCheckError();

    }
    else{

        // 1D FFT to Az
        cufftHandleError( cufftExecZ2Z(plan_dim3, wfc, wfc, CUFFT_FORWARD) );
        scalarMult<<<grid,threads>>>(wfc, renorm_factor_z, wfc);
        cudaCheckError();

        if(par.bval("Az_time")){
            EqnNode_gpu* Az_eqn = par.astval("Az");
            int e_num = par.ival("Az_num");
            ast_cmult<<<grid,threads>>>(wfc, wfc, Az_eqn, dx, dy, dz,
                                        time, e_num);
            cudaCheckError();
        }
        else{
            cMult<<<grid,threads>>>(wfc, (cufftDoubleComplex*) Az, wfc);
            cudaCheckError();
        }

        cufftHandleError( cufftExecZ2Z(plan_dim3, wfc, wfc, CUFFT_INVERSE) );
        scalarMult<<<grid,threads>>>(wfc, renorm_factor_z, wfc);
        cudaCheckError();


        // loop to multiply by Ay
        for (int i = 0; i < yDim; i++){
            cufftHandleError(cufftExecZ2Z(plan_dim2,  &wfc[i*size],
                             &wfc[i*size], CUFFT_FORWARD));
        }

        scalarMult<<<grid,threads>>>(wfc, renorm_factor_y, wfc);
        cudaCheckError();
        if(par.bval("Ay_time")){
            EqnNode_gpu* Ay_eqn = par.astval("Ay");
            int e_num = par.ival("Ay_num");
            ast_cmult<<<grid,threads>>>(wfc, wfc, Ay_eqn, dx, dy, dz,
                                        time, e_num);
            cudaCheckError();
        }
        else{
            cMult<<<grid,threads>>>(wfc, (cufftDoubleComplex*) Ay, wfc);
            cudaCheckError();
        }

        for (int i = 0; i < yDim; i++){
            //size = xDim * zDim;
            cufftHandleError(cufftExecZ2Z(plan_dim2, &wfc[i*size],
                             &wfc[i*size], CUFFT_INVERSE));
        }
        scalarMult<<<grid,threads>>>(wfc, renorm_factor_y, wfc);
        cudaCheckError();

        // 1d forward / mult by Ax
        cufftHandleError( cufftExecZ2Z(plan_1d, wfc, wfc, CUFFT_FORWARD) );
        scalarMult<<<grid,threads>>>(wfc, renorm_factor_x, wfc);
        cudaCheckError();
        if(par.bval("Ax_time")){
            EqnNode_gpu* Ax_eqn = par.astval("Ax");
            int e_num = par.ival("Ax_num");
            ast_cmult<<<grid,threads>>>(wfc, wfc, Ax_eqn, dx, dy, dz,
                                        time, e_num);
            cudaCheckError();
        }
        else{
            cMult<<<grid,threads>>>(wfc, (cufftDoubleComplex*) Ax, wfc);
            cudaCheckError();
        }
        cufftHandleError( cufftExecZ2Z(plan_1d, wfc, wfc, CUFFT_INVERSE) );
        scalarMult<<<grid,threads>>>(wfc, renorm_factor_x, wfc);
        cudaCheckError();

    }

}

// 2D
void apply_gauge(Grid &par, double2 *wfc, double2 *Ax, double2 *Ay,
                 double renorm_factor_x, double renorm_factor_y, bool flip,
                 cufftHandle plan_1d, cufftHandle plan_dim2, double dx,
                 double dy, double dz, double time){

    dim3 grid = par.grid;
    dim3 threads = par.threads;

    if (flip){

        // 1d forward / mult by Ax
        cufftHandleError( cufftExecZ2Z(plan_1d, wfc, wfc, CUFFT_FORWARD) );
        scalarMult<<<grid,threads>>>(wfc, renorm_factor_y ,wfc);
        cudaCheckError();
        if(par.bval("Ax_time")){
            EqnNode_gpu* Ax_eqn = par.astval("Ax");
            int e_num = par.ival("Ax_num");
            ast_cmult<<<grid,threads>>>(wfc, wfc, Ax_eqn, dx, dy, dz,
                                        time, e_num);
            cudaCheckError();
        }
        else{
            cMult<<<grid,threads>>>(wfc, (cufftDoubleComplex*) Ax, wfc);
            cudaCheckError();
        }
        cufftHandleError( cufftExecZ2Z(plan_1d, wfc, wfc, CUFFT_INVERSE) );
        scalarMult<<<grid,threads>>>(wfc, renorm_factor_x, wfc);
        cudaCheckError();


        // 1D FFT to wfc_pAy
        cufftHandleError( cufftExecZ2Z(plan_dim2, wfc, wfc, CUFFT_FORWARD) );
        scalarMult<<<grid,threads>>>(wfc, renorm_factor_x, wfc);
        cudaCheckError();
        if(par.bval("Ay_time")){
            EqnNode_gpu* Ay_eqn = par.astval("Ay");
            int e_num = par.ival("Ay_num");
            ast_cmult<<<grid,threads>>>(wfc, wfc, Ay_eqn, dx, dy, dz,
                                        time, e_num);
            cudaCheckError();
        }
        else{
            cMult<<<grid,threads>>>(wfc, (cufftDoubleComplex*) Ay, wfc);
            cudaCheckError();
        }

        cufftHandleError( cufftExecZ2Z(plan_dim2, wfc, wfc, CUFFT_INVERSE) );
        scalarMult<<<grid,threads>>>(wfc, renorm_factor_y, wfc);
        cudaCheckError();
    }
    else{

        // 1D FFT to wfc_pAy
        cufftHandleError( cufftExecZ2Z(plan_dim2, wfc, wfc, CUFFT_FORWARD) );
        scalarMult<<<grid,threads>>>(wfc, renorm_factor_x, wfc);
        cudaCheckError();
        if(par.bval("Ay_time")){
            EqnNode_gpu* Ay_eqn = par.astval("Ay");
            int e_num = par.ival("Ay_num");
            ast_cmult<<<grid,threads>>>(wfc, wfc, Ay_eqn, dx, dy, dz,
                                        time, e_num);
            cudaCheckError();
        }
        else{
            cMult<<<grid,threads>>>(wfc, (cufftDoubleComplex*) Ay, wfc);
            cudaCheckError();
        }

        cufftHandleError( cufftExecZ2Z(plan_dim2, wfc, wfc, CUFFT_INVERSE) );
        scalarMult<<<grid,threads>>>(wfc, renorm_factor_y, wfc);
        cudaCheckError();

        // 1d forward / mult by Ax
        cufftHandleError( cufftExecZ2Z(plan_1d, wfc, wfc, CUFFT_FORWARD) );
        scalarMult<<<grid,threads>>>(wfc, renorm_factor_y ,wfc);
        cudaCheckError();

        if(par.bval("Ax_time")){
            EqnNode_gpu* Ax_eqn = par.astval("Ax");
            int e_num = par.ival("Ax_num");
            ast_cmult<<<grid,threads>>>(wfc, wfc, Ax_eqn, dx, dy, dz,
                                        time, e_num);
            cudaCheckError();
        }
        else{
            cMult<<<grid,threads>>>(wfc, (cufftDoubleComplex*) Ax, wfc);
            cudaCheckError();
        }
        cufftHandleError( cufftExecZ2Z(plan_1d, wfc, wfc, CUFFT_INVERSE) );
        scalarMult<<<grid,threads>>>(wfc, renorm_factor_x, wfc);
        cudaCheckError();

    }

}

void evolve(Grid &par,
            int numSteps){

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
    std::vector<double *> V = par.dsvecval("V");
    std::vector<double *> Phi = par.dsvecval("Phi");
    std::vector<double2 *> gpu1dpAx = par.d2svecval("pAx_gpu");
    std::vector<double2 *> gpu1dpAy;
    std::vector<double2 *> gpu1dpAz;
    std::vector<double *> Phi_gpu = par.dsvecval("Phi_gpu");
    bool write_it = par.bval("write_it");
    bool graph = par.bval("graph");
    int N = par.ival("atoms");
    int printSteps = par.ival("printSteps");
    bool nonlin = par.bval("gpe");
    bool lz = par.bval("corotating");
    bool ramp = par.bval("ramp");
    int energy_calc_steps = par.ival("energy_calc_steps") == 0 ? printSteps : par.ival("energy_calc_steps");
    double energy_calc_threshold = par.dval("energy_calc_threshold");
    int ramp_type = par.ival("ramp_type");
    int xDim = par.ival("xDim");
    int yDim = 1;
    int zDim = 1;
    bool gstate = par.bval("gstate");

    int wfc_num = par.ival("wfc_num");

    std::vector<double> energy(wfc_num);

    bool energy_escape = false;

    std::vector<double2 *> wfc_array = par.d2svecval("wfc_array");
    std::vector<double2 *> gpuWfc_array = par.d2svecval("wfc_gpu_array");
    double2 **device_wfc_array;
    cudaHandleError(cudaMalloc((void **) &device_wfc_array,
                              sizeof(double2*)*wfc_num));
    cudaHandleError(cudaMemcpy(device_wfc_array,
                   gpuWfc_array.data(), sizeof(double2*)*wfc_num,
                   cudaMemcpyHostToDevice));

    std::vector<double2 *> K_gpu = par.d2svecval("K_gpu");
    std::vector<double2 *> V_gpu = par.d2svecval("V_gpu");

    double *gpu_interactions;
    cudaHandleError(cudaMalloc((void **) &gpu_interactions,
                    sizeof(double)*wfc_num*wfc_num));

    if (dimnum > 1){
        dy = par.dval("dy");
        y = par.dsval("y");
        gpu1dpAy = par.d2svecval("pAy_gpu");
        yDim = par.ival("yDim");
    }
    if (dimnum > 2){
        dz = par.dval("dz");
        gpu1dpAz = par.d2svecval("pAz_gpu");
        zDim = par.ival("zDim");
    }

    int gridSize = xDim * yDim * zDim;

    // getting data from Cuda class
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
    int iterations;
    if(gstate){
        iterations = par.ival("g_i");
        Dt = gdt;
        printf("Timestep for groundstate solver set as: %E\n",Dt);
    }
    else{
        iterations = par.ival("e_i");
        Dt = dt;
        printf("Timestep for evolution set as: %E\n",Dt);
    }
    begin = clock();
    double omega_0=omega*omegaX;

    // 2D VORTEX TRACKING

    double mask_2d = par.dval("mask_2d");
    int x0_shift = par.dval("x0_shift");
    int y0_shift = par.dval("y0_shift");
    int charge = par.ival("charge");
    int kill_idx = par.ival("kill_idx");
    std::vector<double2 *> EV_opt = par.d2svecval("EV_opt");
    int kick_it = par.ival("kick_it");
    std::vector<double *> V_opt = par.dsvecval("V_opt");
    // Double buffering and will attempt to thread free and calloc operations to
    // hide time penalty. Or may not bother.
    int num_vortices[2] = {0,0};

    std::vector<double *> edges(wfc_num);

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


    std::shared_ptr<Vtx::VtxList> vortCoords 
        = std::make_shared<Vtx::VtxList>(7);
    //std::vector<std::shared_ptr<Vtx::Vortex> vortCoords;

    //Previous array of vortex coordinates from vortexLocation 1's
    //struct Vtx::Vortex *vortCoordsP = NULL;
    //std::vector<struct Vtx::Vortex> vortCoordsP;
    std::shared_ptr<Vtx::VtxList> vortCoordsP 
        = std::make_shared<Vtx::VtxList>(7);


    LatticeGraph::Lattice lattice; //Vortex lattice graph.
    double* adjMat;

    // Assuming triangular lattice at rotatio


    //std::cout << "numSteps is: " << numSteps << '\n';
    // Iterating through all of the steps in either g or esteps.
    for (int i=iterations; i < numSteps+iterations; ++i){

        if (par.bval("energy_calc") && (i % energy_calc_steps == 0)) {
            energy_escape = true;
        }

        for (int w = 0; w < wfc_array.size(); ++w){
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
                } else{
                    if (i == 0){
                        omega_0=(double)omega/(double)(numSteps);
                    }
                    else{
                        omega_0 = (double)(i+1) / (double)i;
                    }
                }
            }
    
            cudaHandleError(cudaMemcpy(wfc_array[w], gpuWfc_array[w],
                            sizeof(cufftDoubleComplex)*xDim*yDim*zDim, 
                            cudaMemcpyDeviceToHost));

            // Print-out at pre-determined rate.
            // Vortex & wfc analysis performed here also.
            if(i % printSteps == 0) {
                // If the unit_test flag is on, we need a special case
                printf("Step: %d    Omega: %lf\n", i, omega_0);

                // Printing out time of iteration
                end = clock();
                time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
                printf("Time spent: %lf\n", time_spent);
                if (ramp == 0 && !gstate){ // Real-time evolution, constant Omega value.
                    if (dimnum == 3){
                        // Note: In the case of 3d, we need to think about
                        //       vortex tracking in a new way.
                        //       It may be as simple as splitting the
                        //       problem into 2D elements and working from 
                        //       there, but let's look into it when we need
                        //       it in the future.
                        std::cout << "commencing 3d vortex tracking\n";

                        // Creating the necessary double* values
                        edges[w] = (double*)malloc(sizeof(double)*gridSize);

                        find_edges(par, wfc_array[w], edges);
                        // Output happens all at once, after the wfc loop
                    } else if (dimnum == 2 && mask_2d > 0){
                        vortexLocation = (int *) calloc(xDim*yDim,
                                                        sizeof(int));
                        num_vortices[w] = Tracker::findVortex(
                            vortexLocation,wfc_array[w],mask_2d,xDim,x,i);
                        // If initial step, locate vortices, least-squares
                        // to find exact centre, calculate lattice angle,
                        // generate optical lattice.
                        if (i == 0) {
                            if(num_vortices[w] > 0){
                                // Reserve enough space for the vortices
                                // reserve(num_vortices[w]);
                                vortCoords = std::make_shared<Vtx::VtxList>
                                                (num_vortices[w]);
                                vortCoordsP = std::make_shared<Vtx::VtxList>
                                                  (num_vortices[w]);

                                // Locate the vortex positions to the nearest 
                                // grid then perform a least-squares fit to 
                                // determine the location to sub-grid reolution
                                Tracker::vortPos(vortexLocation,
                                                vortCoords->getVortices(),
                                                xDim,wfc_array[w]);
                                Tracker::lsFit(vortCoords->getVortices(),
                                                wfc_array[w],
                                                xDim);

                                // Find the centre-most vortex in the lattice
                                central_vortex = Tracker::vortCentre(
                                        vortCoords->getVortices(), xDim);
                                // Determine the Angle formed by the lattice
                                // relative to the x-axis
                                vort_angle = Tracker::vortAngle(vortCoords->
                                        getVortices(), central_vortex);
        
                                // Store the vortex angle in the parameter file
                                par.store("Vort_angle", vort_angle);

                                // Determine average lattice spacing.
                                double sepAvg = Tracker::vortSepAvg(
                                        vortCoords->getVortices(),
                                        central_vortex);
        
                                par.store("Central_vort_x",
                                    (double)central_vortex->getCoords().x);
                                par.store("Central_vort_y",
                                    (double)central_vortex->getCoords().y);
                                par.store("Central_vort_winding",
                                        (double)central_vortex->getWinding());
                                par.store("Num_vort", (double) vortCoords->
                                        getVortices().size());
        
                                // Setup the optical lattice to match the spacing
                                // and angle+angle_sweep of the vortex lattice.
                                // Amplitude matched by setting laser_power
                                // parameter switch.
                                optLatSetup(central_vortex, V[w], 
                                            vortCoords->getVortices(),
                                            vort_angle + PI * angle_sweep/180.0,
                                            laser_power*HBAR*sqrt(omegaX * omegaY),
                                            V_opt[w], x, y, par);
                            }
                            // If kick_it param is 2, perform a single kick
                            // of the optical lattice for the first timestep
                            // only. This is performed by loading the
                            // EV_opt exp(V + V_opt) array into GPU memory
                            // for the potential.
                            if (kick_it == 2) {
                                printf("Kicked it 1\n");
                                cudaHandleError(cudaMemcpy(V_gpu[w],
                                                EV_opt[w],
                                                sizeof(double2)
                                                    *xDim*yDim*wfc_num,
                                                cudaMemcpyHostToDevice));
                            }
                        } else {
                            // If i!=0 and the number of vortices changes
                            // if num_vortices[1] < num_vortices[w],
                            // Fewer vortices
                            if (num_vortices[w] > 0){
                                Tracker::vortPos(vortexLocation, 
                                        vortCoords->getVortices(), xDim,
                                        wfc_array[w]);
                                Tracker::lsFit(vortCoords->
                                                    getVortices(), 
                                                    wfc_array[w], xDim);
                                Tracker::vortArrange(vortCoords->
                                                              getVortices(),
                                                        vortCoordsP->
                                                            getVortices());
                                if(write_it){
                                    FileIO::writeOutInt(data_dir + "vLoc_",
                                                        vortexLocation,
                                                        xDim * yDim,
                                                        i);
                                }
                            }
                        }
    
                        // Used to also defined for vortex elimination using
                        // graph positions and UID numbers.
                        if (graph && num_vortices[w] > 0) {
                            for (int ii = 0;
                                  ii < vortCoords->getVortices().size();
                                  ++ii) {
                                std::shared_ptr<LatticeGraph::Node>
                                    n(new LatticeGraph::Node(
                                        *vortCoords->
                                            getVortices().at(ii).get()));
                                lattice.addVortex(std::move(n));
                            }
                            unsigned int *uids = (unsigned int *) malloc(
                                    sizeof(unsigned int) *
                                    lattice.getVortices().size());
                            for (size_t a=0;
                                  a < lattice.getVortices().size();
                                  ++a){
                                uids[a] = lattice.getVortexIdx(a)->getUid();
                            }
                            if(i==0) {
                                //Lambda for vortex annihilation/creation.
                                auto killIt=[&](int idx, int winding, 
                                                double delta_x,
                                                double delta_y) {
                                    if (abs(delta_x) > 0 ||
                                        abs(delta_y) > 0){
                                        // Killing initial vortex and then 
                                        // imprinting new one
                                        WFC::phaseWinding(Phi[w], 1,
                                            x, y, dx,dy,
                                            lattice.getVortexUid(idx)->
                                                getData().getCoordsD().x,
                                            lattice.getVortexUid(idx)->
                                                getData().getCoordsD().y,
                                            xDim);
    
                                        cudaHandleError(cudaMemcpy(
                                            Phi_gpu[w], Phi[w], 
                                            sizeof(double)
                                                *xDim*yDim*wfc_num,
                                            cudaMemcpyHostToDevice));
                                        cMultPhi <<<grid, threads>>>(
                                            gpuWfc_array[w],Phi_gpu[w],
                                            gpuWfc_array[w]);
                                        cudaCheckError();
    
                                        // Imprinting new one
                                        int cval = -winding;
                                        WFC::phaseWinding(Phi[w], cval,
                                            x,y, dx,dy,
                                            lattice.getVortexUid(idx)->
                                                getData().getCoordsD().x
                                                + delta_x,
                                            lattice.getVortexUid(idx)->
                                                getData().getCoordsD().y
                                                + delta_y,
                                            xDim);
    
                                        // Sending to device for imprinting
                                        cudaHandleError(cudaMemcpy(
                                            Phi_gpu[w], Phi[w], 
                                            sizeof(double)
                                                *xDim*yDim*wfc_num,
                                            cudaMemcpyHostToDevice));
                                        cMultPhi <<<grid, threads>>>(
                                            gpuWfc_array[w],Phi_gpu[w],
                                            gpuWfc_array[w]);
                                        cudaCheckError();
                                    }
                                    else{
                                        int cval = -(winding-1);
                                        WFC::phaseWinding(Phi[w], cval,
                                            x,y,dx,dy,
                                            lattice.getVortexUid(idx)->
                                                getData().getCoordsD().x,
                                            lattice.getVortexUid(idx)->
                                                getData().getCoordsD().y,
                                            xDim);
                                        cudaHandleError(cudaMemcpy(
                                            Phi_gpu[w], Phi[w], 
                                            sizeof(double)
                                                *xDim*yDim*wfc_num,
                                            cudaMemcpyHostToDevice));
                                        cMultPhi <<<grid, threads>>>(
                                            gpuWfc_array[w],Phi_gpu[w],
                                            gpuWfc_array[w]);
                                        cudaCheckError();
                                    }
                                };
                                if (kill_idx > 0){
                                    killIt(kill_idx,charge,x0_shift,
                                            y0_shift);
                                }
                            }
                            lattice.createEdges(1.5 * 2e-5 / dx);
    
                            // Assumes that vortices only form edges when
                            // is up to 1.5*2e-5. Replace with delaunay
                            // delaunay triangulation determined edges 
                            // for better computational scaling (and sanity)

                            // O(n^2) -> terrible implementation.
                            // It works for now.
                            // Generates the adjacency matrix from the
                            // graph and outputs to a Mathematica 
                            // compatible format.
                            adjMat = (double *)calloc(
                                lattice.getVortices().size()
                                * lattice.getVortices().size(),
                                                        sizeof(double));
                            lattice.genAdjMat(adjMat);
                            if (write_it){
                                FileIO::writeOutAdjMat(data_dir+"graph",
                                    adjMat, uids,
                                    lattice.getVortices().size(),
                                    i);
                            }
    
                            //Free and clear all memory blocks
                            free(adjMat);
                            free(uids);
                            lattice.getVortices().clear();
                            lattice.getEdges().clear();
                        }

                        //Write out the vortex locations
                        if(write_it){
                            FileIO::writeOutVortex(data_dir+"vort_arr",
                                vortCoords->getVortices(),i);
                        }
                        printf("Located %lu vortices\n", 
                                vortCoords->getVortices().size());
    
                        //Free memory block for now.
                        free(vortexLocation);

                        //Current values become previous values.
                        num_vortices[1] = num_vortices[w];
                        vortCoords->getVortices().swap(
                            vortCoordsP->getVortices());
                        vortCoords->getVortices().clear();
    
                        }
                    }
                }
    
            // U_r(dt/2)*wfc
            if(nonlin == 1){
                if (wfc_num > 1){
                    double* gpu_interactions = par.dsval("gpu_interactions");
                    if(par.bval("V_time")){
                        EqnNode_gpu* V_eqn = par.astval("V");
                        int e_num = par.ival("V_num");
                        cMultDensity_ast<<<grid,threads>>>(V_eqn,
                            gpuWfc_array[w],
                            gpuWfc_array[w],
                            dx, dy, dz, time, e_num, 0.5*Dt,
                            gstate,interaction*gDenConst);
                        cudaCheckError();
                    } else{
                        cMultDensity_multicomp<<<grid,threads>>>(V_gpu[w],
                            gpuWfc_array[w],
                            gpuWfc_array[w],
                            device_wfc_array,
                            gpu_interactions,
                            0.5*Dt,gstate,gDenConst, wfc_num, w);
                        cudaCheckError();
                    }
                } else{
                    if(par.bval("V_time")){
                        EqnNode_gpu* V_eqn = par.astval("V");
                        int e_num = par.ival("V_num");
                        cMultDensity_ast<<<grid,threads>>>(V_eqn,
                                                          gpuWfc_array[w],
                            gpuWfc_array[w],
                            dx, dy, dz, time, e_num, 0.5*Dt,
                            gstate,interaction*gDenConst);
                        cudaCheckError();
                    } else{
                        cMultDensity<<<grid,threads>>>(V_gpu[w],gpuWfc_array[w],
                            gpuWfc_array[w],
                            0.5*Dt,gstate,interaction*gDenConst);
                        cudaCheckError();
                    }
                }
            } else {
                if(par.bval("V_time")){ 
                    EqnNode_gpu* V_eqn = par.astval("V");
                    int e_num = par.ival("V_num");
                    ast_op_mult<<<grid,threads>>>(gpuWfc_array[w],
                        gpuWfc_array[w],
                        V_eqn, dx, dy, dz, time, e_num, gstate, Dt);
                        cudaCheckError();
                } else{
                    cMult<<<grid,threads>>>(V_gpu[w],gpuWfc_array[w],
                                            gpuWfc_array[w]);
                        cudaCheckError();
                }
            }

            // U_p(dt)*fft2(wfc)
            cufftHandleError(cufftExecZ2Z(plan_3d,gpuWfc_array[w],
                                  gpuWfc_array[w], CUFFT_FORWARD));

            // Normalise
            scalarMult<<<grid,threads>>>(gpuWfc_array[w],renorm_factor_nd,
                                         gpuWfc_array[w]);
            cudaCheckError();
            if (par.bval("K_time")){
                EqnNode_gpu* k_eqn = par.astval("k");
                int e_num = par.ival("k_num");
                ast_op_mult<<<grid,threads>>>(gpuWfc_array[w],gpuWfc_array[w],
                    k_eqn, dx, dy, dz, time, e_num, gstate, Dt);
                cudaCheckError();
            } else{
                cMult<<<grid,threads>>>(K_gpu[w],gpuWfc_array[w],
                                        gpuWfc_array[w]);
                cudaCheckError();
            }
            cufftHandleError(cufftExecZ2Z(plan_3d,gpuWfc_array[w],
                                  gpuWfc_array[w], CUFFT_INVERSE));
    
            // Normalise
            scalarMult<<<grid,threads>>>(gpuWfc_array[w],renorm_factor_nd,
                                         gpuWfc_array[w]);
            cudaCheckError();
    
            // U_r(dt/2)*wfc
            if(nonlin == 1){
                if (wfc_num > 1){
                    double* gpu_interactions = par.dsval("gpu_interactions");
                    if(par.bval("V_time")){
                        EqnNode_gpu* V_eqn = par.astval("V");
                        int e_num = par.ival("V_num");
                        cMultDensity_ast<<<grid,threads>>>(V_eqn,
                            gpuWfc_array[w], gpuWfc_array[w], dx, dy, dz,
                            time, e_num, 0.5*Dt,
                            gstate,interaction*gDenConst);
                        cudaCheckError();
                    } else{
                        cMultDensity_multicomp<<<grid,threads>>>(V_gpu[w],
                            gpuWfc_array[w], gpuWfc_array[w],
                            device_wfc_array, gpu_interactions, 0.5*Dt, gstate,
                            interaction*gDenConst, wfc_num, w);
                        cudaCheckError();
                    }
                } else{
                    if(par.bval("V_time")){
                        EqnNode_gpu* V_eqn = par.astval("V");
                        int e_num = par.ival("V_num");
                        cMultDensity_ast<<<grid,threads>>>(V_eqn,
                            gpuWfc_array[w], gpuWfc_array[w], dx, dy, dz,
                            time, e_num, 0.5*Dt,
                            gstate,interaction*gDenConst);
                        cudaCheckError();
                    } else{
                        cMultDensity<<<grid,threads>>>(V_gpu[w],gpuWfc_array[w],
                            gpuWfc_array[w], 0.5*Dt, gstate,
                            interaction*gDenConst);
                        cudaCheckError();
                    }
                }
            } else {
                if(par.bval("V_time")){  
                    EqnNode_gpu* V_eqn = par.astval("V");
                    int e_num = par.ival("V_num");
                    ast_op_mult<<<grid,threads>>>(gpuWfc_array[w],
                        gpuWfc_array[w],
                        V_eqn, dx, dy, dz, time, e_num, gstate, Dt);
                        cudaCheckError();
                } else{
                    cMult<<<grid,threads>>>(V_gpu[w],gpuWfc_array[w],
                                            gpuWfc_array[w]);
                        cudaCheckError();
                }
            }
    
            // Angular momentum pAy-pAx (if engaged)  //
            if(lz == true){
                // Multiplying by ramping factor if necessary
                // Note: using scalarPow to do the scaling inside of the exp
                if (ramp){
                    scalarPow<<<grid,threads>>>((double2*) gpu1dpAy[w], 
                                                omega_0,
                                                (double2*) gpu1dpAy[w]);
                    cudaCheckError();
                    if (dimnum > 1){
                        scalarPow<<<grid,threads>>>((double2*) gpu1dpAx[w],
                                                    omega_0,
                                                    (double2*) gpu1dpAx[w]);
                        cudaCheckError();
                    }
                    if (dimnum > 2){
                        scalarPow<<<grid,threads>>>((double2*) gpu1dpAz[w],
                                                    omega_0,
                                                    (double2*) gpu1dpAz[w]);
                        cudaCheckError();
                    }
                }
                int size = xDim*zDim;
                if (dimnum == 3){
                    apply_gauge(par, gpuWfc_array[w], gpu1dpAx[w],
                                gpu1dpAy[w], gpu1dpAz[w], renorm_factor_x,
                                renorm_factor_y, renorm_factor_z,
                                i%2, plan_1d, plan_dim2, plan_dim3,
                                dx, dy, dz, time, yDim, size);
                } else if (dimnum == 2){
                    apply_gauge(par, gpuWfc_array[w], gpu1dpAx[w], gpu1dpAy[w],
                                renorm_factor_x, renorm_factor_y, i%2, plan_1d,
                                plan_other2d, dx, dy, dz, time);
                } else if (dimnum == 1){
                    cufftHandleError(cufftExecZ2Z(plan_1d,gpuWfc_array[w],
                                     gpuWfc_array[w],CUFFT_FORWARD));
                    scalarMult<<<grid,threads>>>(gpuWfc_array[w],
                                                 renorm_factor_x,
                                                 gpuWfc_array[w]);
                    cudaCheckError();
                    if(par.bval("Ax_time")){
                        EqnNode_gpu* Ax_eqn = par.astval("Ax");
                        int e_num = par.ival("Ax_num");
                        ast_cmult<<<grid,threads>>>(gpuWfc_array[w],
                            gpuWfc_array[w], Ax_eqn, dx, dy, dz, time, e_num);
                        cudaCheckError();
                    } else{
                        cMult<<<grid,threads>>>(gpuWfc_array[w],
                            (cufftDoubleComplex*) gpu1dpAx[w], gpuWfc_array[w]);
                        cudaCheckError();
                    }

                    cufftHandleError(cufftExecZ2Z(plan_1d,gpuWfc_array[w],
                                     gpuWfc_array[w], CUFFT_INVERSE));
                    scalarMult<<<grid,threads>>>(gpuWfc_array[w],
                                                 renorm_factor_x,
                                                 gpuWfc_array[w]);
                    cudaCheckError();
                }
            }

            if(gstate){
                parSum(gpuWfc_array[w], par);
            }

            if (par.bval("energy_calc") && (i % energy_calc_steps == 0)) {

                double oldEnergy = energy[w];
                energy[w] = energy_calc(par, gpuWfc_array[w]);

                printf("Energy[t@%d wfc@%d]=%E\n",i, w, energy[w]);

                if (i == iterations 
                    || fabs(oldEnergy - energy[w]) >= energy_calc_threshold * oldEnergy
                    || !gstate) {
                    energy_escape = false;
                }
            }
        }

        // Execute instructions that don't depend on the individual wfc
        if (i % printSteps == 0) {
            par.store(gstate ? "g_i" : "e_i", i);
            if (write_it) {
                FileIO::writeOutParams(par);
                FileIO::writeOutWfc(par, wfc_array, i);
            }

            if (ramp == 0 && !gstate){
                if (dimnum == 3) {
                    if (write_it) {
                        FileIO::writeOutEdges(par, edges, i);
                        for (int w = 0; w < wfc_num; w++) {
                            free(edges[w]);
                        }
                    }
                }
            }
        }

        if (par.bval("energy_calc") && (i % energy_calc_steps == 0)) {
            par.store("energy", energy);
            if (write_it) {
                FileIO::writeOutEnergy(par, energy, i);
            }
        }

        if (energy_escape) {
            printf("Stopping early at step %d\n", i);
            break;
        }
    }
    par.store("energy", energy);
    par.store("wfc_array", wfc_array);
    par.store("wfc_gpu_array", gpuWfc_array);

    cudaHandleError(cudaFree(device_wfc_array));
    cudaCheckError();
}
