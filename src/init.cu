#include "../include/init.h"
#include "../include/dynamic.h"

void check_memory(Grid &par){
    int xDim = par.ival("xDim");
    int yDim = par.ival("yDim");
    int zDim = par.ival("zDim");

    bool energy_calc = par.bval("energy_calc");

    int gSize = xDim*yDim*zDim;
    size_t free = 0;
    size_t total = 0;

    cudaMemGetInfo(&free, &total);

    // Note that this check is specifically for the case where we need to keep
    // 8 double2* values on the GPU. This is not the case for dynamic fields
    // and the test should be updated accordingly as these are used more.
    size_t req_memory = 16*8*(size_t)gSize;
    if (energy_calc){
        req_memory += 4*16*(size_t)gSize;
    }
    if (free < req_memory){
        std::cout << "Not enough GPU memory for gridsize!\n";
        std::cout << "Free memory is: " << free << '\n';
        std::cout << "Required memory is: " << req_memory << '\n';
        if (energy_calc){
            std::cout << "Required memory for energy calc is: "
                      << 4*16*(size_t)gSize << '\n';
        }
        std::cout << "xDim is: " << xDim << '\n';
        std::cout << "yDim is: " << yDim << '\n';
        std::cout << "zDim is: " << zDim << '\n';
        std::cout << "gSize is: " << gSize << '\n';
        exit(1);
    }
}

int init(Grid &par){

    check_memory(par);
    set_fns(par);

    // Re-establishing variables from parsed Grid class
    // Initializes uninitialized variables to 0 values
    std::string data_dir = par.sval("data_dir");
    int dimnum = par.ival("dimnum");
    int N = par.ival("atoms");
    int xDim = par.ival("xDim");
    int yDim = par.ival("yDim");
    int zDim = par.ival("zDim");
    bool write_file = par.bval("write_file");
    bool cyl_coord = par.bval("cyl_coord");
    bool corotating = par.bval("corotating");
    dim3 threads;
    unsigned int gSize = xDim;
    if (dimnum > 1){
        gSize *= yDim;
    }
    if (dimnum > 2){
        gSize *= zDim;
    }
    double omega = par.dval("omega");
    double gdt = par.dval("gdt");
    double dt = par.dval("dt");
    double omegaX = par.dval("omegaX");
    double omegaY = par.dval("omegaY");
    double omegaZ = par.dval("omegaZ");
    double gammaY = par.dval("gammaY"); //Aspect ratio of trapping geometry.
    double winding = par.dval("winding");
    double box_size = par.dval("box_size");
    double *Energy;
    double *r;
    double *V_opt;
    double *Energy_gpu;
    cufftDoubleComplex *wfc;
    if (par.bval("read_wfc") == true){
        wfc = par.cufftDoubleComplexval("wfc");
    }
    cufftDoubleComplex *EV_opt;
    cufftDoubleComplex *wfc_backup;
    cufftDoubleComplex *EappliedField;

    std::cout << "gSize is: " << gSize << '\n';
    cufftResult result;
    cufftHandle plan_1d;
    cufftHandle plan_2d;
    cufftHandle plan_3d;
    cufftHandle plan_other2d;
    cufftHandle plan_dim2;
    cufftHandle plan_dim3;

    std::string buffer;
    double Rxy; //Condensate scaling factor.
    double a0x, a0y, a0z; //Harmonic oscillator length in x and y directions

    generate_grid(par);

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//

    double mass = 1.4431607e-25; //Rb 87 mass, kg
    par.store("mass",mass);
    double a_s = 4.76e-9;
    par.store("a_s",a_s);

    double sum = 0.0;

    a0x = sqrt(HBAR/(2*mass*omegaX));
    a0y = sqrt(HBAR/(2*mass*omegaY));
    a0z = sqrt(HBAR/(2*mass*omegaZ));
    par.store("a0x",a0x);
    par.store("a0y",a0y);
    par.store("a0z",a0z);

    // Let's go ahead and define the gDensConst here
    // N*4*HBAR*HBAR*PI*(4.67e-9/mass)*sqrt(mass*(omegaZ)/(2*PI*HBAR)
    double gDenConst = N*4*HBAR*HBAR*PI*(a_s/mass);
    if (dimnum == 2){
        gDenConst*= sqrt(mass*(omegaZ)/(2*PI*HBAR));
    }
    par.store("gDenConst", gDenConst);

    Rxy = pow(15,0.2)*pow(N*a_s*sqrt(mass*omegaZ/HBAR),0.2);
    par.store("Rxy",Rxy);

    //std::cout << "Rxy is: " << Rxy << '\n';
    double xMax, yMax, zMax;
    if (box_size > 0){
        xMax = box_size;
        yMax = box_size;
        zMax = box_size;
    }
    else{
        xMax = 6*Rxy*a0x;
        yMax = 6*Rxy*a0y;
        zMax = 6*Rxy*a0z;
    }
    par.store("xMax",xMax);
    par.store("yMax",yMax);
    par.store("zMax",zMax);

    double pxMax, pyMax, pzMax;
    pxMax = (PI/xMax)*(xDim>>1);
    pyMax = (PI/yMax)*(yDim>>1);
    pzMax = (PI/zMax)*(zDim>>1);
    par.store("pyMax",pyMax);
    par.store("pxMax",pxMax);
    par.store("pzMax",pzMax);

    double dx = xMax/(xDim>>1);
    double dy = yMax/(yDim>>1);
    double dz = zMax/(zDim>>1);
    if (dimnum < 3){
        dz = 1;
    }
    if (dimnum < 2){
        dy = 1;
    }
    par.store("dx",dx);
    par.store("dy",dy);
    par.store("dz",dz);

    double dpx, dpy, dpz;
    dpx = PI/(xMax);
    dpy = PI/(yMax);
    dpz = PI/(zMax);
    //std::cout << "yMax is: " << yMax << '\t' << "xMax is: " << xMax << '\n';
    //std::cout << "dpx and dpy are:" << '\n';
    //std::cout << dpx << '\t' << dpy << '\n';
    par.store("dpx",dpx);
    par.store("dpy",dpy);
    par.store("dpz",dpz);


    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//

    /* Initialise wavefunction, momentum, position, angular momentum,
       imaginary and real-time evolution operators . */
    Energy = (double*) malloc(sizeof(double) * gSize);
    r = (double *) malloc(sizeof(double) * gSize);
    V_opt = (double *) malloc(sizeof(double) * gSize);
    EV_opt = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * gSize);
    EappliedField = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) *
                                                         gSize);

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//

/*
    #ifdef __linux
    int cores = omp_get_num_procs();
    par.store("Cores_Total",cores);

    // Assuming dev system specifics (Xeon with HT -> cores detected / 2)
    par.store("Cores_Max",cores/2);
    omp_set_num_threads(cores/2);

    //#pragma omp parallel for private(j)
    #endif
*/

    par.store("gSize", xDim*yDim*zDim);
    if (par.bval("use_param_file")){
        parse_param_file(par);
    }
    generate_fields(par);
    double *K = par.dsval("K");
    double *Ax = par.dsval("Ax");
    double *Ay = par.dsval("Ay");
    double *Az = par.dsval("Az");
    double *V = par.dsval("V");

    double *pAx = par.dsval("pAx");
    double *pAy = par.dsval("pAy");
    double *pAz = par.dsval("pAz");

    double *x = par.dsval("x");
    double *y = par.dsval("y");
    double *z = par.dsval("z");

    double2 *GpAx = par.cufftDoubleComplexval("GpAx");
    double2 *GpAy = par.cufftDoubleComplexval("GpAy");
    double2 *GpAz = par.cufftDoubleComplexval("GpAz");
    double2 *EpAx = par.cufftDoubleComplexval("EpAx");
    double2 *EpAy = par.cufftDoubleComplexval("EpAy");
    double2 *EpAz = par.cufftDoubleComplexval("EpAz");

    double2 *GV = par.cufftDoubleComplexval("GV");
    double2 *EV = par.cufftDoubleComplexval("EV");
    double2 *GK = par.cufftDoubleComplexval("GK");
    double2 *EK = par.cufftDoubleComplexval("EK");

    wfc = par.cufftDoubleComplexval("wfc");

    int index = 0;
    for(int i=0; i < gSize; i++ ){
        sum+=sqrt(wfc[i].x*wfc[i].x + wfc[i].y*wfc[i].y);
    }

    if (write_file){
        double *Bz;
        double *Bx;
        double *By;
        if (dimnum == 2){
            Bz = curl2d(par, Ax, Ay);
        }
        if (dimnum == 3){
            std::cout << "Calculating the 3d curl..." << '\n';
                    Bx = curl3d_x(par, Ax, Ay, Az);
                    By = curl3d_y(par, Ax, Ay, Az);
                    Bz = curl3d_z(par, Ax, Ay, Az);
                    std::cout << "Finished calculating Curl" << '\n';
        }
        std::cout << "writing initial variables to file..." << '\n';
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
        //hdfWriteDouble(xDim, V, 0, "V_0"); //HDF COMING SOON!
        //hdfWriteComplex(xDim, wfc, 0, "wfc_0");
        if (cyl_coord && dimnum > 2){
            double *Br = curl3d_r(par, Bx, By);
            double *Bphi = curl3d_phi(par, Bx, By);

            FileIO::writeOutDouble(buffer, data_dir + "Br",Br,gSize,0);
            FileIO::writeOutDouble(buffer, data_dir + "Bphi",Bphi,gSize,0);
            FileIO::writeOutDouble(buffer, data_dir + "Bz",Bz,gSize,0);

            free(Br);
            free(Bx);
            free(By);
            free(Bz);
            free(Bphi);
        }
        else{
            if (dimnum > 1){
                FileIO::writeOutDouble(buffer, data_dir + "Bz",Bz,gSize,0);
                free(Bz);
            }
            if (dimnum > 2){
                FileIO::writeOutDouble(buffer, data_dir + "Bx",Bx,gSize,0);
                FileIO::writeOutDouble(buffer, data_dir + "By",By,gSize,0);
                free(Bx);
                free(By);
            }
        }

        FileIO::writeOutDouble(buffer, data_dir + "V",V,gSize,0);
        FileIO::writeOutDouble(buffer, data_dir + "K",K,gSize,0);
        FileIO::writeOutDouble(buffer, data_dir + "pAy",pAy,gSize,0);
        FileIO::writeOutDouble(buffer, data_dir + "pAx",pAx,gSize,0);
        FileIO::writeOutDouble(buffer, data_dir + "Ax",Ax,gSize,0);
        FileIO::writeOutDouble(buffer, data_dir + "Ay",Ay,gSize,0);
        FileIO::writeOutDouble(buffer, data_dir + "Az",Az,gSize,0);
        FileIO::writeOutDouble(buffer, data_dir + "x",x,xDim,0);
        FileIO::writeOutDouble(buffer, data_dir + "y",y,yDim,0);
        FileIO::writeOutDouble(buffer, data_dir + "z",z,zDim,0);
        FileIO::writeOut(buffer, data_dir + "WFC",wfc,gSize,0);
        FileIO::writeOut(buffer, data_dir + "EpAz",EpAz,gSize,0);
        FileIO::writeOut(buffer, data_dir + "EpAy",EpAy,gSize,0);
        FileIO::writeOut(buffer, data_dir + "EpAx",EpAx,gSize,0);
        FileIO::writeOut(buffer, data_dir + "GK",GK,gSize,0);
        FileIO::writeOut(buffer, data_dir + "GV",GV,gSize,0);
        FileIO::writeOut(buffer, data_dir + "GpAx",GpAx,gSize,0);
        FileIO::writeOut(buffer, data_dir + "GpAy",GpAy,gSize,0);
        FileIO::writeOut(buffer, data_dir + "GpAz",GpAz,gSize,0);
    }

    if (par.bval("read_wfc") == false){
        sum=sqrt(sum*dx*dy*dz);
        for (int i = 0; i < gSize; i++){
            wfc[i].x = (wfc[i].x)/(sum);
            wfc[i].y = (wfc[i].y)/(sum);
        }
    }

    result = cufftPlan2d(&plan_2d, xDim, yDim, CUFFT_Z2Z);
    if(result != CUFFT_SUCCESS){
        printf("Result:=%d\n",result);
        printf("Error: Could not execute cufftPlan2d(%s, %d, %d).\n", "plan_2d",
                (unsigned int)xDim, (unsigned int)yDim);
        exit(1);
    }

    generate_plan_other3d(&plan_1d, par, 0);
    if (dimnum == 2){
        generate_plan_other2d(&plan_other2d, par);
    }
    if (dimnum == 3){
        generate_plan_other3d(&plan_dim3, par, 2);
        generate_plan_other3d(&plan_dim2, par, 1);
    }
    result = cufftPlan3d(&plan_3d, xDim, yDim, zDim, CUFFT_Z2Z);
    if(result != CUFFT_SUCCESS){
        printf("Result:=%d\n",result);
        printf("Error: Could not execute cufftPlan3d(%s, %d, %d, %d).\n", 
                "plan_3d",
                (unsigned int)xDim, (unsigned int)yDim, (unsigned int) zDim);
        exit(1);
    }

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//

    //std::cout << GV[0].x << '\t' << GK[0].x << '\t'
    //          << pAy[0] << '\t' << pAx[0] << '\n';

    //std::cout << "storing variables..." << '\n';

    // Storing variables that have been initialized
    // Re-establishing variables from parsed Grid class
    // Initializes uninitialized variables to 0 values
    par.store("Energy", Energy);
    par.store("r", r);
    par.store("Energy_gpu", Energy_gpu);
    par.store("wfc", wfc);
    par.store("EV_opt", EV_opt);
    par.store("V_opt", V_opt);
    par.store("wfc_backup", wfc_backup);
    par.store("EappliedField", EappliedField);

    par.store("result", result);
    par.store("plan_1d", plan_1d);
    par.store("plan_2d", plan_2d);
    par.store("plan_other2d", plan_other2d);
    par.store("plan_3d", plan_3d);
    par.store("plan_dim2", plan_dim2);
    par.store("plan_dim3", plan_dim3);

    // Parameters for time-depd variables.
    par.store("K_time", false);
    par.store("V_time", false);
    par.store("Ax_time", false);
    par.store("Ay_time", false);
    par.store("Az_time", false);


    std::cout << "variables stored" << '\n';

    return 0;
}

void set_variables(Grid &par, bool ev_type){
    // Re-establishing variables from parsed Grid class
    // Note that 3d variables are set to nullptr's unless needed
    //      This might need to be fixed later
    double dx = par.dval("dx");
    double dy = par.dval("dy");
    double *V_opt = par.dsval("V_opt");
    double *pAy = par.dsval("pAy");
    double *pAx = par.dsval("pAx");
    double2 *pAy_gpu;
    double2 *pAx_gpu;
    double2 *pAz_gpu;
    double2 *V_gpu;
    double2 *K_gpu;
    cufftDoubleComplex *wfc = par.cufftDoubleComplexval("wfc");
    cufftDoubleComplex *wfc_gpu = par.cufftDoubleComplexval("wfc_gpu");
    cudaError_t err;
    int dimnum = par.ival("dimnum");
    int xDim = par.ival("xDim");
    int yDim = par.ival("yDim");
    int zDim = par.ival("zDim");
    int gsize = xDim;

    // Special variables for the 3d case
    if (dimnum > 1){
        gsize *= yDim;
    }
    if (dimnum > 2){
        gsize *= zDim;
    }
    if(!par.bval("V_time")){
        cudaMalloc((void**) &V_gpu, sizeof(double2)*gsize);
    }
    if(!par.bval("K_time")){
        cudaMalloc((void**) &K_gpu, sizeof(double2)*gsize);
    }
    if(!par.bval("Ax_time")){
        cudaMalloc((void**) &pAx_gpu, sizeof(double2)*gsize);
    }
    if(!par.bval("Ay_time") && dimnum > 1){
        cudaMalloc((void**) &pAy_gpu, sizeof(double2)*gsize);
    }
    if(!par.bval("Az_time") && dimnum > 2){
        cudaMalloc((void**) &pAz_gpu, sizeof(double2)*gsize);
    }

    if (ev_type == 0){
        cufftDoubleComplex *GK = par.cufftDoubleComplexval("GK");
        cufftDoubleComplex *GV = par.cufftDoubleComplexval("GV");
        cufftDoubleComplex *GpAx = par.cufftDoubleComplexval("GpAx");
        cufftDoubleComplex *GpAy = nullptr;
        cufftDoubleComplex *GpAz = nullptr;

        if(!par.bval("K_time")){
            err=cudaMemcpy(K_gpu, GK, sizeof(cufftDoubleComplex)*gsize,
                           cudaMemcpyHostToDevice);
            if(err!=cudaSuccess){
                std::cout << "ERROR: Could not copy K_gpu to device" << '\n';
                exit(1);
            }
        }
        if(!par.bval("V_time")){
            err=cudaMemcpy(V_gpu, GV, sizeof(cufftDoubleComplex)*gsize,
                           cudaMemcpyHostToDevice);
            if(err!=cudaSuccess){
                std::cout << "ERROR: Could not copy V_gpu to device" << '\n';
                exit(1);
            }
        }
        if(!par.bval("Ax_time")){
            err=cudaMemcpy(pAx_gpu, GpAx, sizeof(cufftDoubleComplex)*gsize,
                           cudaMemcpyHostToDevice);
            if(err!=cudaSuccess){
                std::cout << "ERROR: Could not copy pAx_gpu to device" << '\n';
                exit(1);
            }
        }
        err=cudaMemcpy(wfc_gpu, wfc, sizeof(cufftDoubleComplex)*gsize,
                       cudaMemcpyHostToDevice);
        if(err!=cudaSuccess){
            std::cout << "ERROR: Could not copy wfc_gpu to device" << '\n';
            exit(1);
        }
        par.store("K_gpu", K_gpu);
        par.store("V_gpu", V_gpu);
        par.store("wfc_gpu", wfc_gpu);
        par.store("pAy_gpu", pAy_gpu);
        par.store("pAx_gpu", pAx_gpu);

        // Special cases for 3d
        if (dimnum > 1 && !par.bval("Ay_time")){
            GpAy = par.cufftDoubleComplexval("GpAy");
            err=cudaMemcpy(pAy_gpu, GpAy, sizeof(cufftDoubleComplex)*gsize,
                           cudaMemcpyHostToDevice);

            if(err!=cudaSuccess){
                std::cout << "ERROR: Could not copy pAy_gpu to device" << '\n';
                exit(1);
            }
            par.store("pAy_gpu", pAy_gpu);

        }
        if (dimnum > 2 && !par.bval("Az_time")){
            GpAz = par.cufftDoubleComplexval("GpAz");
            err=cudaMemcpy(pAz_gpu, GpAz, sizeof(cufftDoubleComplex)*gsize,
                           cudaMemcpyHostToDevice);

            if(err!=cudaSuccess){
                std::cout << "ERROR: Could not copy pAz_gpu to device" << '\n';
                exit(1);
            }
            par.store("pAz_gpu", pAz_gpu);

        }
        free(GV); free(GK); free(GpAy); free(GpAx); free(GpAz);
    }
    else if (ev_type == 1){

        cufftDoubleComplex *EV = par.cufftDoubleComplexval("EV");
        cufftDoubleComplex *EK = par.cufftDoubleComplexval("EK");
        cufftDoubleComplex *EpAx = par.cufftDoubleComplexval("EpAx");
        cufftDoubleComplex *EpAy = nullptr;
        cufftDoubleComplex *EpAz = nullptr;
        if (!par.bval("K_time")){
            err=cudaMemcpy(K_gpu, EK, sizeof(cufftDoubleComplex)*gsize,
                           cudaMemcpyHostToDevice);
            if(err!=cudaSuccess){
                std::cout << "ERROR: Could not copy K_gpu to device" << '\n';
                exit(1);
            }
            par.store("K_gpu", K_gpu);
        }
        if(!par.bval("Ax_time")){
            err=cudaMemcpy(pAx_gpu, EpAx, sizeof(cufftDoubleComplex)*gsize,
                           cudaMemcpyHostToDevice);
            if(err!=cudaSuccess){
                std::cout << "ERROR: Could not copy pAx_gpu to device" << '\n';
                std::cout << err << '\n';
                exit(1);
            }
            par.store("pAx_gpu", pAx_gpu);
        }

        if (!par.bval("V_time")){
            err=cudaMemcpy(V_gpu, EV, sizeof(cufftDoubleComplex)*gsize,
                           cudaMemcpyHostToDevice);
            if(err!=cudaSuccess){
                std::cout << "ERROR: Could not copy V_gpu to device" << '\n';
                exit(1);
            }
            par.store("V_gpu", V_gpu);
        }
        err=cudaMemcpy(wfc_gpu, wfc, sizeof(cufftDoubleComplex)*gsize,
                       cudaMemcpyHostToDevice);
        if(err!=cudaSuccess){
            std::cout << "ERROR: Could not copy wfc_gpu to device" << '\n';
            exit(1);
        }

        par.store("wfc_gpu", wfc_gpu);

        // Special variables / instructions for 2/3d case
        if (dimnum > 1 && !par.bval("Ay_time")){
            EpAy = par.cufftDoubleComplexval("EpAy");
            err=cudaMemcpy(pAy_gpu, EpAy, sizeof(cufftDoubleComplex)*gsize,
                           cudaMemcpyHostToDevice);
            if(err!=cudaSuccess){
                std::cout << "ERROR: Could not copy pAy_gpu to device" << '\n';
                exit(1);
            }
            par.store("pAy_gpu", pAy_gpu);
        }

        if (dimnum > 2 && !par.bval("Az_time")){
            EpAz = par.cufftDoubleComplexval("EpAz");
            err=cudaMemcpy(pAz_gpu, EpAz, sizeof(cufftDoubleComplex)*gsize,
                           cudaMemcpyHostToDevice);
            if(err!=cudaSuccess){
                std::cout << "ERROR: Could not copy pAz_gpu to device" << '\n';
                exit(1);
            }
            par.store("pAz_gpu", pAz_gpu);
        }

        free(EV); free(EK); free(EpAy); free(EpAx); free(EpAz);
    }

}

int main(int argc, char **argv){

    Grid par = parseArgs(argc,argv);
    //Grid par2 = parseArgs(argc,argv);

    int device = par.ival("device");
    int dimnum = par.ival("dimnum");
    cudaSetDevice(device);

    std::string buffer;
    time_t start,fin;
    time(&start);
    printf("Start: %s\n", ctime(&start));

    //************************************************************//
    /*
    * Initialise the Params data structure to track params and variables
    */
    //************************************************************//

    // If we want to read in a wfc, we may also need to imprint a phase. This
    // will be done in the init_2d and init_3d functions
    // We need a number of parameters for now
    int xDim = par.ival("xDim");
    int yDim = par.ival("yDim");
    int zDim = par.ival("zDim");
    if(par.bval("read_wfc") == true){

        // Initializing the wfc
        int gSize = xDim * yDim * zDim;
        cufftDoubleComplex *wfc;

        std::string infile = par.sval("infile");
        std::string infilei = par.sval("infilei");
        printf("Loading wavefunction...");
        wfc=FileIO::readIn(infile,infilei,gSize);
        par.store("wfc",wfc);
        printf("Wavefunction loaded.\n");
        //std::string data_dir = par.sval("data_dir");
        //FileIO::writeOut(buffer, data_dir + "WFC_CHECK",wfc,gSize,0);
    }

    init(par);

    int gsteps = par.ival("gsteps");
    int esteps = par.ival("esteps");
    std::string data_dir = par.sval("data_dir");
    std::cout << "variables re-established" << '\n';

    if (par.bval("write_file")){
        FileIO::writeOutParam(buffer, par, data_dir + "Params.dat");
    }

    if(gsteps > 0){
        std::cout << "Imaginary-time evolution started..." << '\n';
        set_variables(par, 0);

        evolve(par, gsteps, 0, buffer);
    }

    if(esteps > 0){
        std::cout << "real-time evolution started..." << '\n';
        set_variables(par, 1);
        evolve(par, esteps, 1, buffer);
    }

    std::cout << "done evolving" << '\n';
    time(&fin);
    printf("Finish: %s\n", ctime(&fin));
    printf("Total time: %ld seconds\n ",(long)fin-start);
    std::cout << '\n';
    return 0;
}

