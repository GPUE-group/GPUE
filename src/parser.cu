#include "../include/parser.h"
#include "../include/unit_test.h"
#include "../include/operators.h"

struct stat st = {0};

std::vector<std::string> split_string(std::string input_string,
                                      char split_char){

    // Finding element number to reserve 
    std::vector<int> char_indices = {0};
    for (size_t i = 0; i < input_string.size(); ++i){
        if (input_string[i] == split_char){
            char_indices.push_back(i+1);
        }
    }
     
    std::vector<std::string> result(char_indices.size());
    std::cout << "The result vector size is: " << result.size() << '\n';
    for (size_t i = 0; i < char_indices.size(); ++i){
        result[i] = input_string.substr(char_indices[i],
                                        char_indices[i+1] - char_indices[i]-1);
    }
    for (size_t i = 0; i < result.size(); ++i){
        std::cout << result[i] << '\n';
    }
    return result;
}

void set_interactions(Grid &par, std::string input, int wfc_num,
                      double* interactions, double* gpu_interactions){

    // Setting string to be consistent with rest of function
    if (input[0] == '{'){
        input = input.substr(1, input.size()-2);
    }

    // Stripping unnecessary spaces
    // TODO: does not skip all spaces if multiple spaces are next to each other
    int last_space = -1;
    int i = 0;
    while (i < input.size()){
        if (i == last_space + 1 && input[i] == ' '){
            input.erase(input.begin() + i);
            last_space = i;
            i--;
        }
        if (input[i] == ';' || input[i] == ' '){
            last_space = i;
        }
        i++;
    }

    // splitting along ";"
    std::vector<std::string> rows;
    rows = split_string(input, ';');

    for (size_t i = 0; i < rows.size(); ++i){
        // count when skipping blank columns after splitting
        int count = 0;
        std::vector<std::string> cols = split_string(rows[i], ' ');
        if (cols.size() != wfc_num){
            std::cout << "Improper interaction terms set!\n";
            std::cout << "Number of columns should equal number of wfcs";
            assert(cols.size() == wfc_num);
        }
        for (size_t j = 0; j < cols.size(); ++j){
            if (cols[j].size() > 0){
                int index = wfc_num*i + count;
                interactions[index] = std::stod(cols[j].c_str());
            }
            else{
                count--;
            }
            count++;
        }
    }

    cudaMemcpy(gpu_interactions, interactions, sizeof(double)*wfc_num*wfc_num,
               cudaMemcpyHostToDevice);

    par.store("interactions", interactions);
    par.store("gpu_interactions", gpu_interactions);
}

Grid parseArgs(int argc, char** argv){

    // Creates initial grid class for the parameters
    Grid par;
    int opt;

    // Setting default values
    par.store("xDim", 256);
    par.store("yDim", 256);
    par.store("zDim", 256);
    par.store("omega", 0.0);
    par.store("gammaY", 1.0);
    par.store("gsteps", 1);
    par.store("esteps", 1);
    par.store("gdt", 1e-4);
    par.store("dt", 1e-4);
    par.store("device", 0);
    par.store("atoms", 1);
    par.store("read_wfc", false);
    par.store("printSteps", 100);
    par.store("winding", 0.0);
    par.store("corotating", false);
    par.store("gpe", false);
    par.store("omegaZ", 6.283);
    par.store("interaction",1.0);
    par.store("laser_power",0.0);
    par.store("angle_sweep",0.0);
    par.store("kick_it", 0);
    par.store("write_it", false);
    par.store("x0_shift",0.0);
    par.store("y0_shift",0.0);
    par.store("z0_shift",0.0);
    par.store("sepMinEpsilon",0.0);
    par.store("graph", false);
    par.store("unit_test",false);
    par.store("omegaX", 6.283);
    par.store("omegaY", 6.283);
    par.store("data_dir", (std::string)"data/");
    par.store("ramp", false);
    par.store("ramp_type", 1);
    par.store("dimnum", 2);
    par.store("write_file", true);
    par.store("fudge", 0.0);
    par.store("kill_idx", -1);
    par.store("mask_2d", 0.0);
    par.store("box_size", -0.01);
    par.store("found_sobel", false);
    par.store("energy_calc", false);
    par.store("energy_calc_steps", 0);
    par.store("energy_calc_threshold", -1.0);
    par.store("use_param_file", false);
    par.store("param_file","param.cfg");
    par.store("cyl_coord",false);
    par.store("wfc_num",1);
    par.store("step_offset",0);
    par.Afn = "rotation";
    par.Kfn = "rotation_K";
    par.Vfn = "2d";
    par.Wfcfn = "2d";
    par.store("conv_type", (std::string)"FFT");
    par.store("charge", 0);
    par.store("flip", false);
    par.store("thresh_const", 1.0);

    optind = 1;

    while ((opt = getopt (argc, argv, 
           "b:d:D:C:x:y:w:m:G:g:e:T:t:n:p:rQ:L:E::lsi:P:X:Y:O:k:WU:V:S:ahz:H:uA:v:Z:fc:F:K:R:q:I:j:JN:;")) !=-1)
    {
        switch (opt)
        {
            case 'x':
            {
                int xDim = atoi(optarg);
                printf("Argument for x is given as %d\n",xDim);
                par.store("xDim",(int)xDim);
                break;
            }
            case 'b':
            {
                double box_size = atof(optarg);
                printf("Argument for box_size is given as %E\n",box_size);
                par.store("box_size",(double)box_size);
                break;
            }
            case 'j':
            {
                double thresh_const = atof(optarg);
                printf("Threshold constant is given as %E\n", thresh_const);
                par.store("thresh_const",(double)thresh_const);
                break;
            }
            case 'y':
            {
                int yDim = atoi(optarg);
                printf("Argument for y is given as %d\n",yDim);
                par.store("yDim",(int)yDim);
                break;
            }
            case 'z':
            {
                int zDim = atoi(optarg);
                printf("Argument for z is given as %d\n",zDim);
                par.store("zDim",(int)zDim);
                break;
            }
            case 'w':
            {
                double omega = atof(optarg);
                printf("Argument for OmegaRotate is given as %E\n",omega);
                par.store("omega",(double)omega);
                break;
            }
            case 'm':
            {
                double mask_2d = atof(optarg);
                printf("Argument for mask_2d is given as %E\n",mask_2d);
                par.store("mask_2d",(double)mask_2d);
                break;
            }
            case 'G':
            {
                double gammaY = atof(optarg);
                printf("Argument for gamma is given as %E\n",gammaY);
                par.store("gammaY",(double)gammaY);
                break;
            }
            case 'g':
            {
                double gsteps = atof(optarg);
                printf("Argument for Groundsteps is given as %E\n",gsteps);
                par.store("gsteps",(int)gsteps);
                break;
            }
            case 'e':
            {
                double esteps = atof(optarg);
                printf("Argument for EvSteps is given as %E\n",esteps);
                par.store("esteps",(int)esteps);
                break;
            }
            case 'F':
            {
                double fudge = atof(optarg);
                printf("Argument for Fudge Factor is given as %E\n",fudge);
                par.store("fudge",fudge);
                break;
            }
            case 'T':
            {
                double gdt = atof(optarg);
                printf("Argument for groundstate Timestep is given as %E\n",
                       gdt);
                par.store("gdt",(double)gdt);
                break;
            }
            case 't':
            {
                double dt = atof(optarg);
                printf("Argument for Timestep is given as %E\n",dt);
                par.store("dt",(double)dt);
                break;
            }
            case 'o':
            {
                int step_offset = atoi(optarg);
                printf("Argument for step offset is given as %d\n",step_offset);
                par.store("step_offset",(int)step_offset);
                break;
            }
            case 'C':
            {
                int device = atoi(optarg);
                printf("Argument for device (Card) is given as %d\n",device);
                par.store("device",(int)device);
                break;
            }
            case 'N':
            {
                int wfc_num = atoi(optarg);

                // Initializing the interaction terms
                double *interactions, *gpu_interactions;
                interactions = (double *)malloc(sizeof(double)*wfc_num*wfc_num);
                cudaMalloc((void **) &gpu_interactions,
                           sizeof(double)*wfc_num*wfc_num);

                for (int i = 0; i < wfc_num; ++i){
                    for (int j = 0; j < wfc_num; ++j){
                        int index = wfc_num*i +j;
                        double val;
                        if (i == j){
                            val = 1;
                        }
                        else{
                            val = 0.5;
                        }
                        interactions[index] = val;
                    }
                }

                cudaMemcpy(gpu_interactions, interactions,
                           sizeof(double)*wfc_num*wfc_num,
                           cudaMemcpyHostToDevice);

                printf("Argument for condensate number is given as %d\n",
                        wfc_num);
                par.store("wfc_num",(int)wfc_num);
                par.store("interactions", interactions);
                par.store("gpu_interactions", gpu_interactions);
                break;
            }
            case 'n':
            {
                double atoms = atof(optarg);
                printf("Argument for atoms is given as %E\n",atoms);
                par.store("atoms",(int)atoms);
                break;
            }
            case 'I':
            {
                par.store("param_file", (std::string)optarg);
                par.store("use_param_file", true);
                break;
            }
            case 'R':
            {
                int ramp_type = atoi(optarg);
                printf("Ramping omega with imaginary time evolution\n");
                par.store("ramp",true);
                par.store("ramp_type", ramp_type);
                break;
            }
            case 'r':
            {
                printf("Reading wavefunction from file.\n");
                par.store("read_wfc",true);
                break;
            }
            case 'p':
            {
                int print = atoi(optarg);
                printf("Argument for Printout is given as %d\n",print);
                par.store("printSteps",(int)print);
                break;
            }
            case 'L':
            {
                double l = atof(optarg);
                printf("Vortex winding is given as : %E\n",l);
                par.store("winding",(double)l);
                break;
            }
            case 'l':
            {
                printf("Angular momentum mode engaged\n");
                par.store("corotating",true);
                break;
            }
            case 'E':
            {
                if (optind >= argc || argv[optind][0] == '-') {
                    printf("Energy tag set but no options given!\n");
                } else {
                    double threshold = atof(argv[optind]);
                    int steps = atoi(argv[optind + 1]);

                    printf("Calculating energy every %d steps, stopping if difference ratio is less than %E\n", steps, threshold);
                    par.store("energy_calc",true);
                    par.store("energy_calc_steps", steps);
                    par.store("energy_calc_threshold", threshold);

                    optind++;
                }
                par.store("energy_calc", true);
                break;
            }
            case 'f':
            {
                printf("No longer writing initial variables to file.\n");
                par.store("write_file", false);
                break;
            }
            case 'J':
            {
                printf("Using cylindrical coordinates for B field\n");
                par.store("cyl_coord", true);
                break;
            }
            case 'q':
            {
                int q = atoi(optarg);
                std::cout << "Imprinting vortex with charge q " << q << '\n';
                par.store("flip", true);
                par.store("charge",(int)q);
                break;
            }
            case 's':
            {
                printf("Non-linear mode engaged\n");
                par.store("gpe",true);
                break;
            }
            case 'Z':
            {
                double omegaZ = atof(optarg);
                printf("Argument for OmegaZ is given as %E\n",omegaZ);
                par.store("omegaZ",(double)omegaZ);
                break;
            }
            case 'h':
            {
                std::string command = "src/print_help.sh ";
                system(command.c_str());
                exit(0);
                break;
            }
            case 'H':
            {
                std::string command = "src/print_help.sh ";
                command.append(optarg);
                system(command.c_str());
                exit(0);
                break;
            }
            case 'i':
            {
                int wfc_num = par.ival("wfc_num");
                if (wfc_num == 1){
                    double interaction = atof(optarg);
                    printf("Argument for interaction scaling is %E\n",
                           interaction);
                    par.store("interaction",interaction);
                }
                if (wfc_num > 1){
                    double* interactions = par.dsval("interactions");
                    double* gpu_interactions = par.dsval("gpu_interactions");
                    set_interactions(par, std::string(optarg), wfc_num,
                                     interactions, gpu_interactions);
                    printf("Arguments set for interaction terms\n");
                }
                break;
            }
            case 'P':
            {
                double laser_power = atof(optarg);
                printf("Argument for laser power is %E\n",laser_power);
                par.store("laser_power",laser_power);
                break;
            }
            case 'X':
            {
                double omegaX = atof(optarg);
                printf("Argument for omegaX is %E\n",omegaX);
                par.store("omegaX",(double)omegaX);
                break;
            }
            case 'Y':
            {
                double omegaY = atof(optarg);
                printf("Argument for omegaY is %E\n",omegaY);
                par.store("omegaY",omegaY);
                break;
            }
            case 'O':
            {
                double angle_sweep = atof(optarg);
                printf("Argument for angle_sweep is %E\n",angle_sweep);
                par.store("angle_sweep",angle_sweep);
                break;
            }
            case 'k':
            {
                int kick_it = atoi(optarg);
                printf("Argument for kick_it is %i\n",kick_it);
                par.store("kick_it",kick_it);
                break;
            }
            case 'W':
            {
                printf("Writing out\n");
                par.store("write_it",true);
                break;
            }
            case 'd':
            {
                std::string data_dir = optarg;
                std::cout << "Data directory is: " << data_dir << '\n';
                par.store("data_dir", data_dir + "/");
                break;
            }
            case 'U':
            {
                double x0_shift = atof(optarg);
                printf("Argument for x0_shift is %lf\n",x0_shift);
                par.store("x0_shift",x0_shift);
                break;
            }
            case 'u':
            {
                std::cout << "performing all unit tests" << '\n';
                par.store("unit_test", true);
                test_all();
                exit(0);
            }
            case 'V':
            {
                double y0_shift = atof(optarg);
                printf("Argument for y0_shift is %lf\n",y0_shift);
                par.store("y0_shift",y0_shift);
                break;
            }
            case 'v':
            {
                std::string pot = optarg;
                std::cout << "Chosen potential is: " << pot << '\n';
                par.Vfn = pot;
                par.Wfcfn = pot;
                break;
            }
            case 'S':
            {
                double sepMinEpsilon = atof(optarg);
                printf("Argument for sepMinEpsilon is %lf\n",sepMinEpsilon);
                par.store("sepMinEpsilon",sepMinEpsilon);
                break;
            }
            case 'Q':
            {
                double z0_shift = atof(optarg);
                printf("Argument for z0_shift is %lf\n",z0_shift);
                par.store("z0_shift",z0_shift);
                break;
            }
            case 'c':
            {
                int dimnum = atoi(optarg);
                printf("Argument for number of coordinates is %d\n",dimnum);

                //setting 3d parameters
                if (dimnum == 3){
                    par.Kfn = "rotation_K3d";
                    par.Vfn = "3d";
                    par.Wfcfn = "3d";
                    par.store("box_size", 2.5e-5);
                }
                par.store("dimnum",(int)dimnum);
                break;
            }

            // this case is special and may require reading input from a file
            // or from cin
            case 'A':
            {
                std::string field = optarg;
                std::cout << "Chosen gauge field is: " << field << '\n';
                par.Afn = field;
                break;
            }
            case 'a':
            {
                printf("Graphing mode engaged\n");
                par.store("graph",true);
                break;
            }

            case 'K':
            {
                int kill_idx = atoi(optarg);
                printf("Argument for kill_idx is %d\n",kill_idx);
                par.store("kill_idx",kill_idx);
                break;
            }

            case '?':
            {
                if (optopt == 'c') {
                    fprintf (stderr, 
                             "Option -%c requires an argument.\n", optopt);
                } 
                else if (isprint (optopt)) {
                    fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                } 
                else {
                    fprintf (stderr,
                             "Unknown option character `\\x%x'.\n",optopt);
                }
                return par;
            default:
                abort ();
            }
        }
    }

    std::string data_dir = par.sval("data_dir");
    int dimnum = par.ival("dimnum");

    // Setting variables
    if (stat(data_dir.c_str(), &st) == -1) {
        mkdir(data_dir.c_str(), 0700);
    }

    if (dimnum < 3){
        par.store("zDim", 1);
    }
    if (dimnum < 2){
        par.store("yDim", 1);
    }

    // Update values which depend on other values, so that they don't need to be entered in order

    if (par.bval("use_param_file")) {
        std::string param_file = filecheck(data_dir 
            + par.sval("param_file"));
        std::cout << "Input parameter file is " <<  param_file << '\n';
        par.store("param_file", (std::string)param_file);
    }

    if (par.bval("read_wfc")) {
        std::string infile = filecheck(data_dir + "wfc_load");
        std::string infilei = filecheck(data_dir + "wfci_load");
        par.store("infile", infile);
        par.store("infilei", infilei);
    }

    // If the file gauge field is chosen, we need to make sure the files exist
    if (par.Afn.compare("file") == 0){
        std::cout << "Finding file for Ax..." << '\n';
        par.Axfile = filecheck(data_dir + "Axgauge");
        std::cout << "Finding file for Ay..." << '\n';
        par.Ayfile = filecheck(data_dir + "Aygauge");
        if (dimnum == 3){
            std::cout << "Finding file for Az..." << '\n';
            par.Azfile = filecheck(data_dir + "Azgauge");
        }
    }

    return par;
}
