
#include "../include/ds.h"
#include "../include/operators.h"
#include "../include/split_op.h"

/*----------------------------------------------------------------------------//
* AUX
*-----------------------------------------------------------------------------*/

// I didn't know where to place these functions for now, so the'll be here
void generate_plan_other2d(cufftHandle *plan_fft1d, Grid &par){
    // We need a number of propertied for this fft / transform
    int xDim = par.ival("xDim");
    int yDim = par.ival("yDim");

    int batch = yDim;
    int rank = 1;
    int n[] = {xDim, yDim};
    int idist = 1;
    int odist = 1;
    int inembed[] = {xDim,yDim};
    int onembed[] = {xDim,yDim};
    int istride = yDim;
    int ostride = yDim;

    cufftHandleError( cufftPlanMany(plan_fft1d, rank, n, inembed, istride, 
                                    idist, onembed, ostride, odist, 
                                    CUFFT_Z2Z, batch) );
}

// other plan for 3d case
// Note that we have 3 cases to deal with here: x, y, and z
void generate_plan_other3d(cufftHandle *plan_fft1d, Grid &par, int axis){
    int xDim = par.ival("xDim");
    int yDim = par.ival("yDim");
    int zDim = par.ival("zDim");

    // Along first dimension (x)
    if (axis == 0){
        cufftHandleError( cufftPlan1d(plan_fft1d, xDim, CUFFT_Z2Z, yDim*zDim) );
    }

    // Along second dimension (y)
    // This one is a bit complicated because of how the data is aligned
    else if (axis == 1){
        int batch = yDim;
        int rank = 1;
        int n[] = {xDim, yDim};
        int idist = 1;
        int odist = 1;
        int inembed[] = {xDim, yDim};
        int onembed[] = {xDim, yDim};
        int istride = yDim;
        int ostride = yDim;
    
        cufftHandleError( cufftPlanMany(plan_fft1d, rank, n, inembed, istride, 
                                        idist, onembed, ostride, odist, 
                                        CUFFT_Z2Z, batch) );
        
    }

    // Along third dimension (z)
    else if (axis == 2){

        int batch = xDim*yDim;
        int rank = 1;
        int n[] = {xDim, yDim, zDim};
        int idist = 1;
        int odist = 1;
        int inembed[] = {xDim, yDim, zDim};
        int onembed[] = {xDim, yDim, zDim};
        int istride = xDim*yDim;
        int ostride = xDim*yDim;
    
        cufftHandleError( cufftPlanMany(plan_fft1d, rank, n, inembed, istride, 
                                        idist, onembed, ostride, odist, 
                                        CUFFT_Z2Z, batch) );
    }
}

// Function to set functionPtrs without an unordered map
void set_fns(Grid &par){

    // There are 3 different function distributions to keep in mind:
    // Vfn, Afn, wfcfn

    // Vfn
    par.set_V_fn(par.Vfn);

    // Afn
    par.set_A_fn(par.Afn);

    // Wfcfn
    par.set_wfc_fn(par.Wfcfn);

}


/*----------------------------------------------------------------------------//
* GRID
*-----------------------------------------------------------------------------*/

// Function to store sobel_fft operators and stuff
void Grid::store(std::string id, cufftDoubleComplex *d2param){
    sobel[id] = d2param;
}

// Function to store integer into Grid->param_int
void Grid::store(std::string id, int iparam){
    param_int[id] = iparam;
}

// Function to store double into Grid->param_double
void Grid::store(std::string id, double dparam){
    param_double[id] = dparam;
}

// Function to store double* into param_dstar
void Grid::store(std::string id, double *dsparam){
    param_dstar[id] = dsparam;
}

// Function to store bool into param_bool
void Grid::store(std::string id, bool bparam){
    param_bool[id] = bparam;
}

// Function to store string into data_dir
void Grid::store(std::string id, std::string sparam){
    param_string[id] = sparam;
}

// Function to store asts into data_dir
void Grid::store(std::string id, EqnNode_gpu *ensparam){
    param_ast[id] = ensparam;
}

// Function to store asts into data_dir
void Grid::store(std::string id, EqnNode astparam){
    param_ast_cpu[id] = astparam;
}

// Two boolean functions to check whether a string exists in 
// param_double or param_dstar
bool Grid::is_double(std::string id){
    auto it = param_double.find(id);
    if (it != param_double.end()){
        return true;
    }
    else {
        return false;
    }
}

bool Grid::is_dstar(std::string id){
    auto it = param_dstar.find(id);
    if (it != param_dstar.end()){
        return true;
    }
    else {
        return false;
    }
}

bool Grid::is_ast_gpu(std::string id){
    auto it = param_ast.find(id);
    if (it != param_ast.end()){
        return true;
    }
    else {
        return false;
    }
}

bool Grid::is_ast_cpu(std::string id){
    auto it = param_ast_cpu.find(id);
    if (it != param_ast_cpu.end()){
        return true;
    }
    else {
        return false;
    }
}

// Function to retrieve integer from Grid->param_int
int Grid::ival(std::string id){
    return param_int[id];
}

// Function to retrieve double from Grid->param_double
double Grid::dval(std::string id){
    auto it = param_double.find(id);
    if (it == param_double.end()){
        std::cout << "ERROR: could not find string " << id 
                  << " in Grid::param_double." << '\n';
        assert(it != param_double.end());
    }
    return it->second;
}

// Function to retrieve double star values from param_dstar
double *Grid::dsval(std::string id){
    auto it = param_dstar.find(id);
    if (it == param_dstar.end()){
        std::cout << "ERROR: could not find string " << id
                  << " in Grid::param_dstar." << '\n';
        assert(it != param_dstar.end());
    }
    return it->second;
}

// Function to retrieve bool values from param_bool
bool Grid::bval(std::string id){
    auto it = param_bool.find(id);
    if (it == param_bool.end()){
        std::cout << "ERROR: could not find string " << id 
                  << " in Grid::param_bool." << '\n';
        assert(it != param_bool.end());
    }
    return it->second;
}

// Function to retrieve string from data_dir
std::string Grid::sval(std::string id){
    auto it = param_string.find(id);
    if (it == param_string.end()){
        std::cout << "ERROR: could not find string " << id 
                  << " in Grid::param_string." << '\n';
        assert(it != param_string.end());
    }
    return it->second;
}

// Function to call back the sobel operators
cufftDoubleComplex *Grid::cufftDoubleComplexval(std::string id){
    auto it = sobel.find(id);
    if (it == sobel.end()){
        std::cout << "ERROR: could not find string " << id
                  << " in Grid::sobel." << '\n';
        assert(it != sobel.end());
    }
    return it->second;
}

// Function to call back the ast's
EqnNode_gpu *Grid::astval(std::string id){
    auto it = param_ast.find(id);
    if (it == param_ast.end()){
        std::cout << "ERROR: could not find string " << id
                  << " in Grid::param_ast." << '\n';
        assert(it != param_ast.end());
    }
    return it->second;
}

// Function to call back the ast's
EqnNode Grid::ast_cpuval(std::string id){
    auto it = param_ast_cpu.find(id);
    if (it == param_ast_cpu.end()){
        std::cout << "ERROR: could not find string " << id
                  << " in Grid::param_ast_cpu." << '\n';
        assert(it != param_ast_cpu.end());
    }
    return it->second;
}


// Function for file writing (to replace fileIO::writeOutParam)
void Grid::write(std::string filename){
    std::ofstream output;
    output.open(filename);

    //Needed to recognise Params.dat as .ini format for python post processing
    output << "[Params]" <<"\n";

    // We simply iterate through the int and double param maps
    for (auto item : param_double){
        output << item.first << "=" << item.second << '\n';
        std::cout << item.first << "=" << item.second << '\n';
    }

    for (auto item : param_int){
        output << item.first << "=" << item.second << '\n';
        std::cout << item.first << "=" << item.second << '\n';
    }

    output.close();
}

// Function to print all available variables
void Grid::print_map(){
    for (auto item : param_double){
        std::cout << item.first << '\n';
    }
    for (auto item : param_dstar){
       std::cout << item.first << '\n';
    }
}

void Grid::set_A_fn(std::string id){
    if (id == "rotation"){
        Ax_fn = krotation_Ax;
        Ay_fn = krotation_Ay;
        Az_fn = kconstant_A;
    }
    else if (id == "ring_rotation"){
        Ax_fn = kring_rotation_Ax;
        Ay_fn = kring_rotation_Ay;
        Az_fn = kring_rotation_Az;
    }
    else if (id == "constant"){
        Ax_fn = kconstant_A;
        Ay_fn = kconstant_A;
        Az_fn = kconstant_A;
    }
    else if (id == "ring"){
        Ax_fn = kconstant_A;
        Ay_fn = kconstant_A;
        Az_fn = kring_Az;
    }
    else if (id == "test"){
        Ax_fn = ktest_Ax;
        Ay_fn = ktest_Ay;
        Az_fn = kconstant_A;
    }
    else if (id == "file"){
        Ax_fn = nullptr;
        Ay_fn = nullptr;
        Az_fn = nullptr;
    }
}

void Grid::set_wfc_fn(std::string id){
    if (id == "2d" || id == "3d"){
        wfc_fn = kstd_wfc;
    }
    else if(id == "torus"){
        wfc_fn = ktorus_wfc;
    }
}

void Grid::set_V_fn(std::string id){
    if (id == "2d"){
        V_fn = kharmonic_V;
    }
    else if (id == "torus"){
        V_fn = ktorus_V;
    }
    else if(id == "3d"){
        V_fn = kharmonic_V;
    }
}

