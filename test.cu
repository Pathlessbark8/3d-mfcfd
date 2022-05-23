#include <iostream>
#include <fstream>
#include <random>
#include "split_fluxes_mod.h"
#include <iomanip>
#include "octant_fluxes_mod.h"
#include "wall_flux_dGxneg_mod.h"
#include "wall_flux_dGyneg_mod.h"
#include "wall_flux_dGxpos_mod.h"
#include "wall_flux_dGypos_mod.h"
#include "wall_flux_dGzneg_mod.h"
#include "point_preprocessor_mod.h"
#include "compute_conserved_vector_mod.h"
#include "timestep_delt_mod.h"
#include "generate_connectivity_mod.h"
#include "implicit_aliasing_mod.h"
#include "flux_residual_mod.h"
#include "initial_conditions_mod.h"
#include "interior_flux_dGxneg_mod.h"
#include "interior_flux_dGxpos_mod.h"
#include "interior_flux_dGyneg_mod.h"
#include "interior_flux_dGypos_mod.h"
#include <chrono>
#include <cuda_runtime.h>
#include "fpi_solver_mod.h"
using namespace std;
using namespace std::chrono;
int main()
{

    cout<<setprecision(13)<<scientific;

    read_input_point_data();
    initial_conditions();
    generate_split_stencils();
    //
    points *point_d;
    unsigned long long point_size = sizeof(point);
    cudaStream_t stream;
    cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);
    //
    cudaMalloc(&point_d, point_size);
    cudaMemcpy(point_d, &point, point_size, cudaMemcpyHostToDevice);
    // cudaDeviceSynchronize();
    auto start = high_resolution_clock::now();
    cout << "Starting CUDA excecution\n";
    //
    cout << setprecision(13);
    fpi_solver_cuda(point_d,stream);
    //
    cudaDeviceSynchronize();
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    cout << "Time Taken :" << duration.count() / 1000000.0 << endl;
    //
    cudaMemcpy(&point, point_d, point_size, cudaMemcpyDeviceToHost);
    fstream fout;
    fout.open("output_prim.dat", ios::out);
    for(int i=0;i<max_points;++i){
        fout<<point.prim[0][i]<<" "<<point.prim[1][i]<<" "<<point.prim[2][i]<<" "<<point.prim[3][i]<<" "<<point.prim[4][i]<<endl;
    }
    fout.close();
    cudaFree(point_d);
    //
    cout << "Done\n";
}