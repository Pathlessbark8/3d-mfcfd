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
    
    cout<<"Reading Points\n";
    read_input_point_data();
    cout<<"Setting Initial Conditions\n";
    initial_conditions();
    cout<<"Generating Stencil\n";
    generate_split_stencils();

    // cout<<"Testing"<<endl;
    // for (int j = 0; j < point.zpos_nbhs[21]; j++)
	// //
	// {
    //     int k = point.zpos_conn[21][j];
    //     cout<<21<<" "<<j<<" "<<k<<endl;
    // }
    
    points *point_d;
    unsigned long long point_size = sizeof(point);
    cudaStream_t stream;
    cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);
    //
    cudaMalloc(&point_d, point_size);
    cudaMemcpy(point_d, &point, point_size, cudaMemcpyHostToDevice);
    // for (int j = 0; j < point.zpos_nbhs[21]; j++)
	// //
	// {
    //     int k = point.zpos_conn[21][j];
    //     cout<<k<<endl;
    // }
    cudaDeviceSynchronize();
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
    
    cudaMemcpy(&point, point_d, point_size, cudaMemcpyDeviceToHost);
   
    // fstream fout;
    // fout.open("output_prim.dat", ios::out);
    // for(int i=0;i<max_points;++i){
    //     fout<<point.prim[i][0]<<" "<<point.prim[i][1]<<" "<<point.prim[i][2]<<" "<<point.prim[i][3]<<" "<<point.prim[i][4]<<endl;
    // }
    // fout.close();
    cudaFree(point_d);
    //
    cout << "Done\n";
}