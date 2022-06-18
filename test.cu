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
#include <nccl.h>
#include <mpi.h>

// #include "cuPrintf.cu"
using namespace std;
using namespace std::chrono;

// #define MPICHECK(cmd) do {                          \
//     int e = cmd;                                      \
//     if( e != MPI_SUCCESS ) {                          \
//       printf("Failed: MPI error %s:%d '%d'\n",        \
//           __FILE__,__LINE__, e);   \
//       exit(EXIT_FAILURE);                             \
//     }                                                 \
//   } while(0)

#define CUDACHECK(cmd) do {                         \
    cudaError_t err = cmd;                            \
    if (err != cudaSuccess) {                         \
      printf("Failed: Cuda error %s:%d '%s'\n",       \
          __FILE__,__LINE__,cudaGetErrorString(err)); \
      exit(EXIT_FAILURE);                             \
    }                                                 \
  } while(0)
  
  
  #define NCCLCHECK(cmd) do {                         \
    ncclResult_t res = cmd;                           \
    if (res != ncclSuccess) {                         \
      printf("Failed, NCCL error %s:%d '%s'\n",       \
          __FILE__,__LINE__,ncclGetErrorString(res)); \
      exit(EXIT_FAILURE);                             \
    }                                                 \
  } while(0)


int main()
{

    cout<<setprecision(13)<<scientific;

    read_input_point_data();
    initial_conditions();
    generate_split_stencils();
    //
    fstream fin;
    fin.open("/home/anil/new_3d_code/3d-mfcfd/inputFiles/partFile.dat",ios::in);
    fin>>numDevices;
    int counter;
    for(int i=0;i<max_points;i++){
        fin>>point.counter[i]>>partVector[i];
        numberOfPointsPerDevice[partVector[i]]++;
    }
    fin.close();
    //
    // cout<<numberOfPointsPerDevice[0]<<" "<<numberOfPointsPerDevice[1]<<" "<<numberOfPointsPerDevice[0]+numberOfPointsPerDevice[1]<<endl;
    //
    for(int i=0;i<numDevices;i++){
        splitPoint[i]=new splitPoints[numberOfPointsPerDevice[i]];
    }
    int *pointCounter=new int[numDevices];
    for(int i=0;i<numDevices;i++){
        pointCounter[i]=0;
    }
    for(int i=0;i<max_points;i++){
        assign(splitPoint[partVector[i]][pointCounter[partVector[i]]],i);
        pointCounter[partVector[i]]++;
    }
    //
    //
    // int myRank, nRanks;
    // MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    // MPI_Comm_size(MPI_COMM_WORLD, &nRanks);
    //
    comms=new ncclComm_t[numDevices];
    devs=new int[numDevices];
    for(int i=0;i<numDevices;i++){
        devs[i]=i;
    }
    cout<<"HI1\n";
    //
    splitPoints** splitPoint_d = (splitPoints**)malloc(numDevices * sizeof(splitPoints*));
    cudaStream_t *s = (cudaStream_t*)malloc(sizeof(cudaStream_t)*numDevices);
    cout<<"HI2\n";
    cout<<numDevices<<endl;
    //
    // for(int i=0;i<numDevices;i++){
    //     CUDACHECK(cudaStreamCreate(&s[i]));
    // }
    for (int i = 0; i < numDevices; ++i) {
        CUDACHECK(cudaSetDevice(i));
        CUDACHECK(cudaMalloc(splitPoint_d + i, numberOfPointsPerDevice[i] * sizeof(splitPoints)));
        CUDACHECK(cudaStreamCreate(&s[i]));
        CUDACHECK(cudaMemcpyAsync(splitPoint_d[i], splitPoint[i], numberOfPointsPerDevice[i] * sizeof(splitPoints), cudaMemcpyHostToDevice, s[i]));
    }
    //
    cout<<splitPoint[0][0].x<<" "<<splitPoint[1][1].x<<endl;
    cout<<"HI3\n";
    for (int i = 0; i < numDevices; ++i) {
        CUDACHECK(cudaSetDevice(i));
        CUDACHECK(cudaStreamSynchronize(s[i]));
      }
    NCCLCHECK(ncclCommInitAll(comms, numDevices, devs));
    cout<<"HI4\n";
    // //
    fpi_solver_multi_nccl(splitPoint_d,s);
    // //
    cout<<"HI5\n";
    //
    // for (int i = 0; i < numDevices; ++i) {
    //     CUDACHECK(cudaSetDevice(i));
    //     CUDACHECK(cudaStreamSynchronize(s[i]));
    //   }

    for (int i = 0; i < numDevices; ++i) {
        CUDACHECK(cudaSetDevice(i));
        CUDACHECK(cudaMemcpyAsync(splitPoint[i], splitPoint_d[i], numberOfPointsPerDevice[i] * sizeof(splitPoints), cudaMemcpyDeviceToHost, s[i]));
    }
    cout<<"HI6\n";
    // cout<<splitPoint[0][102].counter<<" "<<splitPoint[0][102].q[0]<<" "<<splitPoint[0][102].q[1]<<endl;
    cout<<splitPoint[0][0].x<<endl;
    cout<<"HI7\n";
    //
    // cudaPrintfDisplay(stdout, true);
    // cudaPrintfEnd();
    //
    for (int i=0; i<numDevices; ++i){
        ncclCommDestroy(comms[i]);
    }
    cout << "Done\n";
}

