// #include <iostream>
// #include <fstream>
// #include <random>
// #include "split_fluxes_mod.h"
// #include <iomanip>
// #include "octant_fluxes_mod.h"
// #include "wall_flux_dGxneg_mod.h"
// #include "wall_flux_dGyneg_mod.h"
// #include "wall_flux_dGxpos_mod.h"
// #include "wall_flux_dGypos_mod.h"
// #include "wall_flux_dGzneg_mod.h"
// #include "point_preprocessor_mod.h"
// #include "compute_conserved_vector_mod.h"
// #include "timestep_delt_mod.h"
// #include "generate_connectivity_mod.h"
// #include "implicit_aliasing_mod.h"
// #include "flux_residual_mod.h"
// #include "initial_conditions_mod.h"
// #include "interior_flux_dGxneg_mod.h"
// #include "interior_flux_dGxpos_mod.h"
// #include "interior_flux_dGyneg_mod.h"
// #include "interior_flux_dGypos_mod.h"
// #include <chrono>
// #include <cuda_runtime.h>
// #include "fpi_solver_mod.h"
// #include <nccl.h>
// using namespace std;
// using namespace std::chrono;

// #define CUDACHECK(cmd) do {                         \
//     cudaError_t err = cmd;                            \
//     if (err != cudaSuccess) {                         \
//       printf("Failed: Cuda error %s:%d '%s'\n",       \
//           __FILE__,__LINE__,cudaGetErrorString(err)); \
//       exit(EXIT_FAILURE);                             \
//     }                                                 \
//   } while(0)
  
  
//   #define NCCLCHECK(cmd) do {                         \
//     ncclResult_t res = cmd;                           \
//     if (res != ncclSuccess) {                         \
//       printf("Failed, NCCL error %s:%d '%s'\n",       \
//           __FILE__,__LINE__,ncclGetErrorString(res)); \
//       exit(EXIT_FAILURE);                             \
//     }                                                 \
//   } while(0)


// int main()
// {

//     cout<<setprecision(13)<<scientific;

//     read_input_point_data();
//     initial_conditions();
//     generate_split_stencils();
//     //
//     fstream fin;
//     fin.open("/home/anil/new_3d_code/3d-mfcfd/inputFiles/partFile.dat",ios::in);
//     fin>>numDevices;
//     int counter;
//     for(int i=0;i<max_points;i++){
//         fin>>point.counter[i]>>partVector[i];
//         numberOfPointsPerDevice[partVector[i]]++;
//     }
//     fin.close();
//     //
//     // cout<<numberOfPointsPerDevice[0]<<" "<<numberOfPointsPerDevice[1]<<" "<<numberOfPointsPerDevice[0]+numberOfPointsPerDevice[1]<<endl;
//     //
//     for(int i=0;i<numDevices;i++){
//         splitPoint[i]=new splitPoints[numberOfPointsPerDevice[i]];
//     }
//     // //
//     int *pointCounter=new int[numDevices];
//     for(int i=0;i<numDevices;i++){
//         pointCounter[i]=0;
//     }
//     for(int i=0;i<max_points;i++){
//         assign(splitPoint[partVector[i]][pointCounter[partVector[i]]],i);
//         pointCounter[partVector[i]]++;
//     }
//     //
//     // comms=new ncclComm_t[numDevices];
//     ncclComm_t comms[2];
//     int devs[2]={0,1};
//     // for(int i=0;i<numDevices;i++){
//     //     devs[i]=i;
//     // }
//     cout<<"HI1\n";
    
//     splitPoints** sendbuff = (splitPoints**)malloc(numDevices * sizeof(splitPoints*));
//     splitPoints** recvbuff = (splitPoints**)malloc(numDevices * sizeof(splitPoints*));
//     cudaStream_t* s = (cudaStream_t*)malloc(sizeof(cudaStream_t)*numDevices);
//     cout<<"HI2\n";
//     cout<<numDevices<<endl;
//     for (int i = 0; i < numDevices; ++i) {
//         CUDACHECK(cudaSetDevice(i));
//         cout<<i<<endl;
//         CUDACHECK(cudaMalloc(sendbuff + i, numberOfPointsPerDevice[i] * sizeof(splitPoints)));
//         CUDACHECK(cudaMalloc(recvbuff + i, numberOfPointsPerDevice[i] * sizeof(splitPoints)));
//         // CUDACHECK(cudaMemset(sendbuff[i], 1, numberOfPointsPerDevice[i] * sizeof(splitPoints)));
//         // CUDACHECK(cudaMemset(recvbuff[i], 0, numberOfPointsPerDevice[i] * sizeof(splitPoints)));
//         // CUDACHECK(cudaStreamCreate(s+i));
//     }
//     cout<<"HI3\n";
//     NCCLCHECK(ncclCommInitAll(comms, numDevices, devs));
//     cout<<"HI4\n";
//     // points *point_d;
//     // unsigned long long point_size = sizeof(point);
//     // cudaStream_t stream;
//     // cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);
//     // //
//     // cudaMalloc(&point_d, point_size);
//     // cudaMemcpy(point_d, &point, point_size, cudaMemcpyHostToDevice);
//     // // cudaDeviceSynchronize();
//     // auto start = high_resolution_clock::now();
//     // cout << "Starting CUDA excecution\n";
//     // //
//     // cout << setprecision(13);
//     // fpi_solver_cuda(point_d,stream);
//     // //
//     // cudaDeviceSynchronize();
//     // auto stop = high_resolution_clock::now();
//     // auto duration = duration_cast<microseconds>(stop - start);
//     // cout << "Time Taken :" << duration.count() / 1000000.0 << endl;
//     // //
//     // cudaMemcpy(&point, point_d, point_size, cudaMemcpyDeviceToHost);
//     // fstream fout;
//     // fout.open("output_prim.dat", ios::out);
//     // for(int i=0;i<max_points;++i){
//     //     fout<<point.prim[0][i]<<" "<<point.prim[1][i]<<" "<<point.prim[2][i]<<" "<<point.prim[3][i]<<" "<<point.prim[4][i]<<endl;
//     // }
//     // fout.close();
//     // cudaFree(point_d);
//     //
//     for (int i=0; i<4; i++){
//         ncclCommDestroy(comms[i]);
//     }
//     cout << "Done\n";
// }

#include <stdio.h>
#include "cuda_runtime.h"
#include "nccl.h"

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


int main(int argc, char* argv[])
{
  ncclComm_t comms[4];


  //managing 4 devices
  int nDev = 2;
  int size = 32*1024*1024;
  int devs[2] = { 0, 1 };


  //allocating and initializing device buffers
  float** sendbuff = (float**)malloc(nDev * sizeof(float*));
  float** recvbuff = (float**)malloc(nDev * sizeof(float*));
  cudaStream_t* s = (cudaStream_t*)malloc(sizeof(cudaStream_t)*nDev);


  for (int i = 0; i < nDev; ++i) {
    CUDACHECK(cudaSetDevice(i));
    CUDACHECK(cudaMalloc(sendbuff + i, size * sizeof(float)));
    CUDACHECK(cudaMalloc(recvbuff + i, size * sizeof(float)));
    // CUDACHECK(cudaMemset(sendbuff[i], 1, size * sizeof(float)));
    // CUDACHECK(cudaMemset(recvbuff[i], 0, size * sizeof(float)));
    CUDACHECK(cudaStreamCreate(s+i));
  }


  //initializing NCCL
  NCCLCHECK(ncclCommInitAll(comms, nDev, devs));


   //calling NCCL communication API. Group API is required when using
   //multiple devices per thread
//   NCCLCHECK(ncclGroupStart());
//   for (int i = 0; i < nDev; ++i)
//     NCCLCHECK(ncclAllReduce((const void*)sendbuff[i], (void*)recvbuff[i], size, ncclFloat, ncclSum,
//         comms[i], s[i]));
//   NCCLCHECK(ncclGroupEnd());


//   //synchronizing on CUDA streams to wait for completion of NCCL operation
//   for (int i = 0; i < nDev; ++i) {
//     CUDACHECK(cudaSetDevice(i));
//     CUDACHECK(cudaStreamSynchronize(s[i]));
//   }


//   //free device buffers
//   for (int i = 0; i < nDev; ++i) {
//     CUDACHECK(cudaSetDevice(i));
//     CUDACHECK(cudaFree(sendbuff[i]));
//     CUDACHECK(cudaFree(recvbuff[i]));
//   }


  //finalizing NCCL
  for(int i = 0; i < nDev; ++i)
      ncclCommDestroy(comms[i]);


  printf("Success \n");
  return 0;
}
