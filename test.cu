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

#define MPICHECK(cmd) do {                          \
    int e = cmd;                                      \
    if( e != MPI_SUCCESS ) {                          \
      printf("Failed: MPI error %s:%d '%d'\n",        \
          __FILE__,__LINE__, e);   \
      exit(EXIT_FAILURE);                             \
    }                                                 \
  } while(0)

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
    int myRank, nRanks, localRank = 0;
  
  
    //initializing MPI
    MPICHECK(MPI_Init(&argc, &argv));
    MPICHECK(MPI_Comm_rank(MPI_COMM_WORLD, &myRank));
    MPICHECK(MPI_Comm_size(MPI_COMM_WORLD, &nRanks));
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
    for (int i = 0; i < numDevices; ++i) {
        CUDACHECK(cudaSetDevice(i));
        CUDACHECK(cudaMalloc(splitPoint_d + i, numberOfPointsPerDevice[i] * sizeof(splitPoints)));
        CUDACHECK(cudaStreamCreate(&s[i]));
        CUDACHECK(cudaMemcpyAsync(splitPoint_d[i], splitPoint[i], numberOfPointsPerDevice[i] * sizeof(splitPoints), cudaMemcpyHostToDevice, s[i]));
    }
    //
    cout<<splitPoint[0][0].x<<" "<<splitPoint[1][1].x<<endl;
    cout<<"HI3\n";
    if (myRank == 0) ncclGetUniqueId(&id);
    MPICHECK(MPI_Bcast((void *)&id, sizeof(id), MPI_BYTE, 0, MPI_COMM_WORLD));
  
  
    // initializing NCCL, group API is required around ncclCommInitRank as it is
    // called across multiple GPUs in each thread/process
    NCCLCHECK(ncclGroupStart());
    for (int i=0; i<numDevices; i++) {
       CUDACHECK(cudaSetDevice(localRank*numDevices + i));
       NCCLCHECK(ncclCommInitRank(comms+i, nRanks*numDevices, id, myRank*numDevices + i));
    }
    NCCLCHECK(ncclGroupEnd());

    for (int i=0; i<numDevices; i++)
      CUDACHECK(cudaStreamSynchronize(s[i]));
    cout<<"HI4\n";
    // //
    fpi_solver_multi_nccl(splitPoint_d,s);
    // //
    cout<<"HI5\n";
    //
    for (int i = 0; i < numDevices; ++i) {
        CUDACHECK(cudaSetDevice(i));
        CUDACHECK(cudaMemcpyAsync(splitPoint[i], splitPoint_d[i], numberOfPointsPerDevice[i] * sizeof(splitPoints), cudaMemcpyDeviceToHost, s[i]));
    }
    cout<<"HI6\n";
    cout<<splitPoint[0][0].x<<endl;
    cout<<"HI7\n";
    //
    for (int i=0; i<numDevices; i++) {
      CUDACHECK(cudaFree(splitPoint_d[i]));
    }
    //
    for (int i=0; i<numDevices; ++i){
        ncclCommDestroy(comms[i]);
    }
    cout << "Done with process "<<myRank<<endl;
}

// #include <stdio.h>
// #include "cuda_runtime.h"
// #include "nccl.h"
// #include "mpi.h"
// #include <unistd.h>
// #include <stdint.h>
// #include <iostream>


// #define MPICHECK(cmd) do {                          \
//   int e = cmd;                                      \
//   if( e != MPI_SUCCESS ) {                          \
//     printf("Failed: MPI error %s:%d '%d'\n",        \
//         __FILE__,__LINE__, e);   \
//     exit(EXIT_FAILURE);                             \
//   }                                                 \
// } while(0)


// #define CUDACHECK(cmd) do {                         \
//   cudaError_t e = cmd;                              \
//   if( e != cudaSuccess ) {                          \
//     printf("Failed: Cuda error %s:%d '%s'\n",             \
//         __FILE__,__LINE__,cudaGetErrorString(e));   \
//     exit(EXIT_FAILURE);                             \
//   }                                                 \
// } while(0)


// #define NCCLCHECK(cmd) do {                         \
//   ncclResult_t r = cmd;                             \
//   if (r!= ncclSuccess) {                            \
//     printf("Failed, NCCL error %s:%d '%s'\n",             \
//         __FILE__,__LINE__,ncclGetErrorString(r));   \
//     exit(EXIT_FAILURE);                             \
//   }                                                 \
// } while(0)


// static uint64_t getHostHash(const char* string) {
//   // Based on DJB2a, result = result * 33 ^ char
//   uint64_t result = 5381;
//   for (int c = 0; string[c] != '\0'; c++){
//     result = ((result << 5) + result) ^ string[c];
//   }
//   return result;
// }


// static void getHostName(char* hostname, int maxlen) {
//   gethostname(hostname, maxlen);
//   for (int i=0; i< maxlen; i++) {
//     if (hostname[i] == '.') {
//         hostname[i] = '\0';
//         return;
//     }
//   }
// }


// int main(int argc, char* argv[])
// {
//   int size = 32*1024*1024;


//   int myRank, nRanks, localRank = 0;


//   //initializing MPI
//   MPICHECK(MPI_Init(&argc, &argv));
//   MPICHECK(MPI_Comm_rank(MPI_COMM_WORLD, &myRank));
//   MPICHECK(MPI_Comm_size(MPI_COMM_WORLD, &nRanks));

//   printf("Rank %d out of %d processors\n", myRank, nRanks);

//   //calculating localRank which is used in selecting a GPU
//   // uint64_t hostHashs[nRanks];
//   // char hostname[1024];
//   // getHostName(hostname, 1024);
//   // hostHashs[myRank] = getHostHash(hostname);
//   // std::cout<<"HostHash: "<<hostHashs[myRank]<<std::endl;
//   // MPICHECK(MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, hostHashs, sizeof(uint64_t), MPI_BYTE, MPI_COMM_WORLD));
//   // for (int i=0; i<nRanks; i++) {
//   //   // std::cout<<hostHashs[i]<<std::endl;
//   // }
//   // for (int p=0; p<nRanks; p++) {
//   //   // std::cout<<"p: "<<p<<std::endl;
//   //    if (p == myRank) break;
//   //    std::cout<<"p: "<<p<<" HostHash: "<<hostHashs[p]<<" Rank "<<myRank<<std::endl;
//   //    if (hostHashs[p] == hostHashs[myRank]) {
//   //     localRank++;
//   //     std::cout<<localRank<<std::endl;
//   //    }
//   // }


//   //each process is using two GPUs
//   int nDev = 2;


//   float** sendbuff = (float**)malloc(nDev * sizeof(float*));
//   float** recvbuff = (float**)malloc(nDev * sizeof(float*));
//   cudaStream_t* s = (cudaStream_t*)malloc(sizeof(cudaStream_t)*nDev);


//   //picking GPUs based on localRank
//   for (int i = 0; i < nDev; ++i) {
//     // std::cout<<myRank<<" "<<localRank<<" "<<i<<std::endl;
//     CUDACHECK(cudaSetDevice(localRank*nDev + i));
//     CUDACHECK(cudaMalloc(sendbuff + i, size * sizeof(float)));
//     CUDACHECK(cudaMalloc(recvbuff + i, size * sizeof(float)));
//     CUDACHECK(cudaMemset(sendbuff[i], 1, size * sizeof(float)));
//     CUDACHECK(cudaMemset(recvbuff[i], 0, size * sizeof(float)));
//     CUDACHECK(cudaStreamCreate(s+i));
//   }


//   ncclUniqueId id;
//   ncclComm_t comms[nDev];

//   //generating NCCL unique ID at one process and broadcasting it to all
//   if (myRank == 0) ncclGetUniqueId(&id);
//   MPICHECK(MPI_Bcast((void *)&id, sizeof(id), MPI_BYTE, 0, MPI_COMM_WORLD));


//   //initializing NCCL, group API is required around ncclCommInitRank as it is
//   //called across multiple GPUs in each thread/process
//   NCCLCHECK(ncclGroupStart());
//   for (int i=0; i<nDev; i++) {
//      CUDACHECK(cudaSetDevice(localRank*nDev + i));
//      NCCLCHECK(ncclCommInitRank(comms+i, nRanks*nDev, id, myRank*nDev + i));
//   }
//   NCCLCHECK(ncclGroupEnd());


//   //calling NCCL communication API. Group API is required when using
//   //multiple devices per thread/process
//   NCCLCHECK(ncclGroupStart());
//   for (int i=0; i<nDev; i++)
//      NCCLCHECK(ncclAllReduce((const void*)sendbuff[i], (void*)recvbuff[i], size, ncclFloat, ncclSum,
//            comms[i], s[i]));
//   NCCLCHECK(ncclGroupEnd());


//   //synchronizing on CUDA stream to complete NCCL communication
//   for (int i=0; i<nDev; i++)
//       CUDACHECK(cudaStreamSynchronize(s[i]));


//   //freeing device memory
//   for (int i=0; i<nDev; i++) {
//      CUDACHECK(cudaFree(sendbuff[i]));
//      CUDACHECK(cudaFree(recvbuff[i]));
//   }


//   //finalizing NCCL
//   for (int i=0; i<nDev; i++) {
//      ncclCommDestroy(comms[i]);
//   }


//   //finalizing MPI
//   MPICHECK(MPI_Finalize());


//   printf("[MPI Rank %d] Success \n", myRank);
//   return 0;
// }