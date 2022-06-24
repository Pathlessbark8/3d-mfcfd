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
// #include <mpi.h>

// // #include "cuPrintf.cu"
// using namespace std;
// using namespace std::chrono;

// #define MPICHECK(cmd) do {                          \
//     int e = cmd;                                      \
//     if( e != MPI_SUCCESS ) {                          \
//       printf("Failed: MPI error %s:%d '%d'\n",        \
//           __FILE__,__LINE__, e);   \
//       exit(EXIT_FAILURE);                             \
//     }                                                 \
//   } while(0)

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

//   int main(int argc, char* argv[])
//   {
  
//       cout<<setprecision(13)<<scientific;
  
//       read_input_point_data();
//       initial_conditions();
//       generate_split_stencils();
//       //
//       fstream fin;
//       fin.open("/home/anil/new_3d_code/3d-mfcfd/inputFiles/partFile.dat",ios::in);
//       fin>>numDevices;
//       int counter;
//       for(int i=0;i<max_points;i++){
//           fin>>point.counter[i]>>partVector[i];
//           numberOfPointsPerDevice[partVector[i]]++;
//       }
//       fin.close();
//       //
//       for(int i=0;i<numDevices;i++){
//           splitPoint[i]=new splitPoints[numberOfPointsPerDevice[i]];
//       }
//       int *pointCounter=new int[numDevices];
//       for(int i=0;i<numDevices;i++){
//           pointCounter[i]=0;
//       }
//       for(int i=0;i<max_points;i++){
//           assign(splitPoint[partVector[i]][pointCounter[partVector[i]]],i);
//           pointCounter[partVector[i]]++;
//       }
//       //
//       comms=new ncclComm_t[numDevices];
//       devs=new int[numDevices];
//       for(int i=0;i<numDevices;i++){
//           devs[i]=i;
//       }
//       cout<<"HI1\n";
//       int nGPUs;
//       cudaGetDeviceCount(&nGPUs);
//       cout<<nGPUs<<endl;
//       //
//       splitPoints** splitPoint_d = (splitPoints**)malloc(numDevices * sizeof(splitPoints*));
//       cudaStream_t *s = (cudaStream_t*)malloc(sizeof(cudaStream_t)*numDevices);
//       cout<<"HI2\n";
//       cout<<numDevices<<endl;
//       //
//       for (int i = 0; i < numDevices; ++i) {
//           CUDACHECK(cudaSetDevice(i));
//           CUDACHECK(cudaMalloc(splitPoint_d + i, numberOfPointsPerDevice[i] * sizeof(splitPoints)));
//           CUDACHECK(cudaStreamCreate(&s[i]));
//           CUDACHECK(cudaMemcpyAsync(splitPoint_d[i], splitPoint[i], numberOfPointsPerDevice[i] * sizeof(splitPoints), cudaMemcpyHostToDevice, s[i]));
//       }
//       //
//       cout<<splitPoint[0][0].x<<" "<<splitPoint[1][1].x<<endl;
//       cout<<"HI3\n";
//       // if (myRank == 0) ncclGetUniqueId(&id);
//       // MPICHECK(MPI_Bcast((void *)&id, sizeof(id), MPI_BYTE, 0, MPI_COMM_WORLD));
    
    
//       // initializing NCCL, group API is required around ncclCommInitRank as it is
//       // called across multiple GPUs in each thread/process
//       // NCCLCHECK(ncclGroupStart());
//       // for (int i=0; i<numDevices; i++) {
//       //    CUDACHECK(cudaSetDevice(localRank*numDevices + i));
//       //    NCCLCHECK(ncclCommInitRank(comms+i, nRanks*numDevices, id, myRank*numDevices + i));
//       // }
//       // NCCLCHECK(ncclGroupEnd());
  
//       for (int i=0; i<numDevices; i++)
//         CUDACHECK(cudaStreamSynchronize(s[i]));
//       cout<<"HI4\n";
//       // //
//       fpi_solver_multi_nccl(splitPoint_d,s);
//       // //
//       cout<<"HI5\n";
//       //
//       for (int i = 0; i < numDevices; ++i) {
//           CUDACHECK(cudaSetDevice(i));
//           CUDACHECK(cudaMemcpyAsync(splitPoint[i], splitPoint_d[i], numberOfPointsPerDevice[i] * sizeof(splitPoints), cudaMemcpyDeviceToHost, s[i]));
//       }
//       cout<<"HI6\n";
//       cout<<splitPoint[0][10].x<<endl;
//       cout<<"HI7\n";
//       //
//       for (int i=0; i<numDevices; i++) {
//         CUDACHECK(cudaFree(splitPoint_d[i]));
//       }
//       //
//       for (int i=0; i<numDevices; ++i){
//           ncclCommDestroy(comms[i]);
//       }
//       cout << "Done with process "<<endl;
//   }


// // int main(int argc, char* argv[])
// // {

// //     cout<<setprecision(13)<<scientific;

// //     read_input_point_data();
// //     initial_conditions();
// //     generate_split_stencils();
// //     //
// //     fstream fin;
// //     fin.open("/home/anil/new_3d_code/3d-mfcfd/inputFiles/partFile.dat",ios::in);
// //     fin>>numDevices;
// //     int counter;
// //     for(int i=0;i<max_points;i++){
// //         fin>>point.counter[i]>>partVector[i];
// //         numberOfPointsPerDevice[partVector[i]]++;
// //     }
// //     fin.close();
// //     //
// //     // cout<<numberOfPointsPerDevice[0]<<" "<<numberOfPointsPerDevice[1]<<" "<<numberOfPointsPerDevice[0]+numberOfPointsPerDevice[1]<<endl;
// //     //
// //     for(int i=0;i<numDevices;i++){
// //         splitPoint[i]=new splitPoints[numberOfPointsPerDevice[i]];
// //     }
// //     int *pointCounter=new int[numDevices];
// //     for(int i=0;i<numDevices;i++){
// //         pointCounter[i]=0;
// //     }
// //     for(int i=0;i<max_points;i++){
// //         assign(splitPoint[partVector[i]][pointCounter[partVector[i]]],i);
// //         pointCounter[partVector[i]]++;
// //     }
// //     //
// //     //
// //     // int myRank, nRanks, localRank = 0;
  
  
// //     //initializing MPI
// //     // MPICHECK(MPI_Init(&argc, &argv));
// //     // MPICHECK(MPI_Comm_rank(MPI_COMM_WORLD, &myRank));
// //     // MPICHECK(MPI_Comm_size(MPI_COMM_WORLD, &nRanks));
// //     //
// //     comms=new ncclComm_t[numDevices];
// //     devs=new int[numDevices];
// //     for(int i=0;i<numDevices;i++){
// //         devs[i]=i;
// //     }
// //     cout<<"HI1\n";
// //     //
// //     splitPoints** splitPoint_d = (splitPoints**)malloc(numDevices * sizeof(splitPoints*));
// //     cudaStream_t *s = (cudaStream_t*)malloc(sizeof(cudaStream_t)*numDevices);
// //     cout<<"HI2\n";
// //     cout<<numDevices<<endl;
// //     //
// //     for (int i = 0; i < numDevices; ++i) {
// //         CUDACHECK(cudaSetDevice(i));
// //         CUDACHECK(cudaMalloc(splitPoint_d + i, numberOfPointsPerDevice[i] * sizeof(splitPoints)));
// //         CUDACHECK(cudaStreamCreate(&s[i]));
// //         CUDACHECK(cudaMemcpyAsync(splitPoint_d[i], splitPoint[i], numberOfPointsPerDevice[i] * sizeof(splitPoints), cudaMemcpyHostToDevice, s[i]));
// //     }
// //     //
// //     cout<<splitPoint[0][0].x<<" "<<splitPoint[1][1].x<<endl;
// //     cout<<"HI3\n";
// //     // if (myRank == 0) ncclGetUniqueId(&id);
// //     // MPICHECK(MPI_Bcast((void *)&id, sizeof(id), MPI_BYTE, 0, MPI_COMM_WORLD));
  
  
// //     // initializing NCCL, group API is required around ncclCommInitRank as it is
// //     // called across multiple GPUs in each thread/process
// //     // NCCLCHECK(ncclGroupStart());
// //     // for (int i=0; i<numDevices; i++) {
// //     //    CUDACHECK(cudaSetDevice(localRank*numDevices + i));
// //     //    NCCLCHECK(ncclCommInitRank(comms+i, nRanks*numDevices, id, myRank*numDevices + i));
// //     // }
// //     // NCCLCHECK(ncclGroupEnd());

// //     for (int i=0; i<numDevices; i++)
// //       CUDACHECK(cudaStreamSynchronize(s[i]));
// //     cout<<"HI4\n";
// //     // //
// //     fpi_solver_multi_nccl(splitPoint_d,s);
// //     // //
// //     cout<<"HI5\n";
// //     //
// //     for (int i = 0; i < numDevices; ++i) {
// //         CUDACHECK(cudaSetDevice(i));
// //         CUDACHECK(cudaMemcpyAsync(splitPoint[i], splitPoint_d[i], numberOfPointsPerDevice[i] * sizeof(splitPoints), cudaMemcpyDeviceToHost, s[i]));
// //     }
// //     cout<<"HI6\n";
// //     cout<<splitPoint[0][10].x<<endl;
// //     cout<<"HI7\n";
// //     //
// //     for (int i=0; i<numDevices; i++) {
// //       CUDACHECK(cudaFree(splitPoint_d[i]));
// //     }
// //     //
// //     for (int i=0; i<numDevices; ++i){
// //         ncclCommDestroy(comms[i]);
// //     }
// //     cout << "Done with process "<<endl;
// // }

// #include <nccl.h>
// #include <cstdio>
// #include <cstdlib>
// #include <iostream>
// using namespace std;
 
// __global__ void kernel(int *a,int g) 
// {
//   int index = threadIdx.x;

//   a[index] *= (g+1);
//   printf("%d\t", a[index]);

// }/*kernel*/
 

// void print_vector(int *in, int n){

//  for(int i=0; i < n; i++)
//   printf("%d\t", in[i]);

//  printf("\n");

// }/*print_vector*/


// int main(int argc, char* argv[]) {

//   int data_size = 8 ;
//   int nGPUs = 0;
//   cudaGetDeviceCount(&nGPUs);
  
//   int *DeviceList = (int *) malloc (nGPUs     * sizeof(int));
//   int *data       = (int*)  malloc (data_size * sizeof(int));
//   int **d_data    = (int**) malloc (nGPUs     * sizeof(int*));
  
//   for(int i = 0; i < nGPUs; i++)
//       DeviceList[i] = i;
  
//       ncclUniqueId id;
//       ncclGetUniqueId(&id);
//       // printf("ID: %d\n",id);

//   /*Initializing NCCL with Multiples Devices per Thread*/
//   ncclComm_t* comms = (ncclComm_t*)  malloc(sizeof(ncclComm_t)  * nGPUs);  
//   cudaStream_t* s   = (cudaStream_t*)malloc(sizeof(cudaStream_t)* nGPUs);
//   ncclGroupStart();
//   for(int i=0;i<nGPUs;i++){
//     cudaSetDevice(DeviceList[i]);
//     ncclCommInitRank(comms, nGPUs, id,i);
//   }
//   ncclGroupEnd();
  
//   ncclResult_t asyncError;
//   ncclCommGetAsyncError(comms[0], &asyncError);
//   cout<<"Error :"<<asyncError<<endl;

//    int RanksinNCCLComm;
//    ncclCommCount( comms[1], &RanksinNCCLComm);
//    cout<<"RanksinNCCLComm "<<RanksinNCCLComm<<endl;

//    int AssociatedCudaDeviceRank;
//    for(int i=0;i<nGPUs;i++){
//     ncclCommCuDevice( comms[i], &AssociatedCudaDeviceRank);
//     cout<<"AssociatedCudaDeviceRank "<<AssociatedCudaDeviceRank<<endl;
//    }

//    int RankofComm;
//    ncclCommUserRank( comms[1], &RankofComm);
//    cout<<"RankofComm "<<RankofComm<<endl;
//   /*Population the data vector*/
//   // for(int i = 0; i < data_size; i++)
//   //     data[i] = rand()%(10-2)*2;
 
//   // // print_vector(data, data_size);
      
//   // for(int g = 0; g < nGPUs; g++) {
//   //     cudaSetDevice(DeviceList[g]);
//   //     cudaStreamCreate(&s[g]);
//   //     cudaMalloc(&d_data[g], data_size * sizeof(int));
     
//   //     if(g == 0)  /*Copy from Host to Device*/
//   //        cudaMemcpy(d_data[g], data, data_size * sizeof(int), cudaMemcpyHostToDevice);
//   // }
        
//   // ncclGroupStart();
 
//   // 		for(int g = 0; g < nGPUs; g++) {
//   // 	  	    cudaSetDevice(DeviceList[g]);
//   //   	  	    ncclBcast(d_data[g], data_size, ncclInt, 0, comms[g], s[g]); /*Broadcasting it to all*/
//   // 		}

//   // ncclGroupEnd();       

//   // for (int g = 0; g < nGPUs; g++) {
//   //     cudaSetDevice(DeviceList[g]);
//   //     ncclCommCuDevice( comms[g], &currRank);
//   //     printf("Current Rank : %d\n",currRank);
//   //     printf("\nThis is device %d\n", g);
//   //     kernel <<< 1 , data_size >>> (d_data[g],g);/*Call the CUDA Kernel: The code multiple the vector position per 2 on GPUs*/
//   //     cudaDeviceSynchronize();    
//   //     cout<<endl;         
//   // }

//   // printf("\n");

//   // for (int g = 0; g < nGPUs; g++) { /*Synchronizing CUDA Streams*/
//   //     cudaSetDevice(DeviceList[g]);
//   //     cudaStreamSynchronize(s[g]);
//   // }
 
//   // for(int g = 0; g < nGPUs; g++) {  /*Destroy CUDA Streams*/
//   //     cudaSetDevice(DeviceList[g]);
//   //     cudaStreamDestroy(s[g]);
//   // }

//   for(int g = 0; g < nGPUs; g++)    /*Finalizing NCCL*/
//      ncclCommDestroy(comms[g]);
  
//   /*Freeing memory*/
//   free(s);
//   free(data); 
//   free(DeviceList);

//   cudaFree(d_data);

//   return 0;

// }/*main*/

//
// Multiple Devices per Thread
// 


//
// execute command: mpirun -np 2 ./ex3gather.out 
//
/* output result:
myRank: 0 localRank: 0
myRank: 1 localRank: 1
myRank0 sendbuff[0]
 j: 0 hptr[i][j]: 0
 j: 1 hptr[i][j]: 1
 j: 2 hptr[i][j]: 2
myRank1 sendbuff[0]
 j: 0 hptr[i][j]: 0
 j: 1 hptr[i][j]: 1
 j: 2 hptr[i][j]: 2
Root is:0 ncclgather result is :
 j: 0 hptr[i][j]: 0
 j: 1 hptr[i][j]: 1
 j: 2 hptr[i][j]: 2
 j: 3 hptr[i][j]: 0
 j: 4 hptr[i][j]: 1
 j: 5 hptr[i][j]: 2
[MPI Rank 0] Success 
[MPI Rank 1] Success
*/
#include <stdio.h>
#include "cuda_runtime.h"
#include "nccl.h"
#include "mpi.h"
// #include "ncclEnhance.h"
#include <unistd.h>
#include <stdint.h>
#include <iostream>

#define MPICHECK(cmd) do {                          \
  int e = cmd;                                      \
  if( e != MPI_SUCCESS ) {                          \
    printf("Failed: MPI error %s:%d '%d'\n",        \
        __FILE__,__LINE__, e);   \
    exit(EXIT_FAILURE);                             \
  }                                                 \
} while(0)


#define CUDACHECK(cmd) do {                         \
  cudaError_t e = cmd;                              \
  if( e != cudaSuccess ) {                          \
    printf("Failed: Cuda error %s:%d '%s'\n",             \
        __FILE__,__LINE__,cudaGetErrorString(e));   \
    exit(EXIT_FAILURE);                             \
  }                                                 \
} while(0)


#define NCCLCHECK(cmd) do {                         \
  ncclResult_t r = cmd;                             \
  if (r!= ncclSuccess) {                            \
    printf("Failed, NCCL error %s:%d '%s'\n",             \
        __FILE__,__LINE__,ncclGetErrorString(r));   \
    exit(EXIT_FAILURE);                             \
  }                                                 \
} while(0)


static uint64_t getHostHash(const char* string) {
  // Based on DJB2a, result = result * 33 ^ char
  uint64_t result = 5381;
  for (int c = 0; string[c] != '\0'; c++){
    result = ((result << 5) + result) ^ string[c];
  }
  return result;
}


static void getHostName(char* hostname, int maxlen) {
  gethostname(hostname, maxlen);
  for (int i=0; i< maxlen; i++) {
    if (hostname[i] == '.') {
        hostname[i] = '\0';
        return;
    }
  }
}

__global__ void  init(float *dptr,int myRank)
{
  int id = threadIdx.x;
  dptr[id] = id;
}


int main(int argc, char* argv[])
{
    //each process is using two GPUs
    int nDev = 1;
    int root = 0;
    int size = 3;

    int myRank, nRanks, localRank = 0;

    //initializing MPI
    MPICHECK(MPI_Init(&argc, &argv));
    MPICHECK(MPI_Comm_rank(MPI_COMM_WORLD, &myRank));
    MPICHECK(MPI_Comm_size(MPI_COMM_WORLD, &nRanks));

    if (myRank == 0)
    {
        std::cout << "================================================================"
                  << "\n    Executing " << argv[0] << " now!\n"
                  << "================================================================\n";
    }

    //calculating localRank which is used in selecting a GPU
    uint64_t hostHashs[nRanks];
    char hostname[1024];
    getHostName(hostname, 1024);
    hostHashs[myRank] = getHostHash(hostname);
    MPICHECK(MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, hostHashs, sizeof(uint64_t), MPI_BYTE, MPI_COMM_WORLD));
    for (int p = 0; p < nRanks; p++)
    {
      if (p == myRank)
        break;
      if (hostHashs[p] == hostHashs[myRank])
        localRank++;
    }
    std::cout << "myRank: " << myRank << " localRank: " << localRank << "\n";

    float **sendbuff = (float **)malloc(nDev * sizeof(float *));
    float **recvbuff = (float **)malloc(nDev * sizeof(float *));
    float **hptr = (float **)malloc(nDev * sizeof(float *));
    cudaStream_t *s = (cudaStream_t *)malloc(sizeof(cudaStream_t) * nDev);

    //picking GPUs based on localRank
    for (int i = 0; i < nDev; ++i)
    {
      CUDACHECK(cudaSetDevice(localRank * nDev + i)); // 给所有设备编号
      CUDACHECK(cudaMalloc(sendbuff + i, size * sizeof(float)));
      CUDACHECK(cudaMalloc(recvbuff + i, nDev * nRanks * size * sizeof(float)));
      CUDACHECK(cudaMemset(sendbuff[i], 1, size * sizeof(float)));
      CUDACHECK(cudaMemset(recvbuff[i], 0, size * sizeof(float)));
      CUDACHECK(cudaStreamCreate(s + i));
      hptr[i] = (float *)malloc(nDev * nRanks * size * sizeof(float));
  }


  ncclUniqueId id;
  ncclComm_t comms[nDev];


  //generating NCCL unique ID at one process and broadcasting it to all
  if (myRank == 0) ncclGetUniqueId(&id);
  MPICHECK(MPI_Bcast((void *)&id, sizeof(id), MPI_BYTE, 0, MPI_COMM_WORLD));

  //initializing NCCL, group API is required around ncclCommInitRank as it is
  //called across multiple GPUs in each thread/process
  NCCLCHECK(ncclGroupStart());
  for (int i = 0; i < nDev; i++)
  {
    CUDACHECK(cudaSetDevice(localRank * nDev + i));
    init<<<1, size>>>(sendbuff[i], myRank);
    NCCLCHECK(ncclCommInitRank(comms + i, nRanks * nDev, id, myRank * nDev + i));
    cudaMemcpy(hptr[i],sendbuff[i],size*sizeof(float),cudaMemcpyDeviceToHost);
    std::cout<<"myRank"<<myRank<<" sendbuff["<<i<<"]"<<"\n";
    for(int j=0;j<size;++j){
        std::cout<<" j: "<<j<<" hptr[i][j]: "<<hptr[i][j]<<"\n";
    }
  }
  NCCLCHECK(ncclGroupEnd());


  // gather Data
  for(int i=0;i<nDev;++i){
    NCCLGather(sendbuff[i], size, ncclFloat, recvbuff[i], size, ncclFloat, root, comms[i], s[i]);
  }

  for (int i = 0; i < nDev; ++i)
  {
    if(myRank * nDev + i==root)
    {
      cudaMemcpy(hptr[i], recvbuff[i], nDev * nRanks * size * sizeof(float), cudaMemcpyDeviceToHost);
      std::cout << "Root is:" << root << " ncclgather result is :\n";
      for (int j = 0; j < nRanks * nDev * size; ++j)
      {
        std::cout << " j: " << j << " hptr[i][j]: " << hptr[i][j] << "\n";
      }
    }
  }

  //synchronizing on CUDA stream to complete NCCL communication
  for (int i=0; i<nDev; i++)
      CUDACHECK(cudaStreamSynchronize(s[i]));


  //freeing device memory
  for (int i=0; i<nDev; i++) {
     CUDACHECK(cudaFree(sendbuff[i]));
     CUDACHECK(cudaFree(recvbuff[i]));
     free(hptr[i]);
  }


  //finalizing NCCL
  for (int i=0; i<nDev; i++) {
     ncclCommDestroy(comms[i]);
  }


  //finalizing MPI
  MPICHECK(MPI_Finalize());


  printf("[MPI Rank %d] Success \n", myRank);
  return 0;
}