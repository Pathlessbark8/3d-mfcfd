// // #include <iostream>
// // #include <fstream>
// // #include <random>
// // #include "split_fluxes_mod.h"
// // #include <iomanip>
// // #include "octant_fluxes_mod.h"
// // #include "wall_flux_dGxneg_mod.h"
// // #include "wall_flux_dGyneg_mod.h"
// // #include "wall_flux_dGxpos_mod.h"
// // #include "wall_flux_dGypos_mod.h"
// // #include "wall_flux_dGzneg_mod.h"
// // #include "point_preprocessor_mod.h"
// // #include "compute_conserved_vector_mod.h"
// // #include "timestep_delt_mod.h"
// // #include "generate_connectivity_mod.h"
// // #include "implicit_aliasing_mod.h"
// // #include "flux_residual_mod.h"
// // #include "initial_conditions_mod.h"
// // #include "interior_flux_dGxneg_mod.h"
// // #include "interior_flux_dGxpos_mod.h"
// // #include "interior_flux_dGyneg_mod.h"
// // #include "interior_flux_dGypos_mod.h"
// // #include <chrono>
// // #include <cuda_runtime.h>
// // #include "fpi_solver_mod.h"
// // #include <nccl.h>
// // #include <mpi.h>

// // // #include "cuPrintf.cu"
// // using namespace std;
// // using namespace std::chrono;

// // #define MPICHECK(cmd) do {                          \
// //   int e = cmd;                                      \
// //   if( e != MPI_SUCCESS ) {                          \
// //     printf("Failed: MPI error %s:%d '%d'\n",        \
// //         __FILE__,__LINE__, e);   \
// //     exit(EXIT_FAILURE);                             \
// //   }                                                 \
// // } while(0)


// // #define CUDACHECK(cmd) do {                         \
// //   cudaError_t e = cmd;                              \
// //   if( e != cudaSuccess ) {                          \
// //     printf("Failed: Cuda error %s:%d '%s'\n",             \
// //         __FILE__,__LINE__,cudaGetErrorString(e));   \
// //     exit(EXIT_FAILURE);                             \
// //   }                                                 \
// // } while(0)


// // #define NCCLCHECK(cmd) do {                         \
// //   ncclResult_t r = cmd;                             \
// //   if (r!= ncclSuccess) {                            \
// //     printf("Failed, NCCL error %s:%d '%s'\n",             \
// //         __FILE__,__LINE__,ncclGetErrorString(r));   \
// //     exit(EXIT_FAILURE);                             \
// //   }                                                 \
// // } while(0)


// // static uint64_t getHostHash(const char* string) {
// //   // Based on DJB2a, result = result * 33 ^ char
// //   uint64_t result = 5381;
// //   for (int c = 0; string[c] != '\0'; c++){
// //     result = ((result << 5) + result) ^ string[c];
// //   }
// //   return result;
// // }

// // int main(int argc, char* argv[])
// // {

// //     cout<<setprecision(13)<<scientific;

// //     read_input_point_data();
// //     initial_conditions();
// //     generate_split_stencils();
// //     //
// //     int myRank, nRanks, localRank = 0;
  
// //     //initializing MPI
// //     MPICHECK(MPI_Init(&argc, &argv));
// //     MPICHECK(MPI_Comm_rank(MPI_COMM_WORLD, &myRank));
// //     MPICHECK(MPI_Comm_size(MPI_COMM_WORLD, &nRanks));
// //     //
// //     fstream fin;
// //     fin.open("/home/anil/new_3d_code/3d-mfcfd/inputFiles/4DevicePartFile.dat",ios::in);
// //     fin>>numDevices;
// //     int counter;
// //     for(int i=0;i<max_points;i++){
// //         fin>>point.counter[i]>>partVector[i];
// //         if(partVector[i]/numDevices==myRank)
// //         {
// //           numberOfPointsPerDevice[partVector[i]%numDevices]++;
// //         }
// //     }
// //     fin.close();
// //     //
// //     cout<<numberOfPointsPerDevice<<" "<<numberOfPointsPerDevice[1]<<endl;
// //     //
// //     //
// //     //
// //     for(int i=0;i<numDevices;i++){
// //         splitPoint[i]=new splitPoints[numberOfPointsPerDevice[i]];
// //     }
// //     int *pointCounter=new int[numDevices];
// //     for(int i=0;i<numDevices;i++){
// //         pointCounter[i]=0;
// //     }
// //     cout<<"YO\n";
// //     for(int i=0;i<max_points;i++){
// //       if(partVector[i]/numDevices==myRank){
// //         assign(splitPoint[partVector[i]%numDevices][pointCounter[partVector[i]%numDevices]],i);
// //         pointCounter[partVector[i]%numDevices]++;
// //       }
// //     }
// //     //
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
// //     // //
// //     for (int i = 0; i < numDevices; ++i) {
// //         CUDACHECK(cudaSetDevice(i));
// //         CUDACHECK(cudaMalloc(splitPoint_d + i, numberOfPointsPerDevice[i] * sizeof(splitPoints)));
// //         CUDACHECK(cudaStreamCreate(&s[i]));
// //     }
// //     for (int i = 0; i < numDevices; ++i) {
// //         CUDACHECK(cudaSetDevice(i));
// //         CUDACHECK(cudaMemcpyAsync(splitPoint_d[i], splitPoint[i], numberOfPointsPerDevice[i] * sizeof(splitPoints), cudaMemcpyHostToDevice, s[i]));
// //     }
// //     //
// //     cout<<splitPoint[0][0].x<<" "<<splitPoint[1][1].x<<endl;
// //     cout<<"HI3\n";
// //     if (myRank == 0) ncclGetUniqueId(&id);
// //     MPICHECK(MPI_Bcast((void *)&id, sizeof(id), MPI_BYTE, 0, MPI_COMM_WORLD));
  
  
// //     // initializing NCCL, group API is required around ncclCommInitRank as it is
// //     // called across multiple GPUs in each thread/process
// //     // NCCLCHECK(ncclGroupStart());
// //     for (int i=0; i<numDevices; i++) {
// //        CUDACHECK(cudaSetDevice(i));
// //        NCCLCHECK(ncclCommInitRank(comms+i, nRanks*numDevices, id, myRank*numDevices + i));
// //        cout<<"Rank "<<myRank<<" Device "<<i<<endl;
// //     }
// //     // NCCLCHECK(ncclGroupEnd());

// //     int count;
// //     ncclCommCount(comms[0], &count);
// //     cout<<myRank <<" Commcount "<<count<<endl; 

// //     ncclCommCuDevice(comms[0], &count);
// //     cout<<myRank <<" CuDevice "<<count<<endl;

// //     ncclCommUserRank(comms[1], &count);
// //     cout<<myRank <<" Comm user "<<count<<endl;
// //     for (int i=0; i<numDevices; i++)
// //       CUDACHECK(cudaStreamSynchronize(s[i]));
// //     cout<<"HI4\n";
// //     // //
// //     fpi_solver_multi_nccl(splitPoint_d,s);
// //     // //
// //     cout<<"Copying memory back to Host\n";
// //     //
// //     for (int i = 0; i < numDevices; ++i) {
// //         CUDACHECK(cudaSetDevice(i));
// //         CUDACHECK(cudaMemcpyAsync(splitPoint[i], splitPoint_d[i], numberOfPointsPerDevice[i] * sizeof(splitPoints), cudaMemcpyDeviceToHost, s[i]));
// //     }
// //     cout<<"Deallocating memory and Destroying Communicators\n";
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
// //     MPI_Finalize();
// //     cout << "Done with process "<<myRank<<endl;
// // }

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
#include <unistd.h>
#include <stdint.h>
#include <stdlib.h>
#include "mpi-ext.h"

// #include "cuPrintf.cu"
using namespace std;
using namespace std::chrono;



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

int main(int argc, char* argv[])
{

  printf("Compile time check:\n");
  #if defined(MPIX_CUDA_AWARE_SUPPORT) && MPIX_CUDA_AWARE_SUPPORT
      printf("This MPI library has CUDA-aware support.\n", MPIX_CUDA_AWARE_SUPPORT);
  #elif defined(MPIX_CUDA_AWARE_SUPPORT) && !MPIX_CUDA_AWARE_SUPPORT
      printf("This MPI library does not have CUDA-aware support.\n");
  #else
      printf("This MPI library cannot determine if there is CUDA-aware support.\n");
  #endif /* MPIX_CUDA_AWARE_SUPPORT */
   
      printf("Run time check:\n");
  #if defined(MPIX_CUDA_AWARE_SUPPORT)
      if (1 == MPIX_Query_cuda_support()) {
          printf("This MPI library has CUDA-aware support.\n");
      } else {
          printf("This MPI library does not have CUDA-aware support.\n");
      }
  #else /* !defined(MPIX_CUDA_AWARE_SUPPORT) */
      printf("This MPI library cannot determine if there is CUDA-aware support.\n");
  #endif /* MPIX_CUDA_AWARE_SUPPORT */


    //MAIN CODE BEGINS
    cout<<setprecision(13)<<scientific;

    read_input_point_data();
    initial_conditions();
    generate_split_stencils();
    //
    int myRank, nRanks, localRank = 0;
  
    //initializing MPI
    MPICHECK(MPI_Init(&argc, &argv));
    MPICHECK(MPI_Comm_rank(MPI_COMM_WORLD, &myRank));
    MPICHECK(MPI_Comm_size(MPI_COMM_WORLD, &nRanks));
    //
    cout<<nRanks<<endl;
    //
    // HASHING HOSTNAME TO GET LOCALRANKS
    uint64_t hostHashs[nRanks];
    char hostname[1024];
    getHostName(hostname, 1024);
    hostHashs[myRank] = getHostHash(hostname);
    MPICHECK(MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, hostHashs, sizeof(uint64_t), MPI_BYTE, MPI_COMM_WORLD));
    for (int p=0; p<nRanks; p++) {
      if (p == myRank) break;
      if (hostHashs[p] == hostHashs[myRank]) localRank++;
    }
  

    //READ POINTS FOR EACH DEVICE FROM FILE
    fstream fin;
    fin.open("/home/anil/new_3d_code/3d-mfcfd/inputFiles/filesFor"+to_string(nRanks)+"Devices/Device"+to_string(myRank)+".dat",ios::in);
    fin>>numDevices;
    fin>>local_points;
    int counter;
    localToGlobalIndex=new int [local_points];
    globalToGhostIndex=new int*[nRanks];
    for(int i=0;i<nRanks;i++)
    {
        globalToGhostIndex[i]=new int[max_points];
    }
    // partVector=new int [max_points];
    // ghostToGlobalIndex=new int *[nRanks];
    for(int i=0;i<local_points;i++){
        fin>>localToGlobalIndex[i];
        fin>>partVector[localToGlobalIndex[i]];
        globalToLocalIndex[localToGlobalIndex[i]]=i;
        // if(localToGlobalIndex[i]==442481){
        //   cout<<globalToLocalIndex[localToGlobalIndex[i]]<<" serhghtsdfghs \n";
        // }
        // if(partVector[i]==myRank)
        // {
        //   numberOfPointsPerDevice++;
        // }
        numberOfPointsPerDevice++;
    }
    MPICHECK(MPI_Allreduce(MPI_IN_PLACE, &partVector, max_points, MPI_INT, MPI_SUM, MPI_COMM_WORLD));
    // cout<<"CHECK : "<<globalToLocalIndex[432]<<endl;
    //
    ////ALLOCATING MEMORY FOR POINTS
    splitPoint=new splitPoints[numberOfPointsPerDevice];
    //
    //ASSIGNING POINTS FOR EACH DEVICE AND CALCULATING NATURE OF POINTS ON EACH PARTITION
    cout<<"Initialising....\n";
    for(int i=0;i<local_points;i++){
      // if(localToGlobalIndex[i]==442481){
      //   cout<<"Local Index is "<<i<<endl;
      // }
      assign(splitPoint[i],localToGlobalIndex[i],myRank);
      // if(localToGlobalIndex[i]==442481){
      //     cout<<"BIUHASDUCHUISHC "<<splitPoint[199977].globalIndex<<" "<<splitPoint[199977].x<<" "<<splitPoint[199977].y<<" "<<splitPoint[285282].z<<endl;
      // }
      findNatureOfLocalPoints(splitPoint[i]);
    }
    // cout<<"BIUHASDUC÷HUISHC "<<splitPoint[199977].globalIndex<<" "<<splitPoint[199977].x<<" "<<splitPoint[199977].y<<" "<<splitPoint[199977].z<<endl;
    allocateSizeForNatureOfLocalPoints();
    for(int i=0;i<local_points;i++){
      assignNatureOfLocalPoints(splitPoint[i],i);
    }

    // if(myRank==1){
    //   cout<<"Number Of Local Nbhs : "<<splitPoint[6081].numberOfGhostxposNbhs<<" "<<splitPoint[6081].numberOfLocalxposNbhs<<endl;
    // }

    // cout<<"Verification : "<<myRank<<" : "<<interiorPointsLocal<<" "<<wallPointsLocal<<" "<<outerPointsLocal<<" "<<interiorPointsLocal+wallPointsLocal+outerPointsLocal<<endl;
    //
    // cout<< " CHECK : "<< splitPoint[258].globalIndex<<endl;
    //
    sendBuffer=new transferPoints*[nRanks];
    int points_on_gpu_to_send_to;
    int total_points_to_send=0;
    for(int i=0;i<nRanks;i++){
        fin>>points_on_gpu_to_send_to;
        // ghostToGlobalIndex[i]=new int[points_on_gpu_to_send_to];
        total_points_to_send+=points_on_gpu_to_send_to;
        cout<<points_on_gpu_to_send_to<<endl;
        sendBuffer[i]=new transferPoints[points_on_gpu_to_send_to];
        // if(myRank==0)
        //   cout<<i<<" "<<points_on_gpu_to_send_to<<endl;
    }
    
    int currDevice=0;
    int *sendPoints=new int[nRanks];
    for(int i=0;i<nRanks;i++){
        sendPoints[i]=0;
    }
    cout<<endl;

    for(int i=0;i<total_points_to_send;++i){
        fin>>currDevice;
        fin>>setprecision(13)>>sendBuffer[currDevice][sendPoints[currDevice]].globalIndex>>sendBuffer[currDevice][sendPoints[currDevice]].x>>sendBuffer[currDevice][sendPoints[currDevice]].y>>sendBuffer[currDevice][sendPoints[currDevice]].z;
        int currIndex=globalToLocalIndex[sendBuffer[currDevice][sendPoints[currDevice]].globalIndex];
        splitPoint[currIndex].isGhost=true;
        splitPoint[currIndex].ghostIndex[splitPoint[currIndex].numberOfPartitionsToSendTo]=sendPoints[currDevice];
        splitPoint[currIndex].partitions[splitPoint[currIndex].numberOfPartitionsToSendTo]=currDevice;
        splitPoint[currIndex].numberOfPartitionsToSendTo++;
        // if(sendBuffer[currDevice][sendPoints[currDevice]].globalIndex==430334){
        //   cout<<"kjhedfjkehfjewhjf "<<sendPoints[currDevice];
        // }
        // if(sendBuffer[currDevice][sendPoints[currDevice]].globalIndex == 529172){
        //   cout<<"sioefhisdhfisdhflsdhflkshdlfksdlkf "<<sendPoints[currDevice]<<endl;
        // }
        // globalToGhostIndex[sendBuffer[currDevice][sendPoints[currDevice]].globalIndex]=sendPoints[currDevice];
        // if(sendBuffer[currDevice][sendPoints[currDevice]].globalIndex == 529172){
        //   cout<<"sioefhisdhfisdhflsdhflkshdlfksdlkf "<<globalToGhostIndex[sendBuffer[currDevice][sendPoints[currDevice]].globalIndex]<<endl;
        //   cout<<sendBuffer[currDevice][sendPoints[currDevice]].x<<" "<<sendBuffer[currDevice][sendPoints[currDevice]].y<<" "<<sendBuffer[currDevice][sendPoints[currDevice]].z<<endl;
        // }
        globalToGhostIndex[currDevice][sendBuffer[currDevice][sendPoints[currDevice]].globalIndex]=sendPoints[currDevice];
        // if(sendBuffer[currDevice][sendPoints[currDevice]].globalIndex == 333830){
        //   cout<<"sioefhisdhfisdhflsdhflkshdlfksdlkf "<<globalToGhostIndex[currDevice][sendBuffer[currDevice][sendPoints[currDevice]].globalIndex]<<endl;
        //   cout<<sendBuffer[currDevice][sendPoints[currDevice]].x<<" "<<sendBuffer[currDevice][sendPoints[currDevice]].y<<" "<<sendBuffer[currDevice][sendPoints[currDevice]].z<<endl;
        // }
        sendPoints[currDevice]++;
        
    }
    // if(myRank==0){
    //   cout<<"o;iwhfis;dhafidsahfisdhflshdfkdash "<<globalToGhostIndex[1][333830]<<endl;
    // }
    fin.close();
    // if(myRank==0){
    //   cout<<"o;iwhfis;dhafidsahfisdhflshdfkdash "<<globalToGhostIndex[1][333830]<<endl;
    // }
    
    cout<<"Number of Points in Process "<<myRank<<" are: "<<numberOfPointsPerDevice<<endl;
    //
    //
    cout<<"HI1\n";
    // //
    //Sharing Size across All Processes
    int *receivePoints=new int[nRanks];
    for(int i=0;i<nRanks;i++){
        receivePoints[i]=0;
    }
    for(int i=0;i<nRanks;i++){
      if(i!=myRank){
        MPI_Request request;
        MPICHECK(MPI_Isend(sendPoints+i , 1, MPI_INT, i, 0, MPI_COMM_WORLD,&request));
      }
    }
    for(int i=0;i<nRanks;i++){
      if(i!=myRank){
        // MPICHECK(MPI_Recv(receivePoints+i , 1, MPI_INT, MPI , 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE));
        MPICHECK(MPI_Recv(receivePoints+i , 1, MPI_INT, i, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE));
      }
    }
    // cout<<myRank<<" "<<sendPoints[0]<<" "<<sendPoints[1]<<" "<<sendPoints[2]<<" "<<sendPoints[3]<<endl;
    // cout<<myRank<<" "<<receivePoints[0]<<" "<<receivePoints[1]<<" "<<receivePoints[2]<<" "<<receivePoints[3]<<endl;

    cout<<myRank<<" "<<sendPoints[0]<<" "<<sendPoints[1]<<endl;
    cout<<myRank<<" "<<receivePoints[0]<<" "<<receivePoints[1]<<endl;

    receiveBuffer=new transferPoints*[nRanks];
    for(int i=0;i<nRanks;i++){
      receiveBuffer[i]=new transferPoints[receivePoints[i]];
    }


    //CREATE MPI STRUCTURE TO TRANSFER POINTS DATA TO OTHER PROCESSES
    // const int nitems=5;
    // int          blocklengths[5] = {1,1,1,1,5};
    // MPI_Datatype types[5] = {MPI_INT, MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE};
    // MPI_Datatype mpi_transferPoints;
    // MPI_Aint     offsets[5];

    // offsets[0] = offsetof(transferPoints, counter);
    // offsets[1] = offsetof(transferPoints, x);
    // offsets[2] = offsetof(transferPoints, y);
    // offsets[3] = offsetof(transferPoints, z);
    // offsets[4] = offsetof(transferPoints, q);


    // MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_transferPoints);
    // MPI_Type_commit(&mpi_transferPoints);

    // for(int i=0;i<nRanks;i++){
    //   if(i!=myRank){
    //     MPICHECK(MPI_Send(sendBuffer[i] , sendPoints[i], MPI_BYTE, i, 0, MPI_COMM_WORLD));
    //     MPICHECK(MPI_Recv(receiveBuffer[myRank] , receivePoints[myRank], MPI_BYTE, myRank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE));
    //   }
    // }

    // if(myRank==0){
    // cout<<myRank<<" "<<sendBuffer[1][0].x<<" "<<sendPoints[1]<<endl;
    // // cout<<myRank<<" "<<receiveBuffer[0]<<" "<<receivePoints[1]<<endl;
    // }
    cout<<"HI2\n";
    // cout<<partVector
    cout<<numDevices<<endl;
    // //
    int totalPointsToSend=0;
    for(int i=0;i<nRanks;i++){
      totalPointsToSend+=sendPoints[i];
    }
    cout<<"localRank: "<<localRank<<endl;
    cout<<myRank<<" "<<totalPointsToSend<<endl;


    splitPoints *splitPoint_d;
    int *globalToLocalIndex_temp;
    int **globalToGhostIndex_send,**globalToGhostIndex_receive;
    int **globalToGhostIndexSendPointer=(int**)malloc(sizeof(int*)*nRanks);
    int** globalToGhostIndexReceivePointer=(int**)malloc(sizeof(int*)*nRanks);

    int *partVector_d;

    transferPoints** sendBuffer_d,** receiveBuffer_d;
    transferPoints** sendPointer=(transferPoints**)malloc(sizeof(transferPoints*)*nRanks);
    transferPoints** receivePointer=(transferPoints**)malloc(sizeof(transferPoints*)*nRanks);
    CUDACHECK(cudaSetDevice(localRank));

    //POINTER TO POINTER
    CUDACHECK(cudaMalloc(&sendBuffer_d,sizeof(transferPoints*)*nRanks));
    for(int i = 0; i < nRanks; i++)
    {
      transferPoints *darray;
      cudaMalloc(&darray, sizeof(transferPoints) * sendPoints[i]);
      cudaMemcpy(darray,sendBuffer[i], sendPoints[i] * sizeof(transferPoints), cudaMemcpyHostToDevice);
      cudaMemcpy(&sendBuffer_d[i], &darray, sizeof(transferPoints*), cudaMemcpyHostToDevice);
      memcpy(&sendPointer[i],&darray,sizeof(transferPoints*));
      // cudaFree(darray);
      // cudaMalloc(&sendBuffer_d[i],sizeof(transferPoints) * sendPoints[i]);
      // cudaMemcpy(sendBuffer_d[i],sendBuffer[i], sendPoints[i] * sizeof(transferPoints), cudaMemcpyHostToDevice);
    }

    // transferPoints** receiveBuffer_d;
    CUDACHECK(cudaMalloc(&receiveBuffer_d,sizeof(transferPoints*)*nRanks));
    for(int i = 0; i < nRanks; i++)
    {
      transferPoints *darray;
      cudaMalloc(&darray, sizeof(transferPoints) * receivePoints[i]);
      // cout<<"darray "<<&darray<<endl;
      // cudaMemcpy(darray,sendBuffer[i], sendPoints[i] * sizeof(transferPoints), cudaMemcpyHostToDevice);
      cudaMemcpy(&receiveBuffer_d[i], &darray, sizeof(transferPoints*), cudaMemcpyHostToDevice);
      memcpy(&receivePointer[i],&darray,sizeof(transferPoints*));
      // cout<<"darray "<<&receivePointer[i]<<endl;
      // cudaFree(darray);
    }

    CUDACHECK(cudaMalloc(&globalToGhostIndex_send, nRanks * sizeof(int*)));
    for(int i = 0; i < nRanks; i++)
    {
      int *darray;
      cudaMalloc(&darray, sizeof(int) * max_points);
      cudaMemcpy(darray,globalToGhostIndex[i], sizeof(int) * max_points, cudaMemcpyHostToDevice);
      cudaMemcpy(&globalToGhostIndex_send[i], &darray, sizeof(int*), cudaMemcpyHostToDevice);
      memcpy(&globalToGhostIndexSendPointer[i],&darray,sizeof(int*));
    }
   

    CUDACHECK(cudaMalloc(&globalToGhostIndex_receive, nRanks * sizeof(int*)));
    for(int i = 0; i < nRanks; i++)
    {
        int *darray;
        cudaMalloc(&darray, sizeof(int) * max_points);
        cudaMemcpy(&globalToGhostIndex_receive[i], &darray, sizeof(int*), cudaMemcpyHostToDevice);
        memcpy(&globalToGhostIndexReceivePointer[i],&darray,sizeof(int*));
    }

    // CUDACHECK(cudaStreamCreate(&s[0]));
        // CUDACHECK(cudaSetDevice(i));
    // if(myRank==1){
    //   cout<<"o;iwhfis;dhafidsahfisdhflshdfkdash "<<globalToGhostIndex[529172]<<endl;
    // }
    CUDACHECK(cudaMalloc(&splitPoint_d, numberOfPointsPerDevice * sizeof(splitPoints)));
    CUDACHECK(cudaMemcpy(splitPoint_d, splitPoint, numberOfPointsPerDevice * sizeof(splitPoints), cudaMemcpyHostToDevice));
    CUDACHECK(cudaMalloc(&globalToLocalIndex_temp, max_points * sizeof(int)));
    CUDACHECK(cudaMemcpy(globalToLocalIndex_temp, globalToLocalIndex, max_points * sizeof(int), cudaMemcpyHostToDevice));
    CUDACHECK(cudaMalloc(&partVector_d, max_points * sizeof(int)));
    CUDACHECK(cudaMemcpy(partVector_d, &partVector, max_points * sizeof(int), cudaMemcpyHostToDevice));
    
    // //
   
    // cout<<"STUCK "<<myRank<<endl;
    cout<<"HI3\n";
    if (myRank == 0) ncclGetUniqueId(&id);
    MPI_Barrier(MPI_COMM_WORLD);
    cout<<"Hiiiii\n";

    MPICHECK(MPI_Bcast((void *)&id, sizeof(id), MPI_BYTE, 0, MPI_COMM_WORLD));
    // //
    cout<<"Hiiiii\n";
    // //
    ncclComm_t comm;
    NCCLCHECK(ncclCommInitRank(&comm, nRanks, id, myRank));
    cout<<"Rank "<<myRank<<" Device "<<localRank<<endl;

    cudaStream_t stream;
    // cudaStreamCreateWithFlags(stream, cudaStreamNonBlocking);
    // int count;
    // ncclCommCount(comm, &count);
    // cout<<myRank <<" Commcount "<<count<<endl; 

    // ncclCommCuDevice(comm, &count);
    // cout<<myRank <<" CuDevice "<<count<<endl;

    // ncclCommUserRank(comm, &count);
    // cout<<myRank <<" Comm user "<<count<<endl;
    // for (int i=0; i<numDevices; i++)
    //   CUDACHECK(cudaDeviceSynchronize());
    cout<<"HI4\n";
    auto start = high_resolution_clock::now();
    // // //
    // cout<<"BIUHASDUCHUISHC "<<splitPoint[199977].globalIndex<<" "<<splitPoint[199977].x<<" "<<splitPoint[199977].y<<" "<<splitPoint[199977].z<<endl;

    fpi_solver_multi_nccl(splitPoint_d,localRank,sendBuffer_d,receiveBuffer_d,nRanks,myRank,sendPoints,receivePoints,comm,stream,sendPointer,receivePointer,globalToLocalIndex_temp,globalToGhostIndex_receive,globalToGhostIndexSendPointer,globalToGhostIndexReceivePointer,partVector_d);
    // //
    auto stop = high_resolution_clock::now();
    cout<<"Copying memory back to Host\n";
    //
    CUDACHECK(cudaMemcpy(splitPoint, splitPoint_d, numberOfPointsPerDevice * sizeof(splitPoints), cudaMemcpyDeviceToHost));

    
    // TO COPY BACK THE SEND BUFFER TO HOST (POINTER TO POINTER METHOD)
    transferPoints *darray;
    for (int i = 0; i <nRanks; i++){
      cudaMalloc(&darray, sizeof(transferPoints) * sendPoints[i]);
      cudaMemcpy(&darray, &sendBuffer_d[i], sizeof(transferPoints*), cudaMemcpyDeviceToHost);
      cudaMemcpy(sendBuffer[i], darray, sizeof(transferPoints) * sendPoints[i], cudaMemcpyDeviceToHost);
      cudaFree(darray);
    }

    // transferPoints *darray;
    for (int i = 0; i <nRanks; i++){
      cudaMalloc(&darray, sizeof(transferPoints) * receivePoints[i]);
      cudaMemcpy(&darray, &receiveBuffer_d[i], sizeof(transferPoints*), cudaMemcpyDeviceToHost);
      cudaMemcpy(receiveBuffer[i], darray, sizeof(transferPoints) * receivePoints[i], cudaMemcpyDeviceToHost);
      // cudaFree(darray);
    }

    // int *tempArray;
    // for (int i = 0; i <nRanks; i++){
    //   cudaMalloc(&tempArray, sizeof(int) * max_points);
    //   cudaMemcpy(&darray, &receiveBuffer_d[i], sizeof(transferPoints*), cudaMemcpyDeviceToHost);
    //   cudaMemcpy(receiveBuffer[i], darray, sizeof(transferPoints) * receivePoints[i], cudaMemcpyDeviceToHost);
    //   // cudaFree(darray);
    // }

    // if(myRank==0){
    // cout<<myRank<<" "<<sendBuffer[1][1].globalIndex<<" "<<sendBuffer[1][1].x<<endl;
    // for(int r=0;r<5;r++){
    //   cout<<"TEST "<<sendBuffer[1][1].q[r]<<endl;
    // }
    // }

    cout<<"Deallocating memory and Destroying Communicators\n";
    // if(myRank==1){
    // cout<<myRank<<" "<<receiveBuffer[0][1].globalIndex<<" "<<receiveBuffer[0][1].x<<endl;
    // for(int r=0;r<5;r++){
    //   cout<<"TEST recv buffer "<<receiveBuffer[0][1].q[r]<<endl;
    // }
    // }
    // cout<< " CHECK : "<< splitPoint[258].globalIndex<<endl;
    cout<<"HI7\n";
    // if(myRank==1){
    //   cout<<"kahfhdal"<<endl;
    //   cout<<myRank<<" "<<splitPoint[6081].globalIndex<<endl;
    //   // cout<<myRank<<" "<<splitPoint[6081].numberOfLocalxnegNbhs<<endl;
    //   cout<<myRank<<" "<<splitPoint[6081].numberOfGhostxnegNbhs<<endl;
    //   for(int r=0;r<5;r++){
    //     cout<<"TEST "<<splitPoint[6081].dq[0][r]<<endl;
    //   }
    // }
    // // //
    CUDACHECK(cudaFree(splitPoint_d));
    // // //
    ncclCommDestroy(comm);
    MPI_Finalize();
    cout << "Done with process "<<myRank<<endl;
    
    auto duration = duration_cast<microseconds>(stop - start);
    cout << "Time Taken :" << duration.count() / 1000000.0 << endl;
}
