/*  
	MFCFD is a 3D Computational Fluid Dynamics Solver based off q-LSKUM
    Copyright (C) 2022 
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/
// Check read/write permissions
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


    //CHECK IF MPI LIBRARY HAD CUDA SUPPORT
    if(myRank==0){
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

      cout<<"Total Number of Proccess are "<<nRanks<<endl;
    }
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
  
    if(myRank==0){
      cout<<"Reading from File\n";
    }
    //READ POINTS FOR EACH DEVICE FROM FILE
    fstream fin;
    fin.open("/home/nsm/3d-mfcfd/inputFiles/"+to_string(max_points)+"/filesFor"+to_string(nRanks)+"Devices/Device"+to_string(myRank)+".dat",ios::in);
    fin>>numDevices;
    fin>>local_points;
    int counter;
    localToGlobalIndex=new int [local_points];
    globalToGhostIndex=new int*[nRanks];
    for(int i=0;i<nRanks;i++)
    {
        globalToGhostIndex[i]=new int[max_points];
    }

    for(int i=0;i<local_points;i++){
        fin>>localToGlobalIndex[i];
        fin>>partVector[localToGlobalIndex[i]];
        globalToLocalIndex[localToGlobalIndex[i]]=i;
        numberOfPointsPerDevice++;
    }
    cout<<"SIZE is "<<numberOfPointsPerDevice * sizeof(splitPoints)<<endl;
    //Share Partition Value across Proccesses
    MPICHECK(MPI_Allreduce(MPI_IN_PLACE, &partVector, max_points, MPI_INT, MPI_SUM, MPI_COMM_WORLD));    

    //ALLOCATING MEMORY FOR POINTS
    splitPoint=new splitPoints[numberOfPointsPerDevice];

    if(myRank==0){
      cout<<"Determining Nature of Points\n";
    }
    //ASSIGNING POINTS FOR EACH DEVICE AND CALCULATING NATURE OF POINTS ON EACH PARTITION
    for(int i=0;i<local_points;i++){
      assign(splitPoint[i],localToGlobalIndex[i],myRank);
      findNatureOfLocalPoints(splitPoint[i]);
      // printf("%d\n",i);
    }
    allocateSizeForNatureOfLocalPoints();
    for(int i=0;i<local_points;i++){
      assignNatureOfLocalPoints(splitPoint[i],i);
    }
    //Initialising the Send Buffer
    sendBuffer=new transferPoints*[nRanks];
    int points_on_gpu_to_send_to;
    int total_points_to_send=0;
    for(int i=0;i<nRanks;i++){
        fin>>points_on_gpu_to_send_to;
        total_points_to_send+=points_on_gpu_to_send_to;
        sendBuffer[i]=new transferPoints[points_on_gpu_to_send_to];
    }
    int currDevice=0;
    int *sendPoints=new int[nRanks];
    for(int i=0;i<nRanks;i++){
        sendPoints[i]=0;
    }

    //Populating the Send Buffer with Data
    for(int i=0;i<total_points_to_send;++i){
        fin>>currDevice;
        fin>>setprecision(13)>>sendBuffer[currDevice][sendPoints[currDevice]].globalIndex>>sendBuffer[currDevice][sendPoints[currDevice]].x>>sendBuffer[currDevice][sendPoints[currDevice]].y>>sendBuffer[currDevice][sendPoints[currDevice]].z;
        int currIndex=globalToLocalIndex[sendBuffer[currDevice][sendPoints[currDevice]].globalIndex];
        splitPoint[currIndex].isGhost=true;
        splitPoint[currIndex].ghostIndex[splitPoint[currIndex].numberOfPartitionsToSendTo]=sendPoints[currDevice];
        splitPoint[currIndex].partitions[splitPoint[currIndex].numberOfPartitionsToSendTo]=currDevice;
        splitPoint[currIndex].numberOfPartitionsToSendTo++;
        globalToGhostIndex[currDevice][sendBuffer[currDevice][sendPoints[currDevice]].globalIndex]=sendPoints[currDevice];
        sendPoints[currDevice]++;   
    }
    fin.close();
    cout<<"Number of Points in Process "<<myRank<<" are: "<<numberOfPointsPerDevice<<endl;
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
    // //

    int totalPointsToSend=0;
    for(int i=0;i<nRanks;i++){
      totalPointsToSend+=sendPoints[i];
    }

    //Initialising and transfering memory to device Pointers
    splitPoints *splitPoint_d;

    int *globalToLocalIndex_temp;
    int **globalToGhostIndex_send,**globalToGhostIndex_receive;
    int **globalToGhostIndexSendPointer=(int**)malloc(sizeof(int*)*nRanks);
    int** globalToGhostIndexReceivePointer=(int**)malloc(sizeof(int*)*nRanks);
    int *partVector_d;

    transferPoints** sendBuffer_d,** receiveBuffer_d;
    transferPoints** sendPointer=(transferPoints**)malloc(sizeof(transferPoints*)*nRanks);
    transferPoints** receivePointer=(transferPoints**)malloc(sizeof(transferPoints*)*nRanks);
    // CUcontext * temp;
    // cuCtxCreate(temp,CU_CTX_SCHED_AUTO,localRank);
    if(myRank==0){
        cout<<"Setting Device"<<endl;
    }
    CUDACHECK(cudaSetDevice(localRank));

    //POINTER TO POINTER
    CUDACHECK(cudaMalloc(&sendBuffer_d,sizeof(transferPoints*)*nRanks));
    for(int i = 0; i < nRanks; i++)
    {
      transferPoints *darray;
      CUDACHECK(cudaMalloc(&darray, sizeof(transferPoints) * sendPoints[i]));
      CUDACHECK(cudaMemcpy(darray,sendBuffer[i], sendPoints[i] * sizeof(transferPoints), cudaMemcpyHostToDevice));
      CUDACHECK(cudaMemcpy(&sendBuffer_d[i], &darray, sizeof(transferPoints*), cudaMemcpyHostToDevice));
      memcpy(&sendPointer[i],&darray,sizeof(transferPoints*));
    }

    CUDACHECK(cudaMalloc(&receiveBuffer_d,sizeof(transferPoints*)*nRanks));
    for(int i = 0; i < nRanks; i++)
    {
      transferPoints *darray;
      CUDACHECK(cudaMalloc(&darray, sizeof(transferPoints) * receivePoints[i]));
      CUDACHECK(cudaMemcpy(&receiveBuffer_d[i], &darray, sizeof(transferPoints*), cudaMemcpyHostToDevice));
      memcpy(&receivePointer[i],&darray,sizeof(transferPoints*));
    }

    CUDACHECK(cudaMalloc(&globalToGhostIndex_send, nRanks * sizeof(int*)));
    for(int i = 0; i < nRanks; i++)
    {
      int *darray;
      CUDACHECK(cudaMalloc(&darray, sizeof(int) * max_points));
      CUDACHECK(cudaMemcpy(darray,globalToGhostIndex[i], sizeof(int) * max_points, cudaMemcpyHostToDevice));
      CUDACHECK(cudaMemcpy(&globalToGhostIndex_send[i], &darray, sizeof(int*), cudaMemcpyHostToDevice));
      memcpy(&globalToGhostIndexSendPointer[i],&darray,sizeof(int*));
    }

    CUDACHECK(cudaMalloc(&globalToGhostIndex_receive, nRanks * sizeof(int*)));
    for(int i = 0; i < nRanks; i++)
    {
        int *darray;
        CUDACHECK(cudaMalloc(&darray, sizeof(int) * max_points));
        CUDACHECK(cudaMemcpy(&globalToGhostIndex_receive[i], &darray, sizeof(int*), cudaMemcpyHostToDevice));
        memcpy(&globalToGhostIndexReceivePointer[i],&darray,sizeof(int*));
    }

    if(myRank==0){
        cout<<"Copying linear arrays"<<endl;
    }
    CUDACHECK(cudaMalloc(&splitPoint_d, numberOfPointsPerDevice * sizeof(splitPoints)));
    CUDACHECK(cudaMemcpy(splitPoint_d, splitPoint, numberOfPointsPerDevice * sizeof(splitPoints), cudaMemcpyHostToDevice));
    CUDACHECK(cudaMalloc(&globalToLocalIndex_temp, max_points * sizeof(int)));
    CUDACHECK(cudaMemcpy(globalToLocalIndex_temp, globalToLocalIndex, max_points * sizeof(int), cudaMemcpyHostToDevice));
    CUDACHECK(cudaMalloc(&partVector_d, max_points * sizeof(int)));
    CUDACHECK(cudaMemcpy(partVector_d, &partVector, max_points * sizeof(int), cudaMemcpyHostToDevice));
    
    // //
   
    
    if (myRank == 0) {
      cout<<"Getting NCCL Unique ID\n";
      NCCLCHECK(ncclGetUniqueId(&id));
    }
    MPICHECK(MPI_Barrier(MPI_COMM_WORLD));

    MPICHECK(MPI_Bcast((void *)&id, sizeof(id), MPI_BYTE, 0, MPI_COMM_WORLD));
    // //
    if(myRank==0){
      cout<<"ID Succesfully Broadcasted\n";
    }

    // Initialising NCCL Communicator and CUDA stream
    ncclComm_t comm;
    NCCLCHECK(ncclCommInitRank(&comm, nRanks, id, myRank));
    cudaStream_t stream;

    if(myRank==0){
      cout<<"Beginning Solver\n";
    }

    auto start = high_resolution_clock::now();
    // 
    fpi_solver_multi_nccl(splitPoint_d,localRank,sendBuffer_d,receiveBuffer_d,nRanks,myRank,sendPoints,receivePoints,comm,stream,sendPointer,receivePointer,globalToLocalIndex_temp,globalToGhostIndex_receive,globalToGhostIndexSendPointer,globalToGhostIndexReceivePointer,partVector_d);
    // 
    auto stop = high_resolution_clock::now();
    if(myRank==0){
      cout<<"Copying memory back to Host\n";
    }
    // CUDACHECK(cudaMemcpy(splitPoint, splitPoint_d, numberOfPointsPerDevice * sizeof(splitPoints), cudaMemcpyDeviceToHost));

    
    // TO COPY BACK THE SEND BUFFER TO HOST (POINTER TO POINTER METHOD)
    // transferPoints *darray;
    // for (int i = 0; i <nRanks; i++){
    //   cudaMalloc(&darray, sizeof(transferPoints) * sendPoints[i]);
    //   cudaMemcpy(&darray, &sendBuffer_d[i], sizeof(transferPoints*), cudaMemcpyDeviceToHost);
    //   cudaMemcpy(sendBuffer[i], darray, sizeof(transferPoints) * sendPoints[i], cudaMemcpyDeviceToHost);
    //   cudaFree(darray);
    // }

    // for (int i = 0; i <nRanks; i++){
    //   cudaMalloc(&darray, sizeof(transferPoints) * receivePoints[i]);
    //   cudaMemcpy(&darray, &receiveBuffer_d[i], sizeof(transferPoints*), cudaMemcpyDeviceToHost);
    //   cudaMemcpy(receiveBuffer[i], darray, sizeof(transferPoints) * receivePoints[i], cudaMemcpyDeviceToHost);
    // }

    if(myRank==0){
      cout<<"Deallocating memory and Destroying Communicators\n";
    }
    // 
    CUDACHECK(cudaFree(splitPoint_d));
    // 
    NCCLCHECK(ncclCommDestroy(comm));
    MPI_Finalize();
    //
    auto duration = duration_cast<microseconds>(stop - start);
    cout << "Done with process "<<myRank<< ". Time Taken by was:" << duration.count() / 1000000.0 << endl;
}
