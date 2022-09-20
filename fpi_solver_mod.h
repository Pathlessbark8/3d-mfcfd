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
#pragma once

#include "data_structure_mod.h"
#include "compute_conserved_vector_mod.h"
#include "q_variables_mod.h"
#include "q_derivatives_mod.h"
#include "timestep_delt_mod.h"
#include "flux_residual_mod.h"
#include "state_update_mod.h"
#include "cuda.h"
#include <cuda_runtime.h>
#include <cuda_profiler_api.h>
#include <thrust/reduce.h>
#include <thrust/execution_policy.h>
#include <nccl.h>
#include <mpi.h>
#include <nvToolsExt.h>
#include <nvToolsExtCuda.h>

#define NCCLCHECK(cmd)                                         \
    do                                                         \
    {                                                          \
        ncclResult_t r = cmd;                                  \
        if (r != ncclSuccess)                                  \
        {                                                      \
            printf("Failed, NCCL error %s:%d '%s'\n",          \
                   __FILE__, __LINE__, ncclGetErrorString(r)); \
            exit(EXIT_FAILURE);                                \
        }                                                      \
    } while (0)

#define CUDACHECK(cmd)                                         \
    do                                                         \
    {                                                          \
        cudaError_t e = cmd;                                   \
        if (e != cudaSuccess)                                  \
        {                                                      \
            printf("Failed: Cuda error %s:%d '%s'\n",          \
                   __FILE__, __LINE__, cudaGetErrorString(e)); \
            exit(EXIT_FAILURE);                                \
        }                                                      \
    } while (0)

#define MPICHECK(cmd)                                \
    do                                               \
    {                                                \
        int e = cmd;                                 \
        if (e != MPI_SUCCESS)                        \
        {                                            \
            printf("Failed: MPI error %s:%d '%d'\n", \
                   __FILE__, __LINE__, e);           \
            exit(EXIT_FAILURE);                      \
        }                                            \
    } while (0)

void fpi_solver(int t)
//
{
    eval_q_variables();
    eval_q_derivatives();
    //
    timestep_delt();
    eval_flux_residual();
    state_update();
    //
    if (t <= 2)
    {
        res_old = res_new;
        residue = 0.00;
    }
    else
    {
        residue = log10(res_new / res_old);
    }
    //
}

void fpi_solver_cuda(points *point_d, cudaStream_t stream)
{

    dim3 threads(threads_per_block, 1, 1);
    dim3 grid(ceil((max_points / threads.x) + 1), 1, 1);

    fstream fout, fout_1;
    fout.open("output_40k.dat", ios::out);
    fout << setprecision(13);
    fout_1.open("residue.dat", ios::out);
    fout_1 << setprecision(13);
    int *wall_points_index_d;
    unsigned long long wall_size = sizeof(wall_points_index);
    cudaMalloc(&wall_points_index_d, wall_size);
    cudaMemcpy(wall_points_index_d, &wall_points_index, wall_size, cudaMemcpyHostToDevice);

    int *outer_points_index_d;
    unsigned long long outer_size = sizeof(outer_points_index);
    cudaMalloc(&outer_points_index_d, outer_size);
    cudaMemcpy(outer_points_index_d, &outer_points_index, outer_size, cudaMemcpyHostToDevice);

    int *interior_points_index_d;
    unsigned long long interior_size = sizeof(interior_points_index);
    cudaMalloc(&interior_points_index_d, interior_size);
    cudaMemcpy(interior_points_index_d, &interior_points_index, interior_size, cudaMemcpyHostToDevice);

    int *supersonic_inlet_points_index_d;
    unsigned long long supersonic_inlet_points_size = sizeof(supersonic_inlet_points_index);
    cudaMalloc(&supersonic_inlet_points_index_d, supersonic_inlet_points_size);
    cudaMemcpy(supersonic_inlet_points_index_d, &supersonic_inlet_points_index, supersonic_inlet_points_size, cudaMemcpyHostToDevice);

    int *supersonic_outlet_points_index_d;
    unsigned long long supersonic_outlet_points_size = sizeof(supersonic_outlet_points_index);
    cudaMalloc(&supersonic_outlet_points_index_d, supersonic_outlet_points_size);
    cudaMemcpy(supersonic_outlet_points_index_d, &supersonic_outlet_points_index, supersonic_outlet_points_size, cudaMemcpyHostToDevice);

    double *sum_res_sqr_d = NULL;
    cudaMalloc(&sum_res_sqr_d, sizeof(double) * (max_points));
    cudaMemset(sum_res_sqr_d, 0, sizeof(double) * (max_points));

    unsigned long long point_size = sizeof(point);
    cudaDeviceSynchronize();
    for (int t = 1; t <= max_iters; t++)
    {
        // cudaMemcpy(point_d, &point, point_size, cudaMemcpyHostToDevice);
        eval_q_variables_cuda<<<grid, threads>>>(*point_d);
        eval_q_derivatives_cuda<<<grid, threads>>>(*point_d, power);
        for (int r = 0; r < inner_iterations; r++)
        {
            q_inner_loop_cuda<<<grid, threads>>>(*point_d, power);
            update_inner_loop_cuda<<<grid, threads>>>(*point_d);
        }

        timestep_delt_cuda<<<grid, threads>>>(*point_d, CFL);

        wall_dGx_pos_cuda<<<grid, threads>>>(*point_d, power, VL_CONST, pi, wall_points, wall_points_index_d);
        wall_dGx_neg_cuda<<<grid, threads>>>(*point_d, power, VL_CONST, pi, wall_points, wall_points_index_d);
        wall_dGy_pos_cuda<<<grid, threads>>>(*point_d, power, VL_CONST, pi, wall_points, wall_points_index_d);
        wall_dGy_neg_cuda<<<grid, threads>>>(*point_d, power, VL_CONST, pi, wall_points, wall_points_index_d);
        wall_dGz_neg_cuda<<<grid, threads>>>(*point_d, power, VL_CONST, pi, wall_points, wall_points_index_d);
        // // //
        outer_dGx_pos_cuda<<<grid, threads>>>(*point_d, power, VL_CONST, pi, outer_points, outer_points_index_d);
        outer_dGx_neg_cuda<<<grid, threads>>>(*point_d, power, VL_CONST, pi, outer_points, outer_points_index_d);
        outer_dGy_pos_cuda<<<grid, threads>>>(*point_d, power, VL_CONST, pi, outer_points, outer_points_index_d);
        outer_dGy_neg_cuda<<<grid, threads>>>(*point_d, power, VL_CONST, pi, outer_points, outer_points_index_d);
        outer_dGz_pos_cuda<<<grid, threads>>>(*point_d, power, VL_CONST, pi, outer_points, outer_points_index_d);
        //
        interior_dGx_pos_cuda<<<grid, threads>>>(*point_d, power, VL_CONST, pi, interior_points, interior_points_index_d);
        interior_dGx_neg_cuda<<<grid, threads>>>(*point_d, power, VL_CONST, pi, interior_points, interior_points_index_d);
        interior_dGy_pos_cuda<<<grid, threads>>>(*point_d, power, VL_CONST, pi, interior_points, interior_points_index_d);
        interior_dGy_neg_cuda<<<grid, threads>>>(*point_d, power, VL_CONST, pi, interior_points, interior_points_index_d);
        interior_dGz_pos_cuda<<<grid, threads>>>(*point_d, power, VL_CONST, pi, interior_points, interior_points_index_d);
        interior_dGz_neg_cuda<<<grid, threads>>>(*point_d, power, VL_CONST, pi, interior_points, interior_points_index_d);
        // //
        //
        //
        // //
        state_update_wall<<<grid, threads>>>(*point_d, wall_points, wall_points_index_d, sum_res_sqr_d);
        state_update_outer<<<grid, threads>>>(*point_d, outer_points, outer_points_index_d, u1_inf, u2_inf, u3_inf, rho_inf, pi, pr_inf);
        state_update_interior<<<grid, threads>>>(*point_d, interior_points, interior_points_index_d, sum_res_sqr_d);

        // cudaMemcpy(&point, point_d, point_size, cudaMemcpyDeviceToHost);
        // state_update_supersonic_outlet();
        // state_update_supersonic_inlet();
        cudaDeviceSynchronize();
        sum_res_sqr = thrust::reduce(thrust::cuda::par.on(stream), sum_res_sqr_d, sum_res_sqr_d + max_points, (double)0.0, thrust::plus<double>());
        res_new = sqrt(sum_res_sqr) / max_points;
        if (t <= 2)
        {
            res_old = res_new;
            residue = 0.00;
        }
        else
        {
            residue = log10(res_new / res_old);
        }
        cout << t << " " << res_new << " " << residue << " " << sum_res_sqr << endl;
        fout_1 << t << " " << res_new << " " << residue << endl;
        if (t % 500 == 0)
        {
            for (int i = 0; i < max_points; i++)
            {
                fout << point.prim[0][i] << " " << point.prim[1][i] << " " << point.prim[2][i] << " " << point.prim[3][i] << " " << point.prim[4][i] << endl;
            }
        }
    }
    cudaDeviceSynchronize();
    fout.close();
    fout_1.close();
    cudaMemcpy(&wall_points_index, wall_points_index_d, wall_size, cudaMemcpyDeviceToHost);
    cudaMemcpy(&outer_points_index, outer_points_index_d, outer_size, cudaMemcpyDeviceToHost);
    cudaMemcpy(&interior_points_index, interior_points_index_d, interior_size, cudaMemcpyDeviceToHost);
    cudaMemcpy(&supersonic_inlet_points_index, supersonic_inlet_points_index_d, supersonic_inlet_points_size, cudaMemcpyDeviceToHost);
    cudaMemcpy(&supersonic_outlet_points_index, supersonic_outlet_points_index_d, supersonic_outlet_points_size, cudaMemcpyDeviceToHost);

    cudaFree(wall_points_index_d);
    cudaFree(outer_points_index_d);
    cudaFree(interior_points_index_d);
    cudaFree(supersonic_inlet_points_index_d);
    cudaFree(supersonic_outlet_points_index_d);
    cudaFree(sum_res_sqr_d);
}

void fpi_solver_multi_nccl(splitPoints *splitPoint_d, int localRank, transferPoints **sendBuffer_d, transferPoints **receiveBuffer_d, int nRanks, int myRank, int *sendPoints, int *receivePoints, ncclComm_t comm, cudaStream_t stream, transferPoints **sendPointer, transferPoints **receivePointer, int *globalToLocalIndex_temp, int **globalToGhostIndex_d, int **globalToGhostIndexSendPointer, int **globalToGhostIndexReceivePointer, int *partVector_d)
{

    dim3 threads(threads_per_block, 1, 1);
    dim3 grid(ceil((numberOfPointsPerDevice / threads.x) + 1), 1, 1);

    int *wallPointsLocalIndex_d;
    unsigned long long wall_size = wallPointsLocal * sizeof(int);
    cudaMalloc(&wallPointsLocalIndex_d, wall_size);

    int *outerPointsLocalIndex_d;
    unsigned long long outer_size = outerPointsLocal * sizeof(int);
    cudaMalloc(&outerPointsLocalIndex_d, outer_size);

    int *interiorPointsLocalIndex_d;
    unsigned long long interior_size = interiorPointsLocal * sizeof(int);
    cudaMalloc(&interiorPointsLocalIndex_d, interior_size);

    int *symmetryPointsLocalIndex_d;
    unsigned long long symmetry_size = symmetryPointsLocal * sizeof(int);
    cudaMalloc(&symmetryPointsLocalIndex_d, symmetry_size);

    int *supersonicInletPointsLocalIndex_d;
    unsigned long long supersonic_inlet_size = supersonicInletPointsLocal * sizeof(int);
    cudaMalloc(&supersonicInletPointsLocalIndex_d, supersonic_inlet_size);

    int *supersonicOutletPointsLocalIndex_d;
    unsigned long long supersonic_outlet_size = supersonicOutletPointsLocal * sizeof(int);
    cudaMalloc(&supersonicOutletPointsLocalIndex_d, supersonic_outlet_size);

    double *sum_res_sqr_d = NULL;
    cudaMalloc(&sum_res_sqr_d, sizeof(double) * (local_points));
    cudaMemset(sum_res_sqr_d, 0, sizeof(double) * (local_points));

    double *sum_res_sqr_final_d = NULL;
    cudaMalloc(&sum_res_sqr_final_d, sizeof(double));
    cudaMemset(sum_res_sqr_final_d, 0, sizeof(double));

    cudaMemcpyAsync(wallPointsLocalIndex_d, wallPointsLocalIndex, wall_size, cudaMemcpyHostToDevice);
    cudaMemcpyAsync(outerPointsLocalIndex_d, outerPointsLocalIndex, outer_size, cudaMemcpyHostToDevice);
    cudaMemcpyAsync(interiorPointsLocalIndex_d, interiorPointsLocalIndex, interior_size, cudaMemcpyHostToDevice);
    cudaMemcpyAsync(symmetryPointsLocalIndex_d, symmetryPointsLocalIndex, symmetry_size, cudaMemcpyHostToDevice);
    cudaMemcpyAsync(supersonicInletPointsLocalIndex_d, supersonicInletPointsLocalIndex, supersonic_inlet_size, cudaMemcpyHostToDevice);
    cudaMemcpyAsync(supersonicOutletPointsLocalIndex_d, supersonicOutletPointsLocalIndex, supersonic_outlet_size, cudaMemcpyHostToDevice);

   
    ncclGroupStart();
    for (int i = 0; i < nRanks; ++i)
    {
        if (i != myRank)
        {
            NCCLCHECK(ncclSend(globalToGhostIndexSendPointer[i], max_points * sizeof(int), ncclChar, i, comm, stream));
            NCCLCHECK(ncclRecv(globalToGhostIndexReceivePointer[i], max_points * sizeof(int), ncclChar, i, comm, stream));
        }
    }
    ncclGroupEnd();

    int *globalToLocalIndex_d;
    CUDACHECK(cudaMalloc(&globalToLocalIndex_d, max_points * sizeof(int)));
    
    ncclGroupStart();
    NCCLCHECK(ncclAllReduce(globalToLocalIndex_temp, globalToLocalIndex_d, max_points, ncclInt, ncclSum, comm, stream));
    ncclGroupEnd();

    ncclGroupStart();
    for (int i = 0; i < nRanks; ++i)
    {
        if (i != myRank)
        {
            NCCLCHECK(ncclSend(sendPointer[i], sendPoints[i] * sizeof(transferPoints), ncclChar, i, comm, stream));
            NCCLCHECK(ncclRecv(receivePointer[i], receivePoints[i] * sizeof(transferPoints), ncclChar, i, comm, stream));
        }
    }
    ncclGroupEnd();
    // printf("cudaProfilerStart\n");
    cudaProfilerStart();
    nvtxRangePushA("Generate Split Stencil ");
    generate_split_stencils_interior_multi_nccl<<<grid, threads>>>(splitPoint_d, interiorPointsLocalIndex_d, interiorPointsLocal, partVector_d, globalToLocalIndex_d, globalToGhostIndex_d, receiveBuffer_d);
    generate_split_stencils_wall_multi_nccl<<<grid, threads>>>(myRank, splitPoint_d, wallPointsLocalIndex_d, wallPointsLocal, partVector_d, globalToLocalIndex_d, globalToGhostIndex_d, receiveBuffer_d);
    generate_split_stencils_outer_multi_nccl<<<grid, threads>>>(splitPoint_d, outerPointsLocalIndex_d, outerPointsLocal, partVector_d, globalToLocalIndex_d, globalToGhostIndex_d, receiveBuffer_d);
    nvtxRangePop();
    cudaProfilerStop();
    // printf("cudaProfilerStop\n");
    for (int t = 1; t <= max_iters; ++t)
    {
        cudaProfilerStart();
        nvtxRangePushA("q-variables");
        eval_q_variables_multi_nccl<<<grid, threads>>>(myRank, splitPoint_d, numberOfPointsPerDevice, sendBuffer_d);
        nvtxRangePop();
        cudaProfilerStop();
        ncclGroupStart();
        for (int i = 0; i < nRanks; ++i)
        {
            if (i != myRank)
            {
                NCCLCHECK(ncclSend(sendPointer[i], sendPoints[i] * sizeof(transferPoints), ncclChar, i, comm, stream));
                NCCLCHECK(ncclRecv(receivePointer[i], receivePoints[i] * sizeof(transferPoints), ncclChar, i, comm, stream));
            }
        }
        ncclGroupEnd();
        cudaProfilerStart();
        nvtxRangePushA("q-derivatives");
        eval_q_derivatives_multi_nccl<<<grid, threads>>>(myRank, splitPoint_d, power, numberOfPointsPerDevice, globalToLocalIndex_d, globalToGhostIndex_d, receiveBuffer_d, partVector_d, sendBuffer_d);
        nvtxRangePop();
        cudaProfilerStop();
        ncclGroupStart();
        for (int i = 0; i < nRanks; ++i)
        {
            if (i != myRank)
            {
                NCCLCHECK(ncclSend(sendPointer[i], sendPoints[i] * sizeof(transferPoints), ncclChar, i, comm, stream));
                NCCLCHECK(ncclRecv(receivePointer[i], receivePoints[i] * sizeof(transferPoints), ncclChar, i, comm, stream));
            }
        }
        ncclGroupEnd();
        nvtxRangePushA("second-order-derivatives");
        for (int r = 0; r < inner_iterations; r++)
        {
            cudaProfilerStart();
            q_inner_loop_multi_nccl<<<grid, threads>>>(myRank, splitPoint_d, power, numberOfPointsPerDevice, globalToLocalIndex_d, globalToGhostIndex_d, receiveBuffer_d, partVector_d, sendBuffer_d);
            update_inner_loop_multi_nccl<<<grid, threads>>>(myRank,splitPoint_d,numberOfPointsPerDevice,sendBuffer_d);
            cudaProfilerStop();
            ncclGroupStart();
            for (int i = 0; i < nRanks; ++i)
            {
                if (i != myRank)
                {
                    NCCLCHECK(ncclSend(sendPointer[i], sendPoints[i] * sizeof(transferPoints), ncclChar, i, comm, stream));
                    NCCLCHECK(ncclRecv(receivePointer[i], receivePoints[i] * sizeof(transferPoints), ncclChar, i, comm, stream));
                }
            }
            ncclGroupEnd();
        }
        nvtxRangePop();
        cudaProfilerStart();
        timestep_delt_multi_nccl<<<grid, threads>>>(myRank, splitPoint_d, CFL, numberOfPointsPerDevice, globalToLocalIndex_d, globalToGhostIndex_d, receiveBuffer_d, partVector_d, sendBuffer_d);

        nvtxRangePushA("flux-residual");
        wall_dGx_pos_multi_nccl<<<grid, threads>>>(myRank, splitPoint_d, power, VL_CONST, pi, wallPointsLocal, wallPointsLocalIndex_d, globalToLocalIndex_d, globalToGhostIndex_d, receiveBuffer_d, partVector_d);
        wall_dGx_neg_multi_nccl<<<grid, threads>>>(myRank, splitPoint_d, power, VL_CONST, pi, wallPointsLocal, wallPointsLocalIndex_d, globalToLocalIndex_d, globalToGhostIndex_d, receiveBuffer_d, partVector_d);
        wall_dGy_pos_multi_nccl<<<grid, threads>>>(splitPoint_d, power, VL_CONST, pi, wallPointsLocal, wallPointsLocalIndex_d, globalToLocalIndex_d, globalToGhostIndex_d, receiveBuffer_d, partVector_d);
        wall_dGy_neg_multi_nccl<<<grid, threads>>>(myRank,splitPoint_d, power, VL_CONST, pi, wallPointsLocal, wallPointsLocalIndex_d, globalToLocalIndex_d, globalToGhostIndex_d, receiveBuffer_d, partVector_d);
        wall_dGz_neg_multi_nccl<<<grid, threads>>>(splitPoint_d, power, VL_CONST, pi, wallPointsLocal, wallPointsLocalIndex_d, globalToLocalIndex_d, globalToGhostIndex_d, receiveBuffer_d, partVector_d);

        interior_dGx_pos_multi_nccl<<<grid, threads>>>(myRank, splitPoint_d, power, VL_CONST, pi, interiorPointsLocal, interiorPointsLocalIndex_d, globalToLocalIndex_d, globalToGhostIndex_d, receiveBuffer_d, partVector_d);
        interior_dGx_neg_multi_nccl<<<grid, threads>>>(splitPoint_d, power, VL_CONST, pi, interiorPointsLocal, interiorPointsLocalIndex_d, globalToLocalIndex_d, globalToGhostIndex_d, receiveBuffer_d, partVector_d);
        interior_dGy_pos_multi_nccl<<<grid, threads>>>(splitPoint_d, power, VL_CONST, pi, interiorPointsLocal, interiorPointsLocalIndex_d, globalToLocalIndex_d, globalToGhostIndex_d, receiveBuffer_d, partVector_d);
        interior_dGy_neg_multi_nccl<<<grid, threads>>>(splitPoint_d, power, VL_CONST, pi, interiorPointsLocal, interiorPointsLocalIndex_d, globalToLocalIndex_d, globalToGhostIndex_d, receiveBuffer_d, partVector_d);
        interior_dGz_pos_multi_nccl<<<grid, threads>>>(splitPoint_d, power, VL_CONST, pi, interiorPointsLocal, interiorPointsLocalIndex_d, globalToLocalIndex_d, globalToGhostIndex_d, receiveBuffer_d, partVector_d);
        interior_dGz_neg_multi_nccl<<<grid, threads>>>(splitPoint_d, power, VL_CONST, pi, interiorPointsLocal, interiorPointsLocalIndex_d, globalToLocalIndex_d, globalToGhostIndex_d, receiveBuffer_d, partVector_d);

        outer_dGx_pos_multi_nccl<<<grid, threads>>>(myRank, splitPoint_d, power, VL_CONST, pi, outerPointsLocal, outerPointsLocalIndex_d, globalToLocalIndex_d, globalToGhostIndex_d, receiveBuffer_d, partVector_d);
        outer_dGx_neg_multi_nccl<<<grid, threads>>>(splitPoint_d, power, VL_CONST, pi, outerPointsLocal, outerPointsLocalIndex_d, globalToLocalIndex_d, globalToGhostIndex_d, receiveBuffer_d, partVector_d);
        outer_dGy_pos_multi_nccl<<<grid, threads>>>(splitPoint_d, power, VL_CONST, pi, outerPointsLocal, outerPointsLocalIndex_d, globalToLocalIndex_d, globalToGhostIndex_d, receiveBuffer_d, partVector_d);
        outer_dGy_neg_multi_nccl<<<grid, threads>>>(myRank,splitPoint_d, power, VL_CONST, pi, outerPointsLocal, outerPointsLocalIndex_d, globalToLocalIndex_d, globalToGhostIndex_d, receiveBuffer_d, partVector_d);
        outer_dGz_pos_multi_nccl<<<grid, threads>>>(splitPoint_d, power, VL_CONST, pi, outerPointsLocal, outerPointsLocalIndex_d, globalToLocalIndex_d, globalToGhostIndex_d, receiveBuffer_d, partVector_d);
        nvtxRangePop();
        nvtxRangePushA("state-update");
        state_update_wall_multi_nccl<<<grid, threads>>>(myRank, splitPoint_d, wallPointsLocal, wallPointsLocalIndex_d, sum_res_sqr_d);
        state_update_outer_multi_nccl<<<grid, threads>>>(myRank,splitPoint_d, outerPointsLocal, outerPointsLocalIndex_d, u1_inf, u2_inf, u3_inf, rho_inf, pi, pr_inf);
        state_update_interior_multi_nccl<<<grid, threads>>>(splitPoint_d, interiorPointsLocal, interiorPointsLocalIndex_d, sum_res_sqr_d);
        state_update_symmetric_multi_nccl<<<grid,threads>>>(splitPoint_d, power, VL_CONST, pi, symmetryPointsLocal, symmetryPointsLocalIndex_d, globalToLocalIndex_d, globalToGhostIndex_d, receiveBuffer_d, partVector_d);
        nvtxRangePop();
        cudaDeviceSynchronize();
        sum_res_sqr = thrust::reduce(thrust::cuda::par.on(stream), sum_res_sqr_d, sum_res_sqr_d + local_points, (double)0.0, thrust::plus<double>());
        cudaProfilerStop();
        // MPI_Barrier(MPI_COMM_WORLD);
        if (myRank == 0)
        {
            MPICHECK(MPI_Reduce(MPI_IN_PLACE, &sum_res_sqr, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD));
        }
        else
        {
            MPICHECK(MPI_Reduce(&sum_res_sqr, &sum_res_sqr, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD));
        }
        // CUDACHECK(cudaMemcpy(splitPoint, splitPoint_d, numberOfPointsPerDevice * sizeof(splitPoints), cudaMemcpyDeviceToHost));
        if (myRank == 0)
        {
            res_new = sqrt(sum_res_sqr) / max_points;
            if (t <= 2)
            {
                res_old = res_new;
                residue = 0.00;
            }
            else
            {
                residue = log10(res_new / res_old);
            }
            cout << t << " " << res_new << " " << residue << " " << sum_res_sqr << endl;
            // fstream fout;
            // fout.open("error.dat",ios::out);
            // for(int r=0;r<max_points;r++){
            //     fout<<splitPoint[r].prim[0]<<" "<<splitPoint[r].prim[1]<<" "<<splitPoint[r].prim[2]<<" "<<splitPoint[r].prim[3]<<" "<<splitPoint[r].prim[4]<<endl;
            // }
            // fout.close();
        }
    }
    // cudaMemcpy(&wallPointsLocalIndex, wallPointsLocalIndex_d, wall_size, cudaMemcpyDeviceToHost);
    // cudaMemcpy(&outerPointsLocalIndex, outerPointsLocalIndex_d, outer_size, cudaMemcpyDeviceToHost);
    // cudaMemcpy(&interiorPointsLocalIndex, interiorPointsLocalIndex_d, interior_size, cudaMemcpyDeviceToHost);
    // cudaMemcpy(&supersonicInletPointsLocalIndex, supersonicInletPointsLocalIndex_d, supersonic_inlet_size, cudaMemcpyDeviceToHost);
    // cudaMemcpy(&supersonicOutletPointsLocalIndex, supersonicOutletPointsLocalIndex_d, supersonic_outlet_size, cudaMemcpyDeviceToHost);

    cudaFree(wallPointsLocalIndex_d);
    cudaFree(outerPointsLocalIndex_d);
    cudaFree(interiorPointsLocalIndex_d);
    cudaFree(symmetryPointsLocalIndex_d);
    cudaFree(supersonicInletPointsLocalIndex_d);
    cudaFree(supersonicOutletPointsLocalIndex_d);
    cudaFree(sum_res_sqr_d);
}