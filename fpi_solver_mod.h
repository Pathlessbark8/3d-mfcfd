#pragma once

#include "data_structure_mod.h"
#include "compute_conserved_vector_mod.h"
#include "q_variables_mod.h"
#include "q_derivatives_mod.h"
#include "timestep_delt_mod.h"
#include "flux_residual_mod.h"
#include "state_update_mod.h"
#include <cuda_runtime.h>
#include <thrust/reduce.h>
#include <thrust/execution_policy.h>

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

void fpi_solver_cuda( points *point_d,cudaStream_t stream)
{

    dim3 threads(threads_per_block, 1, 1);
    dim3 grid(ceil((max_points / threads.x) + 1), 1, 1);

    fstream fout,fout_1;
    fout.open("output_40k.dat",ios::out);
    fout << setprecision(13);
    fout_1.open("residue.dat",ios::out);
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
    cudaMemcpy(interior_points_index_d, &interior_points_index, wall_size, cudaMemcpyHostToDevice);

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
        cout << t << " " << res_new << " " << residue<<" "<<sum_res_sqr << endl;
        fout_1<< t << " " << res_new << " " << residue << endl;
        if(t%500==0)
        { 
            for(int i=0;i<max_points;i++)
            {
                fout<<point.prim[0][i]<<" "<<point.prim[1][i]<<" "<<point.prim[2][i]<<" "<<point.prim[3][i]<<" "<<point.prim[4][i]<<endl;
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

void fpi_solver_multi_nccl( splitPoints **splitPoint_d,cudaStream_t *stream)
{

    dim3 threads(threads_per_block, 1, 1);
    dim3 grid(ceil((max_points / threads.x) + 1), 1, 1);
    for (int t = 1; t <= max_iters; ++t)
    {
        for(int i=0;i<numDevices;++i){
            {
            cudaSetDevice(i);
            cudaStreamCreate(&stream[i]);
            eval_q_variables_multi_nccl<<<grid,threads,0,stream[i]>>>(splitPoint_d[i],numberOfPointsPerDevice[i]);
            }
        }
    }
}