#pragma once

#include "data_structure_mod.h"
#include "compute_conserved_vector_mod.h"
#include "q_variables_mod.h"
#include "q_derivatives_mod.h"
#include "timestep_delt_mod.h"
#include "flux_residual_mod.h"
#include <cuda_runtime.h>

void fpi_solver()
{
    compute_conserved_vector();
    eval_q_variables();
    eval_q_derivatives();
    timestep_delt();
    eval_flux_residual();
}

void fpi_solver_cuda(points *point_d)
{

    dim3 threads(threads_per_block, 1, 1);
    dim3 grid(ceil((max_points / threads.x) + 1), 1, 1);

    compute_conserved_vector_cuda<<<grid, threads>>>(*point_d);
    eval_q_variables_cuda<<<grid, threads>>>(*point_d);
  
    eval_q_derivatives_cuda<<<grid, threads>>>(*point_d, power,inner_iterations);
    timestep_delt_cuda<<<grid, threads>>>(*point_d, CFL);

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

    wall_dGx_pos_cuda<<<grid, threads>>>(*point_d, power,VL_CONST,pi,wall_points,wall_points_index_d);
    wall_dGx_neg_cuda<<<grid, threads>>>(*point_d, power,VL_CONST,pi,wall_points,wall_points_index_d);
    wall_dGy_pos_cuda<<<grid, threads>>>(*point_d, power,VL_CONST,pi,wall_points,wall_points_index_d);
    wall_dGy_neg_cuda<<<grid, threads>>>(*point_d, power,VL_CONST,pi,wall_points,wall_points_index_d);
    wall_dGz_neg_cuda<<<grid, threads>>>(*point_d, power,VL_CONST,pi,wall_points,wall_points_index_d);
    // //
    outer_dGx_pos_cuda<<<grid, threads>>>(*point_d, power,VL_CONST,pi,outer_points,outer_points_index_d);
    outer_dGx_neg_cuda<<<grid, threads>>>(*point_d, power,VL_CONST,pi,outer_points,outer_points_index_d);
    outer_dGy_pos_cuda<<<grid, threads>>>(*point_d, power,VL_CONST,pi,outer_points,outer_points_index_d);
    outer_dGy_neg_cuda<<<grid, threads>>>(*point_d, power,VL_CONST,pi,outer_points,outer_points_index_d);
    outer_dGz_pos_cuda<<<grid, threads>>>(*point_d, power,VL_CONST,pi,outer_points,outer_points_index_d);
    //
    interior_dGx_pos_cuda<<<grid, threads>>>(*point_d, power,VL_CONST,pi,interior_points,interior_points_index_d);
    interior_dGx_neg_cuda<<<grid, threads>>>(*point_d, power,VL_CONST,pi,interior_points,interior_points_index_d);
    interior_dGy_pos_cuda<<<grid, threads>>>(*point_d, power,VL_CONST,pi,interior_points,interior_points_index_d);
    interior_dGy_neg_cuda<<<grid, threads>>>(*point_d, power,VL_CONST,pi,interior_points,interior_points_index_d);
    interior_dGz_pos_cuda<<<grid, threads>>>(*point_d, power,VL_CONST,pi,interior_points,interior_points_index_d);
    interior_dGz_neg_cuda<<<grid, threads>>>(*point_d, power,VL_CONST,pi,interior_points,interior_points_index_d);
    //
    cudaDeviceSynchronize();
    cudaMemcpy(&wall_points_index, wall_points_index_d,wall_size, cudaMemcpyDeviceToHost);
    cudaMemcpy(&outer_points_index, outer_points_index_d,outer_size, cudaMemcpyDeviceToHost);
    cudaMemcpy(&interior_points_index, interior_points_index_d,interior_size, cudaMemcpyDeviceToHost);
    cudaFree(wall_points_index_d);
    cudaFree(outer_points_index_d);
    cudaFree(interior_points_index_d);
}