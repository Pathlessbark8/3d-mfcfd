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

void fpi_solver_cuda(points *point_d, cudaStream_t stream)
{

    dim3 threads(threads_per_block, 1, 1);
    dim3 grid(ceil((max_points / threads.x) + 1), 1, 1);

    compute_conserved_vector_cuda<<<grid, threads>>>(*point_d);
    // eval_q_variables_cuda<<<grid, threads>>>(point_d);
    // eval_q_derivatives_cuda<<<grid, threads>>>(point_d, power);
    // timestep_delt_cuda<<<grid, threads>>>(point_d, CFL);

    // wall_dGx_pos_cuda<<<grid, threads>>>(point_d, power,VL_CONST);
    // wall_dGx_neg_cuda<<<grid, threads>>>(point_d, power,VL_CONST);
    // wall_dGy_pos_cuda<<<grid, threads>>>(point_d, power,VL_CONST);
    // wall_dGy_neg_cuda<<<grid, threads>>>(point_d, power,VL_CONST);
    // wall_dGz_neg_cuda<<<grid, threads>>>(point_d, power,VL_CONST);
    // //
    // outer_dGx_pos_cuda<<<grid, threads>>>(point_d, power,VL_CONST);
    // outer_dGx_neg_cuda<<<grid, threads>>>(point_d, power,VL_CONST);
    // outer_dGy_pos_cuda<<<grid, threads>>>(point_d, power,VL_CONST);
    // outer_dGy_neg_cuda<<<grid, threads>>>(point_d, power,VL_CONST);
    // outer_dGz_pos_cuda<<<grid, threads>>>(point_d, power,VL_CONST);
    // //
    // interior_dGx_pos_cuda<<<grid, threads>>>(point_d, power,VL_CONST,pi);
    // interior_dGx_neg_cuda<<<grid, threads>>>(point_d, power,VL_CONST,pi);
    // interior_dGy_pos_cuda<<<grid, threads>>>(point_d, power,VL_CONST,pi);
    // interior_dGy_neg_cuda<<<grid, threads>>>(point_d, power,VL_CONST,pi);
    // interior_dGz_pos_cuda<<<grid, threads>>>(point_d, power,VL_CONST,pi);
    // interior_dGz_neg_cuda<<<grid, threads>>>(point_d, power,VL_CONST,pi);
    //
    cudaDeviceSynchronize();
}