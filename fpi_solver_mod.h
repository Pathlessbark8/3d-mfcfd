//
#pragma once
//
//	First written on 06.02.2020
//
//
#include "data_structure_mod.h"
#include "q_variables_mod.h"
#include "q_derivatives_mod.h"
#include "timestep_delt_mod.h"
#include "flux_residual_mod.h"
#include "state_update_mod.h"
//
//
//
//
void fpi_solver(int t)
//
{
    //
    int i, k;
    //
    //
    eval_q_variables();
    eval_q_derivatives();
    //
    timestep_delt();
    eval_flux_residual();
// for(int i=0;i<wall_points;i++)
//     {
//         int k= wall_points_index[i];
//         for(int r=0;r<5;r++)
//         {
//             cout<<point.flux_res[r][k]<<" ";
//         }
//         cout<<endl;
//     }
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
//
//
