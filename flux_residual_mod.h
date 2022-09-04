//
/*  
	MFCFD is a 3D Computational Fluid Dynamics Solver based off q-LSKUM
    Copyright (C) 2022 Dhruv Saxena
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
//
//
#include "data_structure_mod.h"
//
#include "interior_flux_dGxpos_mod.h"
#include "interior_flux_dGxneg_mod.h"
#include "interior_flux_dGypos_mod.h"
#include "interior_flux_dGyneg_mod.h"
#include "interior_flux_dGzpos_mod.h"
#include "interior_flux_dGzneg_mod.h"
//
#include "wall_flux_dGxpos_mod.h"
#include "wall_flux_dGxneg_mod.h"
#include "wall_flux_dGypos_mod.h"
#include "wall_flux_dGyneg_mod.h"
#include "wall_flux_dGzneg_mod.h"
//
#include "outer_flux_dGxpos_mod.h"
#include "outer_flux_dGxneg_mod.h"
#include "outer_flux_dGypos_mod.h"
#include "outer_flux_dGyneg_mod.h"
#include "outer_flux_dGzpos_mod.h"
//
//
//
//
void eval_flux_residual()
//
//
{
    //
    int i, k;
    double Gxp[5] = {0}, Gyp[5] = {0}, Gzp[5] = {0};
    double Gxn[5] = {0}, Gyn[5] = {0}, Gzn[5] = {0};
    //

    for (i = 0; i < wall_points; i++)
    {
        //
        k = wall_points_index[i];

        wall_dGx_pos(Gxp, k);
        wall_dGx_neg(Gxn, k);
        wall_dGy_pos(Gyp, k);
        wall_dGy_neg(Gyn, k);
        wall_dGz_neg(Gzn, k);
        //
        for (int r = 0; r < 5; r++)
        {
            point.flux_res[r][k] = 2.0*point.delt[k]*(Gxp[r] + Gxn[r] + Gyp[r] + Gyn[r] + Gzn[r]);
        }
    }

    for (i = 0; i < outer_points; i++)
    {
        //
        k = outer_points_index[i];
        //
        outer_dGx_pos(Gxp, k);
        outer_dGx_neg(Gxn, k);
        outer_dGy_pos(Gyp, k);
        outer_dGy_neg(Gyn, k);
        outer_dGz_pos(Gzp, k);
        //
        for (int r = 0; r < 5; r++)
        {
            point.flux_res[r][k] = point.delt[k]*(Gxp[r] + Gxn[r] + Gyp[r] + Gyn[r] + Gzp[r]);
        }
        //
    }
    //
    for (i = 0; i < interior_points; i++)
    {
        //
        k = interior_points_index[i];
        //
        interior_dGx_pos(Gxp, k);
        interior_dGx_neg(Gxn, k);
        interior_dGy_pos(Gyp, k);
        interior_dGy_neg(Gyn, k);
        interior_dGz_pos(Gzp, k);
        interior_dGz_neg(Gzn, k);
        //
        //
        for (int r = 0; r < 5; r++)
        {
            point.flux_res[r][k] = point.delt[k]*(Gxp[r] + Gxn[r] + Gyp[r] + Gyn[r] + Gzn[r] + Gzp[r]);
        }
    }
    //
    //
}

//
//
