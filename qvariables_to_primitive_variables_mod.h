//
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
//
#pragma once
#include "data_structure_mod.h"
//
//
//
void qtilde_to_primitive(double *q, double *prim)
//
{
    //
    double beta, temp, temp1, temp2;
    //
    //
    beta = -q[4] * 0.50;
    //
    temp = 0.50 / beta;
    //
    prim[1] = q[1] * temp;
    prim[2] = q[2] * temp;
    prim[3] = q[3] * temp;
    //
    temp1 = q[0] + beta * (prim[1] * prim[1] + prim[2] * prim[2] + prim[3] * prim[3]);
    temp2 = temp1 - (log(beta) * 2.50);
    //
    prim[0] = exp(temp2);
    prim[4] = prim[0] * temp;
}
//
__device__ void qtilde_to_primitive_cuda(double *q, double *prim)
//
{
    //
    double beta, temp, temp1, temp2;
    //
    //
    beta = -q[4] * 0.50;
    //
    temp = 0.50 / beta;
    //
    prim[1] = q[1] * temp;
    prim[2] = q[2] * temp;
    prim[3] = q[3] * temp;
    //
    temp1 = q[0] + beta * (prim[1] * prim[1] + prim[2] * prim[2] + prim[3] * prim[3]);
    temp2 = temp1 - (log(beta) * 2.50);
    //
    prim[0] = exp(temp2);
    prim[4] = prim[0] * temp;
}
