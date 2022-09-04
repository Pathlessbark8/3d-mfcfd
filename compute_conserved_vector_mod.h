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
//
//	First written on 10.04.2021.
//	Revised on 17.04.2021.
//
//
#include "data_structure_mod.h"
//
//
//
//
void compute_conserved_vector()
//
{
	//
	int k;
	double temp1, temp2, temp3;
	//
	for (k = 0; k < max_points; k++)
	{
		point.U[0][k] = point.prim[0][k];
		temp1 = point.prim[0][k] * point.prim[1][k];
		temp2 = point.prim[0][k] * point.prim[2][k];
		temp3 = point.prim[0][k] * point.prim[3][k];
		//
		point.U[1][k] = temp1 * point.tan1[0][k] + temp2 * point.tan1[1][k] + temp3 * point.tan1[2][k];
		point.U[2][k] = temp1 * point.tan2[0][k] + temp2 * point.tan2[1][k] + temp3 * point.tan2[2][k];
		point.U[3][k] = temp1 * point.nor[0][k] + temp2 * point.nor[1][k] + temp3 * point.nor[2][k];
		//
		temp1 = point.prim[1][k] * point.prim[1][k] + point.prim[2][k] * point.prim[2][k] + point.prim[3][k] * point.prim[3][k];
		point.U[4][k] = 2.50 * point.prim[4][k] + 0.50 * point.U[0][k] * temp1;
	}
	//
}
//
__global__ void compute_conserved_vector_cuda(points &point)
//
{
	//
	double temp1, temp2, temp3;
	//

	int bx = blockIdx.x;
	int tx = threadIdx.x;
	int k = bx * blockDim.x + tx;

	if (k < 0 || k >= max_points)
	{
		return;
	}
	point.U[0][k] = point.prim[0][k];
	temp1 = point.prim[0][k] * point.prim[1][k];
	temp2 = point.prim[0][k] * point.prim[2][k];
	temp3 = point.prim[0][k] * point.prim[3][k];
	//
	point.U[1][k] = temp1 * point.tan1[0][k] + temp2 * point.tan1[1][k] + temp3 * point.tan1[2][k];
	point.U[2][k] = temp1 * point.tan2[0][k] + temp2 * point.tan2[1][k] + temp3 * point.tan2[2][k];
	point.U[3][k] = temp1 * point.nor[0][k] + temp2 * point.nor[1][k] + temp3 * point.nor[2][k];
	//
	temp1 = point.prim[1][k] * point.prim[1][k] + point.prim[2][k] * point.prim[2][k] + point.prim[3][k] * point.prim[3][k];
	point.U[4][k] = 2.50 * point.prim[4][k] + 0.50 * point.U[0][k] * temp1;
	//
}
