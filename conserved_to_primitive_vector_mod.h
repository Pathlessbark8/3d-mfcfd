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
//
#include "data_structure_mod.h"
//
//

//
//
void conserved_to_primitive(int k, double *U)
{
	//
	double temp;
	double U2_rot, U3_rot, U4_rot;
	//
	//
	point.prim[0][k] = U[0];
	temp = 1.00 / U[0];
	//
	//
	U2_rot = U[1];
	U3_rot = U[2];
	U4_rot = U[3];
	//
	//
	U[1] = point.tan1[0][k] * U2_rot + point.tan2[0][k] * U3_rot + point.nor[0][k] * U4_rot;
	U[2] = point.tan1[1][k] * U2_rot + point.tan2[1][k] * U3_rot + point.nor[1][k] * U4_rot;
	U[3] = point.tan1[2][k] * U2_rot + point.tan2[2][k] * U3_rot + point.nor[2][k] * U4_rot;
	//
	//
	point.prim[1][k] = temp * U[1];
	point.prim[2][k] = temp * U[2];
	point.prim[3][k] = temp * U[3];
	//
	temp = point.prim[1][k] * point.prim[1][k] + point.prim[2][k] * point.prim[2][k] + point.prim[3][k] * point.prim[3][k];

	point.prim[4][k] = 0.40 * (U[4] - 0.50 * point.prim[0][k] * temp);
	//
	//
}


__device__ void conserved_to_primitive_cuda(points &point,int k, double *U)
//
{
	//
	double temp;
	double U2_rot, U3_rot, U4_rot;
	//
	//
	point.prim[0][k] = U[0];
	temp = 1.00 / U[0];
	//
	//
	U2_rot = U[1];
	U3_rot = U[2];
	U4_rot = U[3];
	//
	//
	U[1] = point.tan1[0][k] * U2_rot + point.tan2[0][k] * U3_rot + point.nor[0][k] * U4_rot;
	U[2] = point.tan1[1][k] * U2_rot + point.tan2[1][k] * U3_rot + point.nor[1][k] * U4_rot;
	U[3] = point.tan1[2][k] * U2_rot + point.tan2[2][k] * U3_rot + point.nor[2][k] * U4_rot;
	//
	//
	point.prim[1][k] = temp * U[1];
	point.prim[2][k] = temp * U[2];
	point.prim[3][k] = temp * U[3];
	//
	temp = point.prim[1][k] * point.prim[1][k] + point.prim[2][k] * point.prim[2][k] + point.prim[3][k] * point.prim[3][k];

	point.prim[4][k] = 0.40 * (U[4] - 0.50 * point.prim[0][k] * temp);
	//
	//
}
//

__device__ void conserved_to_primitive_multi_nccl(splitPoints *splitPoint,int k, double *U)
//
{
	//
	double temp;
	double U2_rot, U3_rot, U4_rot;
	//
	//
	splitPoint[k].prim[0] = U[0];
	temp = 1.00 / U[0];
	//
	//
	U2_rot = U[1];
	U3_rot = U[2];
	U4_rot = U[3];
	//
	//
	U[1] = splitPoint[k].tan1[0] * U2_rot + splitPoint[k].tan2[0] * U3_rot + splitPoint[k].nor[0] * U4_rot;
	U[2] = splitPoint[k].tan1[1] * U2_rot + splitPoint[k].tan2[1] * U3_rot + splitPoint[k].nor[1] * U4_rot;
	U[3] = splitPoint[k].tan1[2] * U2_rot + splitPoint[k].tan2[2] * U3_rot + splitPoint[k].nor[2] * U4_rot;
	//
	//
	splitPoint[k].prim[1] = temp * U[1];
	splitPoint[k].prim[2] = temp * U[2];
	splitPoint[k].prim[3] = temp * U[3];
	//
	temp = splitPoint[k].prim[1] * splitPoint[k].prim[1] + splitPoint[k].prim[2] * splitPoint[k].prim[2] + splitPoint[k].prim[3] * splitPoint[k].prim[3];

	splitPoint[k].prim[4] = 0.40 * (U[4] - 0.50 * splitPoint[k].prim[0] * temp);
	//
	//
}