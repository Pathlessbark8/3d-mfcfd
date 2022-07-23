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
		point.U[k][0] = point.prim[k][0];
		temp1 = point.prim[k][0] * point.prim[k][1];
		temp2 = point.prim[k][0] * point.prim[k][2];
		temp3 = point.prim[k][0] * point.prim[k][3];
		//
		point.U[k][1] = temp1 * point.tan1[k][0] + temp2 * point.tan1[k][1] + temp3 * point.tan1[k][2];
		point.U[k][2] = temp1 * point.tan2[k][0] + temp2 * point.tan2[k][1] + temp3 * point.tan2[k][2];
		point.U[k][3] = temp1 * point.nor[k][0] + temp2 * point.nor[k][1] + temp3 * point.nor[k][2];
		//
		temp1 = point.prim[k][1] * point.prim[k][1] + point.prim[k][2] * point.prim[k][2] + point.prim[k][3] * point.prim[k][3];
		point.U[k][4] = 2.50 * point.prim[k][4] + 0.50 * point.U[k][0] * temp1;
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
	point.U[k][0] = point.prim[k][0];
	temp1 = point.prim[k][0] * point.prim[k][1];
	temp2 = point.prim[k][0] * point.prim[k][2];
	temp3 = point.prim[k][0] * point.prim[k][3];
	//
	point.U[k][1] = temp1 * point.tan1[k][0] + temp2 * point.tan1[k][1] + temp3 * point.tan1[k][2];
	point.U[k][2] = temp1 * point.tan2[k][0] + temp2 * point.tan2[k][1] + temp3 * point.tan2[k][2];
	point.U[k][3] = temp1 * point.nor[k][0] + temp2 * point.nor[k][1] + temp3 * point.nor[k][2];
	//
	temp1 = point.prim[k][1] * point.prim[k][1] + point.prim[k][2] * point.prim[k][2] + point.prim[k][3] * point.prim[k][3];
	point.U[k][4] = 2.50 * point.prim[k][4] + 0.50 * point.U[k][0] * temp1;
	//
}
