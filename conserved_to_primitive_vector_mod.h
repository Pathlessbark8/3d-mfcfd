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
	point.prim[k][0] = U[0];
	temp = 1.00 / U[0];
	//
	//
	U2_rot = U[1];
	U3_rot = U[2];
	U4_rot = U[3];
	//
	//
	U[1] = point.tan1[k][0] * U2_rot + point.tan2[k][0] * U3_rot + point.nor[k][0] * U4_rot;
	U[2] = point.tan1[k][1] * U2_rot + point.tan2[k][1] * U3_rot + point.nor[k][1] * U4_rot;
	U[3] = point.tan1[k][2] * U2_rot + point.tan2[k][2] * U3_rot + point.nor[k][2] * U4_rot;
	//
	//
	point.prim[k][1] = temp * U[1];
	point.prim[k][2] = temp * U[2];
	point.prim[k][3] = temp * U[3];
	//
	temp = point.prim[k][1] * point.prim[k][1] + point.prim[k][2] * point.prim[k][2] + point.prim[k][3] * point.prim[k][3];

	point.prim[k][4] = 0.40 * (U[4] - 0.50 * point.prim[k][0] * temp);
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
	point.prim[k][0] = U[0];
	temp = 1.00 / U[0];
	//
	//
	U2_rot = U[1];
	U3_rot = U[2];
	U4_rot = U[3];
	//
	//
	U[1] = point.tan1[k][0] * U2_rot + point.tan2[k][0] * U3_rot + point.nor[k][0] * U4_rot;
	U[2] = point.tan1[k][1] * U2_rot + point.tan2[k][1] * U3_rot + point.nor[k][1] * U4_rot;
	U[3] = point.tan1[k][2] * U2_rot + point.tan2[k][2] * U3_rot + point.nor[k][2] * U4_rot;
	//
	//
	point.prim[k][1] = temp * U[1];
	point.prim[k][2] = temp * U[2];
	point.prim[k][3] = temp * U[3];
	//
	temp = point.prim[k][1] * point.prim[k][1] + point.prim[k][2] * point.prim[k][2] + point.prim[k][3] * point.prim[k][3];

	point.prim[k][4] = 0.40 * (U[4] - 0.50 * point.prim[k][0] * temp);
	//
	//
}
//
