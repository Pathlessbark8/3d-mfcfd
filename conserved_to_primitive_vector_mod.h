#pragma once
//
//
#include "data_structure_mod.h"
//
//

//
//
void conserved_to_primitive(int k, double *U)
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
//
//
