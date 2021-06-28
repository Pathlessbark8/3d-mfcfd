#pragma once
#include "data_structure_mod.h"

void eval_q_variables()

{

	int k;
	double rho, u1, u2, u3, pr, beta;
	double two_times_beta;
	//
	//
	for (k = 0; k < max_points; k++)
	//
	{
		rho = point.prim[0][k];
		u1 = point.prim[1][k];
		u2 = point.prim[2][k];
		u3 = point.prim[3][k];
		pr = point.prim[4][k];
		//
		beta = 0.50 * rho / pr;
		//
		point.q[0][k] = log(rho) + (log(beta) * 2.50) - beta * (u1 * u1 + u2 * u2 + u3 * u3);
		//
		two_times_beta = 2.00 * beta;
		//
		point.q[1][k] = two_times_beta * u1;
		point.q[2][k] = two_times_beta * u2;
		point.q[3][k] = two_times_beta * u3;
		point.q[4][k] = -two_times_beta;
		//
	}
	//
}
//
//
