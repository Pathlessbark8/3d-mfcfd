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
		rho = point.prim[k][0];
		u1 = point.prim[k][1];
		u2 = point.prim[k][2];
		u3 = point.prim[k][3];
		pr = point.prim[k][4];
		//
		beta = 0.50 * rho / pr;
		//
		point.q[k][0] = log(rho) + (log(beta) * 2.50) - beta * (u1 * u1 + u2 * u2 + u3 * u3);
		//
		two_times_beta = 2.00 * beta;
		//
		point.q[k][1] = two_times_beta * u1;
		point.q[k][2] = two_times_beta * u2;
		point.q[k][3] = two_times_beta * u3;
		point.q[k][4] = -two_times_beta;
		//
	}
	//
}
//
__global__ void eval_q_variables_cuda(points &point)

{

	double rho, u1, u2, u3, pr, beta;
	double two_times_beta;
	//
	int bx = blockIdx.x;
    int tx = threadIdx.x;
    int k = bx * blockDim.x + tx;

	if (k < 0 || k >= max_points){
        return;
    }
	//
		rho = point.prim[k][0];
		u1 = point.prim[k][1];
		u2 = point.prim[k][2];
		u3 = point.prim[k][3];
		pr = point.prim[k][4];
		//
		beta = 0.50 * rho / pr;
		//
		point.q[k][0] = log(rho) + (log(beta) * 2.50) - beta * (u1 * u1 + u2 * u2 + u3 * u3);
		//
		two_times_beta = 2.00 * beta;
		//
		point.q[k][1] = two_times_beta * u1;
		point.q[k][2] = two_times_beta * u2;
		point.q[k][3] = two_times_beta * u3;
		point.q[k][4] = -two_times_beta;
		//
	//
}
//
