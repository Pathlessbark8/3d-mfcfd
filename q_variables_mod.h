#pragma once
#include "data_structure_mod.h"
// #include "cuPrintf.cu"

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
	//
}
//

__global__ void eval_q_variables_multi_nccl(splitPoints splitPoint[],int max_points_on_device)

{

	double rho, u1, u2, u3, pr, beta;
	double two_times_beta;
	//
	int bx = blockIdx.x;
    int tx = threadIdx.x;
    int k = bx * blockDim.x + tx;

	if(k>=max_points_on_device){
		return;
	}
	// splitPoint[0].x=2987;
	if (k < 0 ){
        return;
    }
	// cuPrintf("%d\n",k);
	// 	rho = splitPoint[k].prim[0];
	// 	u1 = splitPoint[k].prim[1];
	// 	u2 = splitPoint[k].prim[2];
	// 	u3 = splitPoint[k].prim[3];
	// 	pr = splitPoint[k].prim[4];
	// 	//
	// 	beta = 0.50 * rho / pr;
	// 	//
	// 	splitPoint[k].q[0] = log(rho) + (log(beta) * 2.50) - beta * (u1 * u1 + u2 * u2 + u3 * u3);
	// 	//
	// 	two_times_beta = 2.00 * beta;
	// 	//
	// 	splitPoint[k].q[1] = two_times_beta * u1;
	// 	splitPoint[k].q[2] = two_times_beta * u2;
	// 	splitPoint[k].q[3] = two_times_beta * u3;
	// 	splitPoint[k].q[4] = -two_times_beta;
		splitPoint[k].x=234;
		//
	//
}
