

//	First written on 10.06.2021.
//
#pragma once
#include "data_structure_mod.h"
#include "split_fluxes_mod.h"
#include "q_variables_mod.h"
#include "q_derivatives_mod.h"
#include "qvariables_to_primitive_variables_mod.h"
#include "limiters_mod.h"

//	This subroutine evaluates the wall flux derivative dGz_neg

void wall_dGz_neg(double *G, int i)
//
//
{

	int j, k;
	double prim[5];
	double x_i, y_i, z_i, x_k, y_k, z_k;
	double tan1[3], tan2[3], nor[3];
	double G_i[5], G_k[5];
	double delx, dely, delz, det;
	double dels, delt, deln;
	//
	double sum_delx_sqr, sum_dely_sqr, sum_delz_sqr;
	double sum_delx_dely, sum_dely_delz, sum_delz_delx;
	double sum_delx_delf[5] = {0}, sum_dely_delf[5] = {0}, sum_delz_delf[5] = {0};
	double dist, weights;
	double temp[5], qtilde[5], phi[5];
	double dels_weights, delt_weights, deln_weights;
	//
	//
	sum_delx_sqr = 0.00;
	sum_dely_sqr = 0.00;
	sum_delz_sqr = 0.00;
	//
	sum_delx_dely = 0.00;
	sum_dely_delz = 0.00;
	sum_delz_delx = 0.00;
	//
	//
	x_i = point.x[i];
	y_i = point.y[i];
	z_i = point.z[i];
	//
	for (int r = 0; r < 3; r++)
	{
		tan1[r] = point.tan1[i][r];
		tan2[r] = point.tan2[i][r];
		nor[r] = point.nor[i][r];
	}
	//
	for (j = 0; j < point.nbhs[i]; j++)
	//
	{
		k = point.conn[i][j];
		//
		x_k = point.x[k];
		y_k = point.y[k];
		z_k = point.z[k];
		//
		delx = x_k - x_i;
		dely = y_k - y_i;
		delz = z_k - z_i;
		//
		dels = delx * tan1[0] + dely * tan1[1] + delz * tan1[2];
		delt = delx * tan2[0] + dely * tan2[1] + delz * tan2[2];
		deln = delx * nor[0] + dely * nor[1] + delz * nor[2];
		//
		dist = sqrt(dels * dels + delt * delt + deln * deln);
		weights = 1.00 / (pow(dist, power));
		//
		dels_weights = dels * weights;
		delt_weights = delt * weights;
		deln_weights = deln * weights;
		//
		sum_delx_sqr = sum_delx_sqr + dels * dels_weights;
		sum_dely_sqr = sum_dely_sqr + delt * delt_weights;
		sum_delz_sqr = sum_delz_sqr + deln * deln_weights;
		//
		sum_delx_dely = sum_delx_dely + dels * delt_weights;
		sum_dely_delz = sum_dely_delz + delt * deln_weights;
		sum_delz_delx = sum_delz_delx + deln * dels_weights;
		//
		for (int r = 0; r < 5; r++)
		{
			temp[r] = delx * point.dq[i][0][r] + dely * point.dq[i][1][r] + delz * point.dq[i][2][r];
			qtilde[r] = point.q[i][r] - 0.50 * temp[r];
		}
		venkat_limiter(qtilde, phi, i);
		for (int r = 0; r < 5; r++)
		{
			qtilde[r] = point.q[i][r] - 0.50 * phi[r] * temp[r];
		}
		qtilde_to_primitive(qtilde, prim);
		flux_Gzn(G_i, tan1, tan2, nor, prim);
		//
		for (int r = 0; r < 5; r++)
		{
			temp[r] = delx * point.dq[k][0][r] + dely * point.dq[k][1][r] + delz * point.dq[k][2][r];
			qtilde[r] = point.q[k][r] - 0.50 * temp[r];
		}
		venkat_limiter(qtilde, phi, k);
		for (int r = 0; r < 5; r++)
		{
			qtilde[r] = point.q[k][r] - 0.50 * phi[r] * temp[r];
		}
		qtilde_to_primitive(qtilde, prim);
		flux_Gzn(G_k, tan1, tan2, nor, prim);
		//
		for (int r = 0; r < 5; r++)
		{
			temp[r] = G_k[r] - G_i[r];
		}
		//
		for (int r = 0; r < 5; r++)
		{
			sum_delx_delf[r] = sum_delx_delf[r] + temp[r] * dels_weights;
			sum_dely_delf[r] = sum_dely_delf[r] + temp[r] * delt_weights;
			sum_delz_delf[r] = sum_delz_delf[r] + temp[r] * deln_weights;
		}
		//
	}
	//
	det = sum_delx_sqr * (sum_dely_sqr * sum_delz_sqr - sum_dely_delz * sum_dely_delz) - sum_delx_dely * (sum_delx_dely * sum_delz_sqr - sum_dely_delz * sum_delz_delx) + sum_delz_delx * (sum_delx_dely * sum_dely_delz - sum_dely_sqr * sum_delz_delx);
	//
	for (int r = 0; r < 5; r++)
	{
		temp[r] = sum_delx_sqr * (sum_dely_sqr * sum_delz_delf[r] - sum_dely_delf[r] * sum_dely_delz) - sum_delx_dely * (sum_delx_dely * sum_delz_delf[r] - sum_delx_delf[r] * sum_dely_delz) + sum_delz_delx * (sum_delx_dely * sum_dely_delf[r] - sum_delx_delf[r] * sum_dely_sqr);
	}
	//
	for (int r = 0; r < 5; r++)
	{
		G[r] = temp[r] / det;
	}
	//
}
//

__global__ void wall_dGz_neg_cuda(points &point, double power, double VL_CONST,double pi,int wall_points,int *wall_points_index)
//
{
	int ind = blockIdx.x * blockDim.x + threadIdx.x;

    if (ind < 0 || ind >= wall_points){
        return;
    }

    int i=wall_points_index[ind];
	int j,k;
	double prim[5];
	double x_i, y_i, z_i, x_k, y_k, z_k;
	double tan1[3], tan2[3], nor[3];
	double G_i[5], G_k[5];
	double delx, dely, delz, det;
	double dels, delt, deln;
	//
	double sum_delx_sqr, sum_dely_sqr, sum_delz_sqr;
	double sum_delx_dely, sum_dely_delz, sum_delz_delx;
	double sum_delx_delf[5] = {0}, sum_dely_delf[5] = {0}, sum_delz_delf[5] = {0};
	double dist, weights;
	double temp[5], qtilde[5], phi[5];
	double dels_weights, delt_weights, deln_weights;
	//
	//
	sum_delx_sqr = 0.00;
	sum_dely_sqr = 0.00;
	sum_delz_sqr = 0.00;
	//
	sum_delx_dely = 0.00;
	sum_dely_delz = 0.00;
	sum_delz_delx = 0.00;
	//
	//
	x_i = point.x[i];
	y_i = point.y[i];
	z_i = point.z[i];
	//
	for (int r = 0; r < 3; r++)
	{
		tan1[r] = point.tan1[i][r];
		tan2[r] = point.tan2[i][r];
		nor[r] = point.nor[i][r];
	}
	//
	for (j = 0; j < point.nbhs[i]; j++)
	//
	{
		k = point.conn[i][j];
		//
		x_k = point.x[k];
		y_k = point.y[k];
		z_k = point.z[k];
		//
		delx = x_k - x_i;
		dely = y_k - y_i;
		delz = z_k - z_i;
		//
		dels = delx * tan1[0] + dely * tan1[1] + delz * tan1[2];
		delt = delx * tan2[0] + dely * tan2[1] + delz * tan2[2];
		deln = delx * nor[0] + dely * nor[1] + delz * nor[2];
		//
		dist = sqrt(dels * dels + delt * delt + deln * deln);
		weights = 1.00 / (pow(dist, power));
		//
		dels_weights = dels * weights;
		delt_weights = delt * weights;
		deln_weights = deln * weights;
		//
		sum_delx_sqr = sum_delx_sqr + dels * dels_weights;
		sum_dely_sqr = sum_dely_sqr + delt * delt_weights;
		sum_delz_sqr = sum_delz_sqr + deln * deln_weights;
		//
		sum_delx_dely = sum_delx_dely + dels * delt_weights;
		sum_dely_delz = sum_dely_delz + delt * deln_weights;
		sum_delz_delx = sum_delz_delx + deln * dels_weights;
		//
		for (int r = 0; r < 5; r++)
		{
			temp[r] = delx * point.dq[i][0][r] + dely * point.dq[i][1][r] + delz * point.dq[i][2][r];
			qtilde[r] = point.q[i][r] - 0.50 * temp[r];
		}
		venkat_limiter_cuda(point, qtilde, phi, i, VL_CONST);
		for (int r = 0; r < 5; r++)
		{
			qtilde[r] = point.q[i][r] - 0.50 * phi[r] * temp[r];
		}
		qtilde_to_primitive_cuda(qtilde, prim);
		flux_Gzn_cuda(G_i, tan1, tan2, nor, prim,pi);
		//
		for (int r = 0; r < 5; r++)
		{
			temp[r] = delx * point.dq[k][0][r] + dely * point.dq[k][1][r] + delz * point.dq[k][2][r];
			qtilde[r] = point.q[k][r] - 0.50 * temp[r];
		}
		venkat_limiter_cuda(point, qtilde, phi, k, VL_CONST);
		for (int r = 0; r < 5; r++)
		{
			qtilde[r] = point.q[k][r] - 0.50 * phi[r] * temp[r];
		}
		qtilde_to_primitive_cuda(qtilde, prim);
		flux_Gzn_cuda(G_k, tan1, tan2, nor, prim,pi);
		//
		for (int r = 0; r < 5; r++)
		{
			temp[r] = G_k[r] - G_i[r];
		}
		//
		for (int r = 0; r < 5; r++)
		{
			sum_delx_delf[r] = sum_delx_delf[r] + temp[r] * dels_weights;
			sum_dely_delf[r] = sum_dely_delf[r] + temp[r] * delt_weights;
			sum_delz_delf[r] = sum_delz_delf[r] + temp[r] * deln_weights;
		}
		//
	}
	//
	det = sum_delx_sqr * (sum_dely_sqr * sum_delz_sqr - sum_dely_delz * sum_dely_delz) - sum_delx_dely * (sum_delx_dely * sum_delz_sqr - sum_dely_delz * sum_delz_delx) + sum_delz_delx * (sum_delx_dely * sum_dely_delz - sum_dely_sqr * sum_delz_delx);
	//
	for (int r = 0; r < 5; r++)
	{
		temp[r] = sum_delx_sqr * (sum_dely_sqr * sum_delz_delf[r] - sum_dely_delf[r] * sum_dely_delz) - sum_delx_dely * (sum_delx_dely * sum_delz_delf[r] - sum_delx_delf[r] * sum_dely_delz) + sum_delz_delx * (sum_delx_dely * sum_dely_delf[r] - sum_delx_delf[r] * sum_dely_sqr);
	}
	//
	for (int r = 0; r < 5; r++)
	{
		point.flux_res[i][r] += 2.00*point.delt[i]*temp[r] / det;
	}
	//
}
