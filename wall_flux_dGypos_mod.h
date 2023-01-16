

//	First written on 10.06.2021.
//
#pragma once
#include "data_structure_mod.h"
#include "octant_fluxes_mod.h"
#include "q_variables_mod.h"
#include "q_derivatives_mod.h"
#include "qvariables_to_primitive_variables_mod.h"
#include "limiters_mod.h"
#include <float.h>

//	This subroutine evaluates the wall flux derivative dGy_pos

void wall_dGy_pos(double *G, int i)
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
	for (j = 0; j < point.ypos_nbhs[i]; j++)
	//
	{
		k = point.ypos_conn[i][j];
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
		flux_Gwyp(G_i, tan1, tan2, nor, prim);
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
		flux_Gwyp(G_k, tan1, tan2, nor, prim);
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
		temp[r] = sum_delx_sqr*(sum_dely_delf[r]*sum_delz_sqr - sum_dely_delz*sum_delz_delf[r]) 
				- sum_delx_dely*(sum_delx_delf[r]*sum_delz_sqr - sum_delz_delx*sum_delz_delf[r]) 
				+ sum_delz_delx*(sum_delx_delf[r]*sum_dely_delz - sum_delz_delx*sum_dely_delf[r]);
	}
	//
	for (int r = 0; r < 5; r++)
	{
		G[r] = temp[r] / det;
	}
	//
}
//

__global__ void wall_dGy_pos_cuda(points &point,double power, double VL_CONST,double pi,int wall_points,int *wall_points_index)
//
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
	for (j = 0; j < point.ypos_nbhs[i]; j++)
	//
	{
		k = point.ypos_conn[i][j];
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
		// if(i==100001){
		// 	printf("x_i=%.15f y_i=%.15f z_i=%.15f\n",x_i,y_i,z_i);
		// 	printf("delx=%.15f dely=%.15f delz=%.15f\n",delx,dely,delz);
		// 	printf("dels=%.15f delt=%.15f deln=%.15f\n",dels,delt,deln);
		// }
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
		// if(i==100001){
		// 	printf("delt_weights=%.15f deln_weights=%.15f dels_weights=%.15f\n",delt_weights,deln_weights,dels_weights);
		// 	printf("sum_delx_sqr=%.15f sum_dely_sqr=%.15f sum_delz_sqr=%.15f\n",sum_delx_sqr,sum_dely_sqr,sum_delz_sqr);
		// 	printf("sum_delx_dely=%.15f sum_dely_delz=%.15f sum_delz_delx=%.15f\n",sum_delx_dely,sum_dely_delz,sum_delz_delx);
		// }
		//
		for (int r = 0; r < 5; r++)
		{
			temp[r] = delx * point.dq[i][0][r] + dely * point.dq[i][1][r] + delz * point.dq[i][2][r];
			qtilde[r] = point.q[i][r] - 0.50 * temp[r];
			// if(ind==0){
			// 	printf("i = %d\n",i);
			// 	printf("k = %d\n",k);
			// 	printf("delx=%.15f\n",delx);
			// 	printf("dely=%.15f\n",dely);
			// 	printf("delz=%.15f\n",delz);
			// 	printf("point.dq[%d][0][%d]=%.15f\n",i,r,point.dq[i][0][r]);
			// 	printf("point.dq[%d][1][%d]=%.15f\n",i,r,point.dq[i][1][r]);
			// 	printf("point.dq[%d][2][%d]=%.15f\n",i,r,point.dq[i][2][r]);
			// 	printf("qtilde[%d]=%.15f\n",r,qtilde[r]);
			// 	printf("temp[%d]=%.15f\n",r,temp[r]);
			// 	printf("point.q[%d][%d]=%.15f\n",i,r,point.q[i][r]);
			// }
		}
		venkat_limiter_cuda(point,qtilde, phi, i,VL_CONST);
		for (int r = 0; r < 5; r++)
		{
			qtilde[r] = point.q[i][r] - 0.50 * phi[r] * temp[r];
		}
		qtilde_to_primitive_cuda(qtilde, prim);
		flux_Gwyp_cuda(G_i, tan1, tan2, nor, prim,pi);
		//
		for (int r = 0; r < 5; r++)
		{
			temp[r] = delx * point.dq[k][0][r] + dely * point.dq[k][1][r] + delz * point.dq[k][2][r];
			qtilde[r] = point.q[k][r] - 0.50 * temp[r];
		}
		venkat_limiter_cuda(point,qtilde, phi, k,VL_CONST);
		for (int r = 0; r < 5; r++)
		{
			qtilde[r] = point.q[k][r] - 0.50 * phi[r] * temp[r];
		}
		qtilde_to_primitive_cuda(qtilde, prim);
		flux_Gwyp_cuda(G_k, tan1, tan2, nor, prim,pi);
		//
		for (int r = 0; r < 5; r++)
		{
			temp[r] = G_k[r] - G_i[r];
			// if(i==100001){
			// 	printf("temp[%d] = %.15f\n",r,temp[r]);
			// 	printf("G_k[%d] = %.15f\n",r,G_k[r]);
			// 	printf("G_i[%d] = %.15f\n",r,G_i[r]);
			// }
		}
		//
		for (int r = 0; r < 5; r++)
		{
			sum_delx_delf[r] = sum_delx_delf[r] + temp[r] * dels_weights;
			sum_dely_delf[r] = sum_dely_delf[r] + temp[r] * delt_weights;
			sum_delz_delf[r] = sum_delz_delf[r] + temp[r] * deln_weights;
			// if(i==100001){
			// 	printf("sum_delx_delf[%d] = %.15f\n",r,sum_delx_delf[r]);
			// 	printf("sum_dely_delf[%d] = %.15f\n",r,sum_dely_delf[r]);
			// 	printf("sum_delz_delf[%d] = %.15f\n",r,sum_delz_delf[r]);
			// 	printf("temp[%d] = %.15f\n",r,temp[r]);
			// 	printf("dels_weights = %.15f\n",dels_weights);
			// 	printf("delt_weights = %.15f\n",delt_weights);
			// 	printf("deln_weights = %.15f\n",deln_weights);
			// }
		}
		//
	}
	//
	det = sum_delx_sqr * (sum_dely_sqr * sum_delz_sqr - sum_dely_delz * sum_dely_delz) - sum_delx_dely * (sum_delx_dely * sum_delz_sqr - sum_dely_delz * sum_delz_delx) + sum_delz_delx * (sum_delx_dely * sum_dely_delz - sum_dely_sqr * sum_delz_delx);
	//
	for (int r = 0; r < 5; r++)
	{
		temp[r] = sum_delx_sqr*(sum_dely_delf[r]*sum_delz_sqr - sum_dely_delz*sum_delz_delf[r]) 
				- sum_delx_dely*(sum_delx_delf[r]*sum_delz_sqr - sum_delz_delx*sum_delz_delf[r]) 
				+ sum_delz_delx*(sum_delx_delf[r]*sum_dely_delz - sum_delz_delx*sum_dely_delf[r]);
		if(abs(temp[r]-0.0)<1e-15)
			temp[r]=0.0;
	}
	//
	// for (int r = 0; r < 5; r++)
	// {
	// 	temp[r] = temp[r]/det;
	// }
	for (int r = 0; r < 5; r++)
	{
		point.flux_res[i][r] += 2.00*point.delt[i]*temp[r]/det;
	}
	// if(i==100001){
    //     printf("i = %d\n",i);
    //     for (int r = 0; r < 5; r++)
    //     {
	// 		printf("New :\n");
	// 		printf("%.22f \n",temp[r]);
    //         printf("%.22f \n",det);
	// 		if(abs(temp[r] - det) < DBL_EPSILON)
	// 			printf("Equal\n");
	// 		else
	// 			printf("Not Equal\n");
	// 		printf("%.22f \n",temp[r]/det);
	// 		// printf("%.100f \n",point.delt[i]);
    //     }
    //     printf("\n");
    //     for(int r=0;r<5;r++)
    //         printf("%.14f ",point.flux_res[i][r]);
    //     printf("\n");
    // }
	//
}
//
//

