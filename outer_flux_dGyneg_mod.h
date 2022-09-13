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

//	First written on 10.06.2021.
//
#pragma once
#include "data_structure_mod.h"
#include "octant_fluxes_mod.h"
#include "q_variables_mod.h"
#include "q_derivatives_mod.h"
#include "qvariables_to_primitive_variables_mod.h"
#include "limiters_mod.h"

//	This subroutine evaluates the wall flux derivative dGy_neg

void outer_dGy_neg(double *G, int i)
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
		tan1[r] = point.tan1[r][i];
		tan2[r] = point.tan2[r][i];
		nor[r] = point.nor[r][i];
	}
	//
	for (j = 0; j < point.yneg_nbhs[i]; j++)
	//
	{
		k = point.yneg_conn[j][i];
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
			temp[r] = delx * point.dq[0][r][i] + dely * point.dq[1][r][i] + delz * point.dq[2][r][i];
			qtilde[r] = point.q[r][i] - 0.50 * temp[r];
		}
		venkat_limiter(qtilde, phi, i);
		for (int r = 0; r < 5; r++)
		{
			qtilde[r] = point.q[r][i] - 0.50 * phi[r] * temp[r];
		}
		qtilde_to_primitive(qtilde, prim);
		flux_Goyn(G_i, tan1, tan2, nor, prim);
		//
		for (int r = 0; r < 5; r++)
		{
			temp[r] = delx * point.dq[0][r][k] + dely * point.dq[1][r][k] + delz * point.dq[2][r][k];
			qtilde[r] = point.q[r][k] - 0.50 * temp[r];
		}
		venkat_limiter(qtilde, phi, k);
		for (int r = 0; r < 5; r++)
		{
			qtilde[r] = point.q[r][k] - 0.50 * phi[r] * temp[r];
		}
		qtilde_to_primitive(qtilde, prim);
		flux_Goyn(G_k, tan1, tan2, nor, prim);
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



__global__ void outer_dGy_neg_cuda(points &point,double power, double VL_CONST,double pi,int outer_points,int *outer_points_index)
//
//
{
	int ind = blockIdx.x * blockDim.x + threadIdx.x;
	if (ind < 0 || ind >= outer_points){
        return;
    }

    int i=outer_points_index[ind];
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
		tan1[r] = point.tan1[r][i];
		tan2[r] = point.tan2[r][i];
		nor[r] = point.nor[r][i];
	}
	//
	for (j = 0; j < point.yneg_nbhs[i]; j++)
	//
	{
		k = point.yneg_conn[j][i];
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
			temp[r] = delx * point.dq[0][r][i] + dely * point.dq[1][r][i] + delz * point.dq[2][r][i];
			qtilde[r] = point.q[r][i] - 0.50 * temp[r];
		}
		venkat_limiter_cuda(point,qtilde, phi, i,VL_CONST);
		for (int r = 0; r < 5; r++)
		{
			qtilde[r] = point.q[r][i] - 0.50 * phi[r] * temp[r];
		}
		qtilde_to_primitive_cuda(qtilde, prim);
		flux_Goyn_cuda(G_i, tan1, tan2, nor, prim,pi);
		//
		for (int r = 0; r < 5; r++)
		{
			temp[r] = delx * point.dq[0][r][k] + dely * point.dq[1][r][k] + delz * point.dq[2][r][k];
			qtilde[r] = point.q[r][k] - 0.50 * temp[r];
		}
		venkat_limiter_cuda(point,qtilde, phi, k,VL_CONST);
		for (int r = 0; r < 5; r++)
		{
			qtilde[r] = point.q[r][k] - 0.50 * phi[r] * temp[r];
		}
		qtilde_to_primitive_cuda(qtilde, prim);
		flux_Goyn_cuda(G_k, tan1, tan2, nor, prim,pi);
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
		point.flux_res[r][i] += temp[r]*point.delt[i] / det;
	}
	//
}
//

__global__ void outer_dGy_neg_multi_nccl(int myRank, splitPoints *splitPoint, double power, double VL_CONST, double pi, int outerPointsLocal, int *outerPointsLocalIndex,int* globalToLocalIndex,int **globalToGhostIndex,transferPoints **receiveBuffer,int *partVector)
//
//
{
	int ind = blockIdx.x * blockDim.x + threadIdx.x;
	if (ind < 0 || ind >= outerPointsLocal){
        return;
    }

    int i=outerPointsLocalIndex[ind];
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
	for(int r=0;r<5;r++){
		sum_delx_delf[r]=0;
		sum_dely_delf[r]=0;
		sum_delz_delf[r]=0;
	}
    //
    x_i = splitPoint[i].x;
    y_i = splitPoint[i].y;
    z_i = splitPoint[i].z;
    //
    for (int r = 0; r < 3; r++)
    {
        tan1[r] = splitPoint[i].tan1[r];
        tan2[r] = splitPoint[i].tan2[r];
        nor[r] = splitPoint[i].nor[r];
    }
    //
	// if(ind==0 && myRank==0){
	// 	printf("GLOBAL INDEX %d\n",splitPoint[i].globalIndex);
	// }
    for (j = 0; j < splitPoint[i].numberOfLocalynegNbhs; j++)
    //
    {
        k = splitPoint[i].localyneg_conn[j];
        k=globalToLocalIndex[k];
        //
        x_k = splitPoint[k].x;
        y_k = splitPoint[k].y;
        z_k = splitPoint[k].z;
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
		// 		if(ind==0 && myRank==0){
        //         printf("INDEX %d\n",splitPoint[k].globalIndex);
        //         printf("sum_delx_sqr %.15f\n",sum_delx_sqr);
        //         printf("sum_dely_sqr %.15f\n",sum_dely_sqr);
        //         printf("sum_delz_sqr %.15f\n",sum_delz_sqr);
		// 		printf("sum_delx_dely %.15f\n",sum_delx_dely);
        //         printf("sum_dely_delz %.15f\n",sum_dely_delz);
        //         printf("sum_delz_delx %.15f\n",sum_delz_delx);
        //         printf("dels %.15f\n",dels);
        //         printf("delt %.15f\n",delt);
        //         printf("deln %.15f\n",deln);
        //         printf("dels_weights %.15f\n",dels_weights);
        //         printf("delt_weights %.15f\n",delt_weights);
        //         printf("deln_weights %.15f\n",deln_weights);
        //         printf("delx %.15f\n",delx);
        //         printf("dely %.15f\n",dely);
        //         printf("delz %.15f\n",delz);
		// 		printf("tan1 %.15f\n",tan2[0]);
        //         printf("tan1 %.15f\n",tan2[1]);
        //         printf("tan1 %.15f\n",tan2[2]);
        // }
		//
        for (int r = 0; r < 5; r++)
        {
            temp[r] = delx * splitPoint[i].dq[0][r] + dely * splitPoint[i].dq[1][r] + delz * splitPoint[i].dq[2][r];
            qtilde[r] = splitPoint[i].q[r] - 0.50 * temp[r];
        }
		venkat_limiter_multi_nccl(splitPoint,qtilde, phi, i,VL_CONST);
        for (int r = 0; r < 5; r++)
        {
            qtilde[r] = splitPoint[i].q[r] - 0.50 * phi[r] * temp[r];
        }
        qtilde_to_primitive_cuda(qtilde, prim);
        flux_Goyn_cuda(G_i, tan1, tan2, nor, prim,pi);
        //
        for (int r = 0; r < 5; r++)
        {
            temp[r] = delx * splitPoint[k].dq[0][r] + dely * splitPoint[k].dq[1][r] + delz * splitPoint[k].dq[2][r];
            qtilde[r] = splitPoint[k].q[r] - 0.50 * temp[r];
			// if(splitPoint[i].globalIndex== 699502){
			// 	printf("1 Local CHECK %d\n",splitPoint[k].globalIndex);
            //     printf("1 Local splitPoint[k].dq[0][r] %.25f\n",splitPoint[k].dq[0][r]);
            //     printf("1 splitPoint[k].dq[1][r] %.15f\n",splitPoint[k].dq[1][r]);
            //     printf("1 splitPoint[k].dq[2][r] %.15f\n",splitPoint[k].dq[2][r]);
            // }
        }
		venkat_limiter_multi_nccl(splitPoint,qtilde, phi, k,VL_CONST);
        for (int r = 0; r < 5; r++)
        {
            qtilde[r] = splitPoint[k].q[r] - 0.50 * phi[r] * temp[r];
			// if(splitPoint[i].globalIndex==699502){
			// 	printf("splitPoint[k].q[r] %.15f\n",splitPoint[k].q[r]);
        	// 	printf("phi[r] %.15f\n",phi[r]);
			// 	printf("temp[r] %.15f\n",temp[r]);
       		// }
        }
        qtilde_to_primitive_cuda(qtilde, prim);
        flux_Goyn_cuda(G_k, tan1, tan2, nor, prim,pi);
        //
        for (int r = 0; r < 5; r++)
        {
            temp[r] = G_k[r] - G_i[r];
			// if(splitPoint[i].globalIndex==699502){
			// 	printf("G_k[r] %.15f\n",G_k[r]);
        	// 	printf("G_i[r] %.15f\n",G_i[r]);
       		// }
        }
        //
        for (int r = 0; r < 5; r++)
        {
            sum_delx_delf[r] = sum_delx_delf[r] + temp[r] * dels_weights;
            sum_dely_delf[r] = sum_dely_delf[r] + temp[r] * delt_weights;
            sum_delz_delf[r] = sum_delz_delf[r] + temp[r] * deln_weights;
			// if(splitPoint[i].globalIndex==699502){
			// 	printf("Local sum_dely_delf %.15f\n",sum_dely_delf[r]);
        	// 	printf("Local sum_delx_delf %.15f\n",sum_delx_delf[r]);
        	// 	printf("Local sum_delz_delf %.15f\n",sum_delz_delf[r]);
			// 	printf("Local dels_weights %.15f\n",dels_weights);
			// 	printf("Local delt_weights %.15f\n",delt_weights);
			// 	printf("Local deln_weights %.15f\n",deln_weights);
			// 	printf("Local temp[r] %.15f\n",temp[r]);
			// 	// printf("Local sum_delx_delf %.15f\n",sum_delx_delf[r]);
       		// }
        }
        //
    }
    //
    for (j = 0; j < splitPoint[i].numberOfGhostynegNbhs; j++)
    //
    {
        k = splitPoint[i].ghostyneg_conn[j];
        int device=partVector[k];
        int ghostIndex=globalToGhostIndex[device][k];
        //
        x_k = receiveBuffer[device][ghostIndex].x;
        y_k = receiveBuffer[device][ghostIndex].y;
        z_k = receiveBuffer[device][ghostIndex].z;
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
		// if(ind==0 && myRank==0){
        //         printf("INDEX %d\n",receiveBuffer[device][ghostIndex].globalIndex);
        //         printf("sum_delx_sqr %.15f\n",sum_delx_sqr);
        //         printf("sum_dely_sqr %.15f\n",sum_dely_sqr);
        //         printf("sum_delz_sqr %.15f\n",sum_delz_sqr);
        //         printf("dels %.15f\n",dels);
        //         printf("delt %.15f\n",delt);
        //         printf("deln %.15f\n",deln);
        //         printf("dels_weights %.15f\n",dels_weights);
        //         printf("delt_weights %.15f\n",delt_weights);
        //         printf("deln_weights %.15f\n",deln_weights);
        //         printf("delx %.15f\n",delx);
        //         printf("dely %.15f\n",dely);
        //         printf("delz %.15f\n",delz);
		// 		printf("tan1 %.15f\n",tan2[0]);
        //         printf("tan1 %.15f\n",tan2[1]);
        //         printf("tan1 %.15f\n",tan2[2]);
        // }
		//
        for (int r = 0; r < 5; r++)
        {
            temp[r] = delx * splitPoint[i].dq[0][r] + dely * splitPoint[i].dq[1][r] + delz * splitPoint[i].dq[2][r];
            qtilde[r] = splitPoint[i].q[r] - 0.50 * temp[r];
        }
		venkat_limiter_multi_nccl(splitPoint,qtilde, phi, i,VL_CONST);
        for (int r = 0; r < 5; r++)
        {
            qtilde[r] = splitPoint[i].q[r] - 0.50 * phi[r] * temp[r];
        }
        qtilde_to_primitive_cuda(qtilde, prim);
        flux_Goyn_cuda(G_i, tan1, tan2, nor, prim,pi);
        //
        for (int r = 0; r < 5; r++)
        {
            temp[r] = delx * receiveBuffer[device][ghostIndex].dq[0][r] + dely * receiveBuffer[device][ghostIndex].dq[1][r] + delz * receiveBuffer[device][ghostIndex].dq[2][r];
            qtilde[r] = receiveBuffer[device][ghostIndex].q[r] - 0.50 * temp[r];
        }
        venkat_limiter_multi_nccl_ghost(receiveBuffer, qtilde, phi, device,ghostIndex, VL_CONST);
        for (int r = 0; r < 5; r++)
        {
            qtilde[r] = receiveBuffer[device][ghostIndex].q[r] - 0.50 * phi[r] * temp[r];
        }
        qtilde_to_primitive_cuda(qtilde, prim);
        flux_Goyn_cuda(G_k, tan1, tan2, nor, prim,pi);
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
			// if(splitPoint[i].globalIndex==699502){
			// 	printf("Ghost sum_dely_delf %.15f\n",sum_dely_delf[r]);
        	// 	printf("Ghost sum_delx_delf %.15f\n",sum_delx_delf[r]);
        	// 	printf("Ghost sum_delx_delf %.15f\n",sum_delx_delf[r]);
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
	// 	if(i==249755 && myRank==1){
	// 		printf("sum_dely_delf %.15f\n",sum_dely_delf[r]);
    //     	printf(" sum_delx_delf %.15f\n",sum_delx_delf[r]);
    //     	printf(" sum_delx_delf %.15f\n",sum_delx_delf[r]);
    //    }
	}
    //
    for (int r = 0; r < 5; r++)
    {
       splitPoint[i].flux_res[r] += temp[r]*splitPoint[i].delt / det;
	//     if(splitPoint[i].globalIndex==699502){
	// 		printf("Global Index %d %d\n",splitPoint[i].globalIndex,i);
    //     	printf(" splitPoint[i].flux_res[r] %.25f\n",splitPoint[i].flux_res[r]);
    //     	printf(" temp %.25f\n",temp[r]);
    //     	printf(" splitPoint[i].delt %.15f\n",splitPoint[i].delt);
    //     	printf(" det %.20f\n",det);
    //    }
    }
    //
}