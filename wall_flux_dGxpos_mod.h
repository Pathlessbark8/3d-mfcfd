

//	First written on 10.06.2021.
//
#pragma once
#include "data_structure_mod.h"
#include "octant_fluxes_mod.h"
#include "q_variables_mod.h"
#include "q_derivatives_mod.h"
#include "qvariables_to_primitive_variables_mod.h"
#include "limiters_mod.h"

//	This subroutine evaluates the wall flux derivative dGx_pos

void wall_dGx_pos(double *G, int i)
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

    for (j = 0; j < point.xpos_nbhs[i]; j++)
    //
    {
        k = point.xpos_conn[j][i];
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
        flux_Gwxp(G_i, tan1, tan2, nor, prim);
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
        flux_Gwxp(G_k, tan1, tan2, nor, prim);
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
        temp[r] = sum_delx_delf[r] * (sum_dely_sqr * sum_delz_sqr - sum_dely_delz * sum_dely_delz) - sum_dely_delf[r] * (sum_delx_dely * sum_delz_sqr - sum_delz_delx * sum_dely_delz) + sum_delz_delf[r] * (sum_delx_dely * sum_dely_delz - sum_delz_delx * sum_dely_sqr);
    }
    //
    for (int r = 0; r < 5; r++)
    {
        G[r] = temp[r] / det;
    }
    //
}
//

__global__ void wall_dGx_pos_cuda(points &point, double power, double VL_CONST, double pi, int wall_points, int *wall_points_index)
//
//
{
    int ind = blockIdx.x * blockDim.x + threadIdx.x;

    if (ind < 0 || ind >= wall_points)
    {
        return;
    }

    int i = wall_points_index[ind];
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
    for (int r = 0; r < 5; r++)
    {
        sum_delx_delf[r] = 0;
        sum_dely_delf[r] = 0;
        sum_delz_delf[r] = 0;
    }
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

    for (j = 0; j < point.xpos_nbhs[i]; j++)
    //
    {
        k = point.xpos_conn[j][i];
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

        venkat_limiter_cuda(point, qtilde, phi, i, VL_CONST);
        for (int r = 0; r < 5; r++)
        {
            qtilde[r] = point.q[r][i] - 0.50 * phi[r] * temp[r];
        }
        qtilde_to_primitive_cuda(qtilde, prim);
        flux_Gwxp_cuda(G_i, tan1, tan2, nor, prim, pi);
        //
        for (int r = 0; r < 5; r++)
        {
            temp[r] = delx * point.dq[0][r][k] + dely * point.dq[1][r][k] + delz * point.dq[2][r][k];
            qtilde[r] = point.q[r][k] - 0.50 * temp[r];
        }
        venkat_limiter_cuda(point, qtilde, phi, k, VL_CONST);
        for (int r = 0; r < 5; r++)
        {
            qtilde[r] = point.q[r][k] - 0.50 * phi[r] * temp[r];
        }
        qtilde_to_primitive_cuda(qtilde, prim);
        flux_Gwxp_cuda(G_k, tan1, tan2, nor, prim, pi);
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
        temp[r] = sum_delx_delf[r] * (sum_dely_sqr * sum_delz_sqr - sum_dely_delz * sum_dely_delz) - sum_dely_delf[r] * (sum_delx_dely * sum_delz_sqr - sum_delz_delx * sum_dely_delz) + sum_delz_delf[r] * (sum_delx_dely * sum_dely_delz - sum_delz_delx * sum_dely_sqr);
    }
    //
    for (int r = 0; r < 5; r++)
    {
        point.flux_res[r][i] = 2.00 *point.delt[i]* temp[r] / det;
    }
    //
}

__global__ void wall_dGx_pos_multi_nccl(int myRank,splitPoints *splitPoint, double power, double VL_CONST, double pi, int wallPointsLocal, int *wallPointsLocalIndex,int* globalToLocalIndex,int **globalToGhostIndex,transferPoints **receiveBuffer,int *partVector)
//
//
{
    int ind = blockIdx.x * blockDim.x + threadIdx.x;

    if (ind < 0 || ind >= wallPointsLocal)
    {
        return;
    }
    
    int i = wallPointsLocalIndex[ind];
    // if(myRank==1 && i==5150){
	// 	printf("Global 2 new Rank %d\n",splitPoint[i].globalIndex);
    //     printf("delt=%.15f\n",splitPoint[i].delt);
    //     printf("Min dist=%.15f\n",splitPoint[i].min_dist);
    //     printf("Status =%d\n",splitPoint[i].status);
	// }
    // if(ind==0){
        // printf("%d\n",i);
        // printf("Global Index=%d\n",splitPoint[i].globalIndex);
    // }

    // if(myRank==1 && i==6081){
    //     printf("%d\n",i);
    //     printf("Global Index=%d\n",splitPoint[i].globalIndex);
    // }
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

    for (j = 0; j < splitPoint[i].numberOfLocalxposNbhs; j++)
    //
    {
        int nbh = splitPoint[i].localxpos_conn[j];
        k=globalToLocalIndex[nbh];
        // printf("k = %d\n",k);

        //
        x_k = splitPoint[k].x;
        y_k = splitPoint[k].y;
        z_k = splitPoint[k].z;
        //
        delx = x_k - x_i;
        dely = y_k - y_i;
        delz = z_k - z_i;
        //
        // if(myRank==1 && i==5150){
        //     printf("Global Index %d\n",splitPoint[i].globalIndex);
        //     printf("count =%d %d\n", splitPoint[i].numberOfLocalxposNbhs,splitPoint[i].numberOfGhostxposNbhs);
        //     printf("Local nbh = %d\n",nbh);
        //     printf("delx = %.15f\n",delx);
        //     printf("dely = %.15f\n",dely);
        //     printf("delz = %.15f\n",delz);
        // }
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
        // if(myRank==1 && i==6081){
        //     printf("sum_delx_sqr = %.15f\n",sum_delx_sqr);
        //     printf("sum_dely_sqr = %.15f\n",sum_dely_sqr);
        //     printf("sum_delz_sqr = %.15f\n",sum_delz_sqr);
        // }
        sum_delx_dely = sum_delx_dely + dels * delt_weights;
        sum_dely_delz = sum_dely_delz + delt * deln_weights;
        sum_delz_delx = sum_delz_delx + deln * dels_weights;
        //
        // if(myRank==1 && i==6081){
        //     printf("sum_delx_dely = %.15f\n",sum_delx_dely);
        //     printf("sum_dely_delz = %.15f\n",sum_dely_delz);
        //     printf("sum_delz_delx = %.15f\n",sum_delz_delx);
        // }
        for (int r = 0; r < 5; r++)
        {
            temp[r] = delx * splitPoint[i].dq[0][r] + dely * splitPoint[i].dq[1][r] + delz * splitPoint[i].dq[2][r];
            qtilde[r] = splitPoint[i].q[r] - 0.50 * temp[r];
            // if(myRank==1 && i==6081)
            // printf("temp[%d]=%f %f %f %f\n",r,temp[r],splitPoint[i].dq[0][r],splitPoint[i].dq[1][r],splitPoint[i].dq[2][r]);
        }

        venkat_limiter_multi_nccl(splitPoint, qtilde, phi, i, VL_CONST);
        for (int r = 0; r < 5; r++)
        {
            qtilde[r] = splitPoint[i].q[r] - 0.50 * phi[r] * temp[r];
        }
        qtilde_to_primitive_cuda(qtilde, prim);
        flux_Gwxp_cuda(G_i, tan1, tan2, nor, prim, pi);
        //
        for (int r = 0; r < 5; r++)
        {
            temp[r] = delx * splitPoint[k].dq[0][r] + dely * splitPoint[k].dq[1][r] + delz * splitPoint[k].dq[2][r];
            qtilde[r] = splitPoint[k].q[r] - 0.50 * temp[r];
            // if(myRank==1 && i==6081 && nbh==421000){
            //     printf("k = %d\n",nbh);
            //     printf("splitPoint.q[%d]=%.15f\n",r,splitPoint[k].q[r]);
            //     printf("splitPoint.dq[0][%d]=%.15f\n",r,splitPoint[k].dq[0][r]);
            //     printf("splitPoint.dq[1][%d]=%.15f\n",r,splitPoint[k].dq[1][r]);
            //     printf("splitPoint.dq[2][%d]=%.15f\n",r,splitPoint[k].dq[2][r]);
            // }
        }
        venkat_limiter_multi_nccl(splitPoint, qtilde, phi, k, VL_CONST);
        for (int r = 0; r < 5; r++)
        {
            qtilde[r] = splitPoint[k].q[r] - 0.50 * phi[r] * temp[r];
            // if(myRank==1 && i==6081 && nbh==421000){
			// 	printf("qtilde[%d] %.15f point[%d] %.15f phi[%d] %.15f temp[%d] %.15f \n",r,qtilde[r],k,splitPoint[k].q[r],r,phi[r],r,temp[r]);			
            // }
        }
        qtilde_to_primitive_cuda(qtilde, prim);
        flux_Gwxp_cuda(G_k, tan1, tan2, nor, prim, pi);
        //
        for (int r = 0; r < 5; r++)
        {
            
            temp[r] = G_k[r] - G_i[r];
            // if(myRank==1 && i==6081 && nbh==421000)
            //    printf("Local G_k[%d] %.15f G_i[%d] %.15f temp[%d] %.15f\n",r,G_k[r],r,G_i[r],r,temp[r]);
        }
        //
        for (int r = 0; r < 5; r++)
        {
            // if(myRank==1 && i==6081){
            //     printf("k = %d\n",nbh);
            //     printf("sum_delx_delf[%d] = %.15f\n",r,sum_delx_delf[r]);
            //     printf("sum_dely_delf[%d] = %.15f\n",r,sum_dely_delf[r]);
            //     printf("sum_delz_delf[%d] = %.15f\n",r,sum_delz_delf[r]);
            // }
            sum_delx_delf[r] = sum_delx_delf[r] + temp[r] * dels_weights;
            sum_dely_delf[r] = sum_dely_delf[r] + temp[r] * delt_weights;
            sum_delz_delf[r] = sum_delz_delf[r] + temp[r] * deln_weights;

            // if(myRank==1 && i==6081){
            //     printf("Local k = %d\n",nbh);
            //     printf("sum_delx_delf[%d] = %.15f\n",r,sum_delx_delf[r]);
            //     printf("sum_dely_delf[%d] = %.15f\n",r,sum_dely_delf[r]);
            //     printf("sum_delz_delf[%d] = %.15f\n",r,sum_delz_delf[r]);
            // }
            // if(myRank==0 && i==53486)
            //     printf("sum_delx_delf[%d]=%f\n",r,sum_delx_delf[r]);
            // if(myRank==1 && i==6081){
            //     printf("k = %d\n",nbh);
            //     printf("sum_delx_delf[%d]=%.15f\n",r,sum_delx_delf[r]);
            //     printf("temp[%d]=%.15f\n",r,temp[r]);
            //     printf("dels_weights=%.15f\n",dels_weights);
            //     printf("delt_weights=%.15f\n",delt_weights);
            //     printf("deln_weights=%.15f\n",deln_weights);
            // }


        }
        //
    }
    
    for (j = 0; j < splitPoint[i].numberOfGhostxposNbhs; j++)
    //
    {
        k = splitPoint[i].ghostxpos_conn[j];
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
        // if(myRank==1 && i==5150){
        //     printf("Global Index %d\n",splitPoint[i].globalIndex);
        //     printf("count =%d %d\n", splitPoint[i].numberOfLocalxposNbhs,splitPoint[i].numberOfGhostxposNbhs);
        //     printf("Ghost nbh = %d\n",receiveBuffer[device][ghostIndex].globalIndex);
        //     printf("delx = %.15f\n",delx);
        //     printf("dely = %.15f\n",dely);
        //     printf("delz = %.15f\n",delz);
        //     printf("x_k = %.15f\n",x_k);
        //     printf("y_k = %.15f\n",y_k);
        //     printf("z_k = %.15f\n",z_k);
        // }
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
            temp[r] = delx * splitPoint[i].dq[0][r] + dely * splitPoint[i].dq[1][r] + delz * splitPoint[i].dq[2][r];
            qtilde[r] = splitPoint[i].q[r] - 0.50 * temp[r];
            //  if(myRank==1 && i==6081)
            // printf("temp[%d]=%f %f %f %f\n",r,temp[r],splitPoint[i].dq[0][r],splitPoint[i].dq[1][r],splitPoint[i].dq[2][r]);
        }

        venkat_limiter_multi_nccl(splitPoint, qtilde, phi, i, VL_CONST);
        for (int r = 0; r < 5; r++)
        {
            qtilde[r] = splitPoint[i].q[r] - 0.50 * phi[r] * temp[r];
        }
        qtilde_to_primitive_cuda(qtilde, prim);
        flux_Gwxp_cuda(G_i, tan1, tan2, nor, prim, pi);
        //
        for (int r = 0; r < 5; r++)
        {
            temp[r] = delx * receiveBuffer[device][ghostIndex].dq[0][r] + dely * receiveBuffer[device][ghostIndex].dq[1][r] + delz * receiveBuffer[device][ghostIndex].dq[2][r];
            qtilde[r] = receiveBuffer[device][ghostIndex].q[r] - 0.50 * temp[r];
            // if(myRank==1 && i==6081){
            //     printf("k = %d\n",k);
            //     printf("receiveBuffer[%d][%d].q[%d]=%.15f\n",device,ghostIndex,r,receiveBuffer[device][ghostIndex].q[r]);
            //     printf("receiveBuffer[%d][%d].dq[0][%d]=%.15f\n",device,ghostIndex,r,receiveBuffer[device][ghostIndex].dq[0][r]);
            //     printf("receiveBuffer[%d][%d].dq[1][%d]=%.15f\n",device,ghostIndex,r,receiveBuffer[device][ghostIndex].dq[1][r]);
            //     printf("receiveBuffer[%d][%d].dq[2][%d]=%.15f\n",device,ghostIndex,r,receiveBuffer[device][ghostIndex].dq[2][r]);
            // }
        }
        venkat_limiter_multi_nccl_ghost(receiveBuffer, qtilde, phi, device,ghostIndex, VL_CONST);
        for (int r = 0; r < 5; r++)
        {
            qtilde[r] = receiveBuffer[device][ghostIndex].q[r] - 0.50 * phi[r] * temp[r];
        }
        qtilde_to_primitive_cuda(qtilde, prim);
        flux_Gwxp_cuda(G_k, tan1, tan2, nor, prim, pi);
        //
        for (int r = 0; r < 5; r++)
        {
            temp[r] = G_k[r] - G_i[r];
            // if(myRank==1 && i==6081)
            //    printf("Ghost G_k[%d] %.15f G_i[%d] %.15f temp[%d] %.15f\n",r,G_k[r],r,G_i[r],r,temp[r]);
        }
        //
        for (int r = 0; r < 5; r++)
        {
            // if(myRank==1 && i==6081){
            //     printf("k = %d\n",k);
            //     printf("sum_delx_delf[%d] = %.15f\n",r,sum_delx_delf[r]);
            //     printf("sum_dely_delf[%d] = %.15f\n",r,sum_dely_delf[r]);
            //     printf("sum_delz_delf[%d] = %.15f\n",r,sum_delz_delf[r]);
            // }
            sum_delx_delf[r] = sum_delx_delf[r] + temp[r] * dels_weights;
            sum_dely_delf[r] = sum_dely_delf[r] + temp[r] * delt_weights;
            sum_delz_delf[r] = sum_delz_delf[r] + temp[r] * deln_weights;

            // if(myRank==1 && i==6081){
            //     printf("Ghost k = %d\n",k);
            //     printf("sum_delx_delf[%d] = %.15f\n",r,sum_delx_delf[r]);
            //     printf("sum_dely_delf[%d] = %.15f\n",r,sum_dely_delf[r]);
            //     printf("sum_delz_delf[%d] = %.15f\n",r,sum_delz_delf[r]);
            // }
            // if(myRank==1 && i==6081){
            //     printf("Ghost sum_delx_delf[%d]=%.15f\n",r,sum_delx_delf[r]);
            // }
            // if(myRank==1 && i==6081){
            //     printf("k = %d\n",k);
            //     printf("sum_delx_delf[%d]=%.15f\n",r,sum_delx_delf[r]);
            //     printf("temp[%d]=%.15f\n",r,temp[r]);
            //     printf("dels_weights=%.15f\n",dels_weights);
            //     printf("delt_weights=%.15f\n",delt_weights);
            //     printf("deln_weights=%.15f\n",deln_weights);
            // }
        }
        //
    }
    
    det = sum_delx_sqr * (sum_dely_sqr * sum_delz_sqr - sum_dely_delz * sum_dely_delz) - sum_delx_dely * (sum_delx_dely * sum_delz_sqr - sum_dely_delz * sum_delz_delx) + sum_delz_delx * (sum_delx_dely * sum_dely_delz - sum_dely_sqr * sum_delz_delx);
    //
    for (int r = 0; r < 5; r++)
    {
        temp[r] = sum_delx_delf[r] * (sum_dely_sqr * sum_delz_sqr - sum_dely_delz * sum_dely_delz) - sum_dely_delf[r] * (sum_delx_dely * sum_delz_sqr - sum_delz_delx * sum_dely_delz) + sum_delz_delf[r] * (sum_delx_dely * sum_dely_delz - sum_delz_delx * sum_dely_sqr);
    }
    //
    for (int r = 0; r < 5; r++)
    {
        splitPoint[i].flux_res[r] = 2.00 *splitPoint[i].delt* temp[r] / det;
        // if(myRank==0 && i==53486){
        // // printf("delt=%f\n",splitPoint[i].delt);
        // // printf("det=%f\n",det);
        // printf("splitPoint[%d].flux_res[%d]=%f\n",i,r,splitPoint[i].flux_res[r]);
        //         // printf("delt[%d]=%f\n",r,sum_delx_delf[r]);
        // }

        // if(myRank==1 && i==5150){
        //     printf("splitPoint[%d].delt=%.15f\n",i,splitPoint[i].delt);
        //     printf("det=%.15f\n",det);
        //     printf("splitPoint[%d].flux_res[%d]=%.15f\n",i,r,splitPoint[i].flux_res[r]);
        // }

    }
    
}