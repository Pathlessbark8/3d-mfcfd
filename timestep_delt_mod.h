#pragma once
//
//
#include "data_structure_mod.h"
//

void timestep_delt()
//
//
{
    //
    int i, k, r;
    double rho, u1, u2, u3, pr;
    double delx, dely, delz, dist;
    double min_delt;
    double mod_u, delt;
    //
    //
    for (i = 0; i < max_points; i++)
    {
        //
        min_delt = 1.00;
        //
        for (r = 0; r < point.nbhs[i]; r++)
        //
        {
            k = point.conn[r][i];
            //
            rho = point.prim[0][k];
            u1 = point.prim[1][k];
            u2 = point.prim[2][k];
            u3 = point.prim[3][k];
            pr = point.prim[4][k];
            //
            delx = point.x[k] - point.x[i];
            dely = point.y[k] - point.y[i];
            delz = point.z[k] - point.z[i];
            //
            dist = sqrt(delx * delx + dely * dely + delz * delz);
            //
            mod_u = sqrt(u1 * u1 + u2 * u2 + u3 * u3);

            delt = CFL * dist / (mod_u + 3.00 * sqrt(pr / rho));
            //
            if (delt < min_delt)
            {
                min_delt = delt;
            }
            //
        }
        //
        point.delt[i] = min_delt;
        //
    }
    //
}
//
__global__ void timestep_delt_cuda(points &point,double CFL)
//
//
{
    //
    int k, r;
    double rho, u1, u2, u3, pr;
    double delx, dely, delz, dist;
    double min_delt;
    double mod_u, delt;
    //
    int bx = blockIdx.x;
    int tx = threadIdx.x;
    int i = bx * blockDim.x + tx;
    if (i < 0 || i >= max_points)
    {
        return;
    }
    //
    //
    min_delt = 1.00;
    //
    for (r = 0; r < point.nbhs[i]; r++)
    //
    {
        k = point.conn[r][i];
        //
        rho = point.prim[0][k];
        u1 = point.prim[1][k];
        u2 = point.prim[2][k];
        u3 = point.prim[3][k];
        pr = point.prim[4][k];
        //
        delx = point.x[k] - point.x[i];
        dely = point.y[k] - point.y[i];
        delz = point.z[k] - point.z[i];
        //
        dist = sqrt(delx * delx + dely * dely + delz * delz);
        //
        mod_u = sqrt(u1 * u1 + u2 * u2 + u3 * u3);

        delt = CFL * dist / (mod_u + 3.00 * sqrt(pr / rho));
        //
        if (delt < min_delt)
        {
            min_delt = delt;
        }
        //
    }
    //
    point.delt[i] = min_delt;
    //
    //
}
//

__global__ void timestep_delt_multi_nccl(int myRank,splitPoints *splitPoint,double CFL,int max_points_on_device,int *globalToLocalIndex,int **globalToGhostIndex,transferPoints** receiveBuffer,int *partVector,transferPoints** sendBuffer)
//
//
{
    //
    int k, r;
    double rho, u1, u2, u3, pr;
    double delx, dely, delz, dist;
    double min_delt;
    double mod_u, delt;
    //
    int bx = blockIdx.x;
    int tx = threadIdx.x;
    int i = bx * blockDim.x + tx;
    if (i < 0 || i >= max_points_on_device)
    {
        return;
    }
    //
    //
    min_delt = 1.00;
    //
    for (r = 0; r < splitPoint[i].numberOfLocalNbhs; r++)
    //
    {
        int nbh = splitPoint[i].localNbhs[r];
        k=globalToLocalIndex[nbh];
        //
        // if(myRank==1 && i==5150){
        //     printf("Local nbh %d \n",nbh);
        // }
        rho = splitPoint[k].prim[0];
        u1 = splitPoint[k].prim[1];
        u2 = splitPoint[k].prim[2];
        u3 = splitPoint[k].prim[3];
        pr = splitPoint[k].prim[4];
        //
        // if(myRank==1 && i==5150){
        //     printf("Local nbh %d \n",nbh);
        //     printf("rho %f \n",rho);
        //     printf("u1 %f \n",u1);
        //     printf("u2 %f \n",u2);
        //     printf("u3 %f \n",u3);
        //     printf("pr %f \n",pr);
        // }
        //
        delx = splitPoint[k].x - splitPoint[i].x;
        dely = splitPoint[k].y - splitPoint[i].y;
        delz = splitPoint[k].z - splitPoint[i].z;
        //
        dist = sqrt(delx * delx + dely * dely + delz * delz);
        //
        mod_u = sqrt(u1 * u1 + u2 * u2 + u3 * u3);

        delt = CFL * dist / (mod_u + 3.00 * sqrt(pr / rho));
        //
        if (delt < min_delt)
        {
            min_delt = delt;
        }
        //
    }
    //
    for (r = 0; r < splitPoint[i].numberOfGhostNbhs; r++)
    //
    {
        k = splitPoint[i].ghostNbhs[r];
        int device=partVector[k];
        int ghostIndex=globalToGhostIndex[device][k];
        //
        // if(myRank==1 && i==5150){
        //     printf("Global nbh %d \n",k);
        // }
        rho = receiveBuffer[device][ghostIndex].prim[0];
        u1 = receiveBuffer[device][ghostIndex].prim[1];
        u2 = receiveBuffer[device][ghostIndex].prim[2];
        u3 = receiveBuffer[device][ghostIndex].prim[3];
        pr = receiveBuffer[device][ghostIndex].prim[4];
        //
        // if(myRank==1 && i==5150){
        //     printf("Global nbh %d \n",k);
        //     printf("rho %f \n",rho);
        //     printf("u1 %f \n",u1);
        //     printf("u2 %f \n",u2);
        //     printf("u3 %f \n",u3);
        //     printf("pr %f \n",pr);
        // }
        //
        delx = receiveBuffer[device][ghostIndex].x - splitPoint[i].x;
        dely = receiveBuffer[device][ghostIndex].y - splitPoint[i].y;
        delz = receiveBuffer[device][ghostIndex].z - splitPoint[i].z;
        //
        dist = sqrt(delx * delx + dely * dely + delz * delz);
        //
        mod_u = sqrt(u1 * u1 + u2 * u2 + u3 * u3);

        delt = CFL * dist / (mod_u + 3.00 * sqrt(pr / rho));
        //
        if (delt < min_delt)
        {
            min_delt = delt;
        }
        //
    }
    splitPoint[i].delt = min_delt;
    //
    // if(myRank==1 && i==5150){
    //     printf("Global 1 Rank %d\n",splitPoint[i].globalIndex);
    //     printf("Timestep delt %f \n",min_delt);
    //     printf("splitPoint[%d].delt=%.15f\n",i,splitPoint[i].delt);
    //     printf("splitPoint[%d].min_dist=%.15f\n",i,splitPoint[i].min_dist);
    //     printf("Status =%d\n",splitPoint[i].status);
    // }
    // if(splitPoint[i].isGhost){
	// 	for(int t=0;t<splitPoint[i].numberOfPartitionsToSendTo;++t){
	// 		sendBuffer[splitPoint[i].partitions[t]][splitPoint[i].ghostIndex[t]].delt=splitPoint[i].delt;
	// 	}
	// }
    
}
