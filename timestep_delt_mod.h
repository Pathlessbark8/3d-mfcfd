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
        rho = splitPoint[k].prim[0];
        u1 = splitPoint[k].prim[1];
        u2 = splitPoint[k].prim[2];
        u3 = splitPoint[k].prim[3];
        pr = splitPoint[k].prim[4];
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
        rho = receiveBuffer[device][ghostIndex].prim[0];
        u1 = receiveBuffer[device][ghostIndex].prim[1];
        u2 = receiveBuffer[device][ghostIndex].prim[2];
        u3 = receiveBuffer[device][ghostIndex].prim[3];
        pr = receiveBuffer[device][ghostIndex].prim[4];
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
    // if(splitPoint[i].globalIndex== 599001){
    //     for(int r=0;r<5;r++){
    //         printf("q[%d] is :%.15f \n",i,splitPoint[i].q[r]);
    //         printf("dq[0][%d] is :%.15f \n",i,splitPoint[i].dq[0][r]);
    //         printf("dq[1][%d] is :%.15f \n",i,splitPoint[i].dq[1][r]);
    //         printf("dq[2][%d] is :%.15f \n",i,splitPoint[i].dq[2][r]);
    //     }
    // }
    
}
