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
//
//	First written on 23.06.2021
//
//
#pragma once

#include "data_structure_mod.h"
//
//

//
void get_interior_point_neighbours(int i)
//
{
    //
    int k, r, count;
    double xi, yi, zi, xk, yk, zk;
    double delx, dely, delz;
    double dels, delt, deln;
    //
    xi = point.x[i];
    yi = point.y[i];
    zi = point.z[i];
    //
    point.xpos_nbhs[i] = 0;
    point.xneg_nbhs[i] = 0;
    point.ypos_nbhs[i] = 0;
    point.yneg_nbhs[i] = 0;
    point.zpos_nbhs[i] = 0;
    point.zneg_nbhs[i] = 0;
    //
    for (r = 0; r < point.nbhs[i]; r++)
    {
        k = point.conn[r][i];
        xk = point.x[k];
        yk = point.y[k];
        zk = point.z[k];
        //
        delx = xk - xi;
        dely = yk - yi;
        delz = zk - zi;
        //
        dels = delx * point.tan1[0][i] + dely * point.tan1[1][i] + delz * point.tan1[2][i];
        delt = delx * point.tan2[0][i] + dely * point.tan2[1][i] + delz * point.tan2[2][i];
        deln = delx * point.nor[0][i] + dely * point.nor[1][i] + delz * point.nor[2][i];
        //
        if (dels <= 0.00)
        {
            point.xpos_nbhs[i] = point.xpos_nbhs[i] + 1;
            count = point.xpos_nbhs[i];
            point.xpos_conn[count-1][i] = k;
        }
        //
        if (dels > 0.00)
        {
            point.xneg_nbhs[i] = point.xneg_nbhs[i] + 1;
            count = point.xneg_nbhs[i];
            point.xneg_conn[count-1][i] = k;
        }
        //
        if (delt <= 0.00)
        {
            point.ypos_nbhs[i] = point.ypos_nbhs[i] + 1;
            count = point.ypos_nbhs[i];
            point.ypos_conn[count-1][i] = k;
        }
        //
        if (delt > 0.00)
        {
            point.yneg_nbhs[i] = point.yneg_nbhs[i] + 1;
            count = point.yneg_nbhs[i];
            point.yneg_conn[count-1][i] = k;
        }
        //
        if (deln <= 0.00)
        {
            point.zpos_nbhs[i] = point.zpos_nbhs[i] + 1;
            count = point.zpos_nbhs[i];
            point.zpos_conn[count-1][i] = k;
        }
        //
        if (deln > 0.00)
        {
            point.zneg_nbhs[i] = point.zneg_nbhs[i] + 1;
            count = point.zneg_nbhs[i];
            point.zneg_conn[count-1][i] = k;
        }
        //
    }
    //
    //
}
//
// void get_interior_point_neighbours_local_points(int i){
//     //
 
//         //
//     }
// }
//
//
void get_wall_point_neighbours(int i)
//
{
    //
    int k, r, count;
    double xi, yi, zi, xk, yk, zk;
    double delx, dely, delz;
    double dels, delt, deln;
    //
    xi = point.x[i];
    yi = point.y[i];
    zi = point.z[i];
    //
    point.xpos_nbhs[i] = 0;
    point.xneg_nbhs[i] = 0;
    point.ypos_nbhs[i] = 0;
    point.yneg_nbhs[i] = 0;
    point.zneg_nbhs[i] = 0;
    //
    for (r = 0; r < point.nbhs[i]; r++)
    {
        k = point.conn[r][i];
        xk = point.x[k];
        yk = point.y[k];
        zk = point.z[k];
        //
        delx = xk - xi;
        dely = yk - yi;
        delz = zk - zi;
        //
        dels = delx * point.tan1[0][i] + dely * point.tan1[1][i] + delz * point.tan1[2][i];
        delt = delx * point.tan2[0][i] + dely * point.tan2[1][i] + delz * point.tan2[2][i];
        deln = delx * point.nor[0][i] + dely * point.nor[1][i] + delz * point.nor[2][i];
        //
        //
        if (dels <= 0.00)
        {
            point.xpos_nbhs[i] = point.xpos_nbhs[i] + 1;
            count = point.xpos_nbhs[i];
            point.xpos_conn[count-1][i] = k;
        }
        //
        if (dels >= 0.00)
        {
            point.xneg_nbhs[i] = point.xneg_nbhs[i] + 1;
            count = point.xneg_nbhs[i];
            point.xneg_conn[count-1][i] = k;
        }
        //
        if (delt <= 0.00)
        {
            point.ypos_nbhs[i] = point.ypos_nbhs[i] + 1;
            count = point.ypos_nbhs[i];
            point.ypos_conn[count-1][i] = k;
        }
        //
        if (delt >= 0.00)
        {
            point.yneg_nbhs[i] = point.yneg_nbhs[i] + 1;
            count = point.yneg_nbhs[i];
            point.yneg_conn[count-1][i] = k;
        }
        //
        point.zneg_nbhs[i] = point.zneg_nbhs[i] + 1;
        count = point.zneg_nbhs[i];
        point.zneg_conn[count-1][i] = k;
    }
    //
    //
}
//
//
//
void get_outer_point_neighbours(int i)
//
//
{
    //
    int k, r, count;
    double xi, yi, zi, xk, yk, zk;
    double delx, dely, delz;
    double dels, delt, deln;
    //
    xi = point.x[i];
    yi = point.y[i];
    zi = point.z[i];
    //
    point.xpos_nbhs[i] = 0;
    point.xneg_nbhs[i] = 0;
    point.ypos_nbhs[i] = 0;
    point.yneg_nbhs[i] = 0;
    point.zneg_nbhs[i] = 0;
    //
    for (r = 0; r < point.nbhs[i]; r++)
    {
        k = point.conn[r][i];
        xk = point.x[k];
        yk = point.y[k];
        zk = point.z[k];
        //
        delx = xk - xi;
        dely = yk - yi;
        delz = zk - zi;
        //
        dels = delx * point.tan1[0][i] + dely * point.tan1[1][i] + delz * point.tan1[2][i];
        delt = delx * point.tan2[0][i] + dely * point.tan2[1][i] + delz * point.tan2[2][i];
        deln = delx * point.nor[0][i] + dely * point.nor[1][i] + delz * point.nor[2][i];
        //
        if (dels <= 0.00)
        {
            point.xpos_nbhs[i] = point.xpos_nbhs[i] + 1;
            count = point.xpos_nbhs[i];
            point.xpos_conn[count-1][i] = k;
        }
        //
        if (dels >= 0.00)
        {
            point.xneg_nbhs[i] = point.xneg_nbhs[i] + 1;
            count = point.xneg_nbhs[i];
            point.xneg_conn[count-1][i] = k;
        }
        //
        if (delt <= 0.00)
        {
            point.ypos_nbhs[i] = point.ypos_nbhs[i] + 1;
            count = point.ypos_nbhs[i];
            point.ypos_conn[count-1][i] = k;
        }
        //
        if (delt >= 0.00)
        {
            point.yneg_nbhs[i] = point.yneg_nbhs[i] + 1;
            count = point.yneg_nbhs[i];
            point.yneg_conn[count-1][i] = k;
        }
        //
        point.zpos_nbhs[i] = point.zpos_nbhs[i] + 1;
        count = point.zpos_nbhs[i];
        point.zpos_conn[count-1][i] = k;
    }
    //
    //
}
//
//

//
void generate_split_stencils()
//
{
    //
    int i, k;

    for (k = 0; k < interior_points; k++)
    {
        i = interior_points_index[k];
        get_interior_point_neighbours(i);
    }
    //
    for (k = 0; k < wall_points; k++)
    {
        i = wall_points_index[k];
        get_wall_point_neighbours(i);
    }
    //
    for (k = 0; k < outer_points; k++)
    {
        i = outer_points_index[k];
        get_outer_point_neighbours(i);
    }
    //
}
//
//
// FOR MULTI NODE
//
__global__ void generate_split_stencils_interior_multi_nccl(splitPoints *splitPoint, int *interiorPointsLocalIndex, int interiorPointsLocal,int*partVector,int *globalToLocalIndex,int **globalToGhostIndex,transferPoints ** receiveBuffer)
//
{
    //
    int i, k;
    int bx = blockIdx.x;
    int tx = threadIdx.x;
    int t = bx * blockDim.x + tx;

	if(t<0 || t>=interiorPointsLocal){
		return;
	}
    i = interiorPointsLocalIndex[t];
    //
    int r, count;
    double xi, yi, zi, xk, yk, zk;
    double delx, dely, delz;
    double dels, delt, deln;
    //
    xi = splitPoint[i].x;
    yi = splitPoint[i].y;
    zi = splitPoint[i].z;
    // //
    splitPoint[i].numberOfLocalxposNbhs = 0;
    splitPoint[i].numberOfLocalxnegNbhs = 0;
    splitPoint[i].numberOfLocalyposNbhs = 0;
    splitPoint[i].numberOfLocalynegNbhs = 0;
    splitPoint[i].numberOfLocalzposNbhs = 0;
    splitPoint[i].numberOfLocalznegNbhs = 0;
    //
    splitPoint[i].numberOfGhostxposNbhs = 0;
    splitPoint[i].numberOfGhostxnegNbhs = 0;
    splitPoint[i].numberOfGhostyposNbhs = 0;
    splitPoint[i].numberOfGhostynegNbhs = 0;
    splitPoint[i].numberOfGhostzposNbhs = 0;
    splitPoint[i].numberOfGhostznegNbhs = 0;
    
    for (r = 0; r < splitPoint[i].numberOfLocalNbhs; r++)
    {
        k = splitPoint[i].localNbhs[r];
        int nbh=globalToLocalIndex[k];
        xk = splitPoint[nbh].x;
        yk = splitPoint[nbh].y;
        zk = splitPoint[nbh].z;
        //
        delx = xk - xi;
        dely = yk - yi;
        delz = zk - zi;
        //
        dels = delx * splitPoint[i].tan1[0] + dely * splitPoint[i].tan1[1] + delz * splitPoint[i].tan1[2];
        delt = delx * splitPoint[i].tan2[0] + dely * splitPoint[i].tan2[1] + delz * splitPoint[i].tan2[2];
        deln = delx * splitPoint[i].nor[0] + dely * splitPoint[i].nor[1] + delz * splitPoint[i].nor[2];
            
        if (dels <= 0.00)
        {
            splitPoint[i].numberOfLocalxposNbhs = splitPoint[i].numberOfLocalxposNbhs + 1;
            count = splitPoint[i].numberOfLocalxposNbhs;
            splitPoint[i].localxpos_conn[count-1] = k;
        }
        //
        if (dels > 0.00)
        {
            splitPoint[i].numberOfLocalxnegNbhs = splitPoint[i].numberOfLocalxnegNbhs + 1;
            count = splitPoint[i].numberOfLocalxnegNbhs;
            splitPoint[i].localxneg_conn[count-1] = k;
        }
        //
        if (delt <= 0.00)
        {
            splitPoint[i].numberOfLocalyposNbhs = splitPoint[i].numberOfLocalyposNbhs + 1;
            count = splitPoint[i].numberOfLocalyposNbhs;
            splitPoint[i].localypos_conn[count-1] = k;
        }
        //
        if (delt > 0.00)
        {
            splitPoint[i].numberOfLocalynegNbhs = splitPoint[i].numberOfLocalynegNbhs + 1;
            count = splitPoint[i].numberOfLocalynegNbhs;
            splitPoint[i].localyneg_conn[count-1] = k;
        }
        //
        if (deln <= 0.00)
        {
            splitPoint[i].numberOfLocalzposNbhs = splitPoint[i].numberOfLocalzposNbhs + 1;
            count = splitPoint[i].numberOfLocalzposNbhs;
            splitPoint[i].localzpos_conn[count-1] = k;
        }
        //
        if (deln > 0.00)
        {
            splitPoint[i].numberOfLocalznegNbhs = splitPoint[i].numberOfLocalznegNbhs + 1;
            count = splitPoint[i].numberOfLocalznegNbhs;
            splitPoint[i].localzneg_conn[count-1] = k;
        }
    }
    //
    for (r = 0; r < splitPoint[i].numberOfGhostNbhs; r++)
    {
        k = splitPoint[i].ghostNbhs[r];
        int device=partVector[k];
        int ghostIndex=globalToGhostIndex[device][k];

        xk = receiveBuffer[device][ghostIndex].x;
        yk = receiveBuffer[device][ghostIndex].y;
        zk = receiveBuffer[device][ghostIndex].z;
        //
        delx = xk - xi;
        dely = yk - yi;
        delz = zk - zi;
        //
        dels = delx * splitPoint[i].tan1[0] + dely * splitPoint[i].tan1[1] + delz * splitPoint[i].tan1[2];
        delt = delx * splitPoint[i].tan2[0] + dely * splitPoint[i].tan2[1] + delz * splitPoint[i].tan2[2];
        deln = delx * splitPoint[i].nor[0] + dely * splitPoint[i].nor[1] + delz * splitPoint[i].nor[2];
            
        if (dels <= 0.00)
        {
            splitPoint[i].numberOfGhostxposNbhs = splitPoint[i].numberOfGhostxposNbhs + 1;
            count = splitPoint[i].numberOfGhostxposNbhs;
            splitPoint[i].ghostxpos_conn[count-1] = k;
            // if(splitPoint[i].globalIndex==529172){
            //     printf("Ghost xpos_conn[%d]=%d\n",count-1,k);
            // }
        }
        //
        if (dels > 0.00)
        {
            splitPoint[i].numberOfGhostxnegNbhs = splitPoint[i].numberOfGhostxnegNbhs + 1;
            count = splitPoint[i].numberOfGhostxnegNbhs;
            splitPoint[i].ghostxneg_conn[count-1] = k;
        }
        //
        if (delt <= 0.00)
        {
            splitPoint[i].numberOfGhostyposNbhs = splitPoint[i].numberOfGhostyposNbhs + 1;
            count = splitPoint[i].numberOfGhostyposNbhs;
            splitPoint[i].ghostypos_conn[count-1] = k;
        }
        //
        if (delt > 0.00)
        {
            splitPoint[i].numberOfGhostynegNbhs = splitPoint[i].numberOfGhostynegNbhs + 1;
            count = splitPoint[i].numberOfGhostynegNbhs;
            splitPoint[i].ghostyneg_conn[count-1] = k;
        }
        //
        if (deln <= 0.00)
        {
            splitPoint[i].numberOfGhostzposNbhs = splitPoint[i].numberOfGhostzposNbhs + 1;
            count = splitPoint[i].numberOfGhostzposNbhs;
            splitPoint[i].ghostzpos_conn[count-1] = k;
        }
        //
        if (deln > 0.00)
        {
            splitPoint[i].numberOfGhostznegNbhs = splitPoint[i].numberOfGhostznegNbhs + 1;
            count = splitPoint[i].numberOfGhostznegNbhs;
            splitPoint[i].ghostzneg_conn[count-1] = k;
        }
    }
}

__global__ void generate_split_stencils_wall_multi_nccl(int myRank,splitPoints *splitPoint, int *wallPointsLocalIndex, int wallPointsLocal,int*partVector,int *globalToLocalIndex,int **globalToGhostIndex,transferPoints ** receiveBuffer)
//
{
    //
    int i, k;
    int bx = blockIdx.x;
    int tx = threadIdx.x;
    int t = bx * blockDim.x + tx;

	if(t<0 || t>=wallPointsLocal){
		return;
	}
    i = wallPointsLocalIndex[t];
    //
    int r, count;
    double xi, yi, zi, xk, yk, zk;
    double delx, dely, delz;
    double dels, delt, deln;
    //
    xi = splitPoint[i].x;
    yi = splitPoint[i].y;
    zi = splitPoint[i].z;
    // //
    splitPoint[i].numberOfLocalxposNbhs = 0;
    splitPoint[i].numberOfLocalxnegNbhs = 0;
    splitPoint[i].numberOfLocalyposNbhs = 0;
    splitPoint[i].numberOfLocalynegNbhs = 0;
    splitPoint[i].numberOfLocalzposNbhs = 0;
    splitPoint[i].numberOfLocalznegNbhs = 0;
    //
    splitPoint[i].numberOfGhostxposNbhs = 0;
    splitPoint[i].numberOfGhostxnegNbhs = 0;
    splitPoint[i].numberOfGhostyposNbhs = 0;
    splitPoint[i].numberOfGhostynegNbhs = 0;
    splitPoint[i].numberOfGhostzposNbhs = 0;
    splitPoint[i].numberOfGhostznegNbhs = 0;
    
    for (r = 0; r < splitPoint[i].numberOfLocalNbhs; r++)
    {
        k = splitPoint[i].localNbhs[r];
        int nbh=globalToLocalIndex[k];
        xk = splitPoint[nbh].x;
        yk = splitPoint[nbh].y;
        zk = splitPoint[nbh].z;
        //
        delx = xk - xi;
        dely = yk - yi;
        delz = zk - zi;
        //
        dels = delx * splitPoint[i].tan1[0] + dely * splitPoint[i].tan1[1] + delz * splitPoint[i].tan1[2];
        delt = delx * splitPoint[i].tan2[0] + dely * splitPoint[i].tan2[1] + delz * splitPoint[i].tan2[2];
        deln = delx * splitPoint[i].nor[0] + dely * splitPoint[i].nor[1] + delz * splitPoint[i].nor[2];
            
        if (dels <= 0.00)
        {
            
            splitPoint[i].numberOfLocalxposNbhs = splitPoint[i].numberOfLocalxposNbhs + 1;
            count = splitPoint[i].numberOfLocalxposNbhs;
            splitPoint[i].localxpos_conn[count-1] = k;
        }
        //
        if (dels >= 0.00)
        {
            splitPoint[i].numberOfLocalxnegNbhs = splitPoint[i].numberOfLocalxnegNbhs + 1;
            count = splitPoint[i].numberOfLocalxnegNbhs;
            splitPoint[i].localxneg_conn[count-1] = k;
        }
        //
        if (delt <= 0.00)
        {
            splitPoint[i].numberOfLocalyposNbhs = splitPoint[i].numberOfLocalyposNbhs + 1;
            count = splitPoint[i].numberOfLocalyposNbhs;
            splitPoint[i].localypos_conn[count-1] = k;
        }
        //
        if (delt >= 0.00)
        {
            splitPoint[i].numberOfLocalynegNbhs = splitPoint[i].numberOfLocalynegNbhs + 1;
            count = splitPoint[i].numberOfLocalynegNbhs;
            splitPoint[i].localyneg_conn[count-1] = k;
        }
        //
        splitPoint[i].numberOfLocalznegNbhs = splitPoint[i].numberOfLocalznegNbhs + 1;
        count = splitPoint[i].numberOfLocalznegNbhs;
        splitPoint[i].localzneg_conn[count-1] = k;
    }
    //
    for (r = 0; r < splitPoint[i].numberOfGhostNbhs; r++)
    {
        k = splitPoint[i].ghostNbhs[r];
        // if(myRank==1 && i==6081)
        // printf("k=%d\n",k);
        int device=partVector[k];
        int ghostIndex=globalToGhostIndex[device][k];

        xk = receiveBuffer[device][ghostIndex].x;
        yk = receiveBuffer[device][ghostIndex].y;
        zk = receiveBuffer[device][ghostIndex].z;
        //
        delx = xk - xi;
        dely = yk - yi;
        delz = zk - zi;
        //
        dels = delx * splitPoint[i].tan1[0] + dely * splitPoint[i].tan1[1] + delz * splitPoint[i].tan1[2];
        delt = delx * splitPoint[i].tan2[0] + dely * splitPoint[i].tan2[1] + delz * splitPoint[i].tan2[2];
        deln = delx * splitPoint[i].nor[0] + dely * splitPoint[i].nor[1] + delz * splitPoint[i].nor[2];
            
        if (dels <= 0.00)
        {
            splitPoint[i].numberOfGhostxposNbhs = splitPoint[i].numberOfGhostxposNbhs + 1;
            count = splitPoint[i].numberOfGhostxposNbhs;
            splitPoint[i].ghostxpos_conn[count-1] = k;
        }
        //
        if (dels > 0.00)
        {
            splitPoint[i].numberOfGhostxnegNbhs = splitPoint[i].numberOfGhostxnegNbhs + 1;
            count = splitPoint[i].numberOfGhostxnegNbhs;
            splitPoint[i].ghostxneg_conn[count-1] = k;
        }
        //
        if (delt <= 0.00)
        {
            splitPoint[i].numberOfGhostyposNbhs = splitPoint[i].numberOfGhostyposNbhs + 1;
            count = splitPoint[i].numberOfGhostyposNbhs;
            splitPoint[i].ghostypos_conn[count-1] = k;
        }
        //
        if (delt > 0.00)
        {
            splitPoint[i].numberOfGhostynegNbhs = splitPoint[i].numberOfGhostynegNbhs + 1;
            count = splitPoint[i].numberOfGhostynegNbhs;
            splitPoint[i].ghostyneg_conn[count-1] = k;
        }
        //
        if (deln <= 0.00)
        {
            splitPoint[i].numberOfGhostzposNbhs = splitPoint[i].numberOfGhostzposNbhs + 1;
            count = splitPoint[i].numberOfGhostzposNbhs;
            splitPoint[i].ghostzpos_conn[count-1] = k;
        }
        //
        splitPoint[i].numberOfGhostznegNbhs = splitPoint[i].numberOfGhostznegNbhs + 1;
        count = splitPoint[i].numberOfGhostznegNbhs;
        splitPoint[i].ghostzneg_conn[count-1] = k;
    }
}

__global__ void generate_split_stencils_outer_multi_nccl(splitPoints *splitPoint, int *outerPointsLocalIndex, int outerPointsLocal,int*partVector,int *globalToLocalIndex,int **globalToGhostIndex,transferPoints ** receiveBuffer)
//
{
    //
    int i, k;
    int bx = blockIdx.x;
    int tx = threadIdx.x;
    int t = bx * blockDim.x + tx;

	if(t<0 || t>=outerPointsLocal){
		return;
	}
    i = outerPointsLocalIndex[t];
    //
    int r, count;
    double xi, yi, zi, xk, yk, zk;
    double delx, dely, delz;
    double dels, delt, deln;
    //
    xi = splitPoint[i].x;
    yi = splitPoint[i].y;
    zi = splitPoint[i].z;
    // //
    splitPoint[i].numberOfLocalxposNbhs = 0;
    splitPoint[i].numberOfLocalxnegNbhs = 0;
    splitPoint[i].numberOfLocalyposNbhs = 0;
    splitPoint[i].numberOfLocalynegNbhs = 0;
    splitPoint[i].numberOfLocalzposNbhs = 0;
    splitPoint[i].numberOfLocalznegNbhs = 0;
    //
    splitPoint[i].numberOfGhostxposNbhs = 0;
    splitPoint[i].numberOfGhostxnegNbhs = 0;
    splitPoint[i].numberOfGhostyposNbhs = 0;
    splitPoint[i].numberOfGhostynegNbhs = 0;
    splitPoint[i].numberOfGhostzposNbhs = 0;
    splitPoint[i].numberOfGhostznegNbhs = 0;
    
    for (r = 0; r < splitPoint[i].numberOfLocalNbhs; r++)
    {
        k = splitPoint[i].localNbhs[r];
        int nbh=globalToLocalIndex[k];
        xk = splitPoint[nbh].x;
        yk = splitPoint[nbh].y;
        zk = splitPoint[nbh].z;
        //
        delx = xk - xi;
        dely = yk - yi;
        delz = zk - zi;
        //
        dels = delx * splitPoint[i].tan1[0] + dely * splitPoint[i].tan1[1] + delz * splitPoint[i].tan1[2];
        delt = delx * splitPoint[i].tan2[0] + dely * splitPoint[i].tan2[1] + delz * splitPoint[i].tan2[2];
        deln = delx * splitPoint[i].nor[0] + dely * splitPoint[i].nor[1] + delz * splitPoint[i].nor[2];
            
        if (dels <= 0.00)
        {
            splitPoint[i].numberOfLocalxposNbhs = splitPoint[i].numberOfLocalxposNbhs + 1;
            count = splitPoint[i].numberOfLocalxposNbhs;
            splitPoint[i].localxpos_conn[count-1] = k;
        }
        //
        if (dels >= 0.00)
        {
            splitPoint[i].numberOfLocalxnegNbhs = splitPoint[i].numberOfLocalxnegNbhs + 1;
            count = splitPoint[i].numberOfLocalxnegNbhs;
            splitPoint[i].localxneg_conn[count-1] = k;
        }
        //
        if (delt <= 0.00)
        {
            splitPoint[i].numberOfLocalyposNbhs = splitPoint[i].numberOfLocalyposNbhs + 1;
            count = splitPoint[i].numberOfLocalyposNbhs;
            splitPoint[i].localypos_conn[count-1] = k;
        }
        //
        if (delt >= 0.00)
        {
            splitPoint[i].numberOfLocalynegNbhs = splitPoint[i].numberOfLocalynegNbhs + 1;
            count = splitPoint[i].numberOfLocalynegNbhs;
            splitPoint[i].localyneg_conn[count-1] = k;
        }
        //
        splitPoint[i].numberOfLocalzposNbhs = splitPoint[i].numberOfLocalzposNbhs + 1;
        count = splitPoint[i].numberOfLocalzposNbhs;
        splitPoint[i].localzpos_conn[count-1] = k;
    }

    //
    for (r = 0; r < splitPoint[i].numberOfGhostNbhs; r++)
    {
        k = splitPoint[i].ghostNbhs[r];
        int device=partVector[k];
        int ghostIndex=globalToGhostIndex[device][k];

        xk = receiveBuffer[device][ghostIndex].x;
        yk = receiveBuffer[device][ghostIndex].y;
        zk = receiveBuffer[device][ghostIndex].z;
        //
        delx = xk - xi;
        dely = yk - yi;
        delz = zk - zi;
        //
        dels = delx * splitPoint[i].tan1[0] + dely * splitPoint[i].tan1[1] + delz * splitPoint[i].tan1[2];
        delt = delx * splitPoint[i].tan2[0] + dely * splitPoint[i].tan2[1] + delz * splitPoint[i].tan2[2];
        deln = delx * splitPoint[i].nor[0] + dely * splitPoint[i].nor[1] + delz * splitPoint[i].nor[2];
            
       if (dels <= 0.00)
        {
            splitPoint[i].numberOfGhostxposNbhs = splitPoint[i].numberOfGhostxposNbhs + 1;
            count = splitPoint[i].numberOfGhostxposNbhs;
            splitPoint[i].ghostxpos_conn[count-1] = k;
        }
        //
        if (dels >= 0.00)
        {
            splitPoint[i].numberOfGhostxnegNbhs = splitPoint[i].numberOfGhostxnegNbhs + 1;
            count = splitPoint[i].numberOfGhostxnegNbhs;
            splitPoint[i].ghostxneg_conn[count-1] = k;
        }
        //
        if (delt <= 0.00)
        {
            splitPoint[i].numberOfGhostyposNbhs = splitPoint[i].numberOfGhostyposNbhs + 1;
            count = splitPoint[i].numberOfGhostyposNbhs;
            splitPoint[i].ghostypos_conn[count-1] = k;
        }
        //
        if (delt >= 0.00)
        {
            splitPoint[i].numberOfGhostynegNbhs = splitPoint[i].numberOfGhostynegNbhs + 1;
            count = splitPoint[i].numberOfGhostynegNbhs;
            splitPoint[i].ghostyneg_conn[count-1] = k;
        }
        //
        splitPoint[i].numberOfGhostzposNbhs = splitPoint[i].numberOfGhostzposNbhs + 1;
        count = splitPoint[i].numberOfGhostzposNbhs;
        splitPoint[i].ghostzpos_conn[count-1] = k;

    }
}

