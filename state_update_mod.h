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
//
//
#include "data_structure_mod.h"
#include "primitive_to_conserved_vector_mod.h"
#include "conserved_to_primitive_vector_mod.h"
#include "conserved_vector_Ubar_mod.h"
#include <float.h>

//
//
void state_update()
{
	//
	int i, k, nbh, p, r;
	double U[5], temp, min_dist;
	double res_sqr, sum_res_sqr;
	double dx, dy, dz, ds;
	//
	sum_res_sqr = 0.00;
	max_res = 0.00;
	//
	//
	for (i = 0; i < wall_points; i++)
	{
		//
		k = wall_points_index[i];
		primitive_to_conserved(k, U);
		temp = U[0];
		for (int r = 0; r < 5; r++)
		{
			U[r] = U[r] - point.flux_res[r][k];
		}
		U[3] = 0.00;
		//
		res_sqr = (U[0] - temp) * (U[0] - temp);
		if (res_sqr > max_res)
		{
			max_res = res_sqr;
			max_res_point = k;
		}
		sum_res_sqr = sum_res_sqr + res_sqr;
		//
		//                                        print*, i, k, point.flux_res[r][k]
		conserved_to_primitive(k, U);
		//
	}
	//
	for (i = 0; i < outer_points; i++)
	{
		//
		k = outer_points_index[i];
		conserved_vector_Ubar(k, U);
		temp = U[0];
		for (int r = 0; r < 5; r++)
		{
			U[r] = U[r] - point.flux_res[r][k];
		}
		//

		conserved_to_primitive(k, U);
		//
	}
	//
	for (i = 0; i < interior_points; i++)
	{
		//
		k = interior_points_index[i];
		primitive_to_conserved(k, U);
		temp = U[0];
		for (int r = 0; r < 5; r++)
		{
			U[r] = U[r] - point.flux_res[r][k];
		}
		//
		res_sqr = (U[0] - temp) * (U[0] - temp);
		if (res_sqr > max_res)
		{
			max_res = res_sqr;
			max_res_point = k;
		}
		sum_res_sqr = sum_res_sqr + res_sqr;
		//
		conserved_to_primitive(k, U);
		//
	}
	//
	for (i = 0; i < supersonic_outlet_points; i++)
	{
		min_dist = 100000.00;
		k = supersonic_outlet_points_index[i];
		for (int r = 0; r < point.nbhs[k]; r++)
		{
			nbh = point.conn[r][k];
			dx = point.x[nbh] - point.x[k];
			dy = point.y[nbh] - point.y[k];
			dz = point.z[nbh] - point.z[k];
			ds = sqrt(dx * dx + dy * dy + dz * dz);
			if (ds < min_dist && point.status[nbh] != 1)
			{
				min_dist = ds;
				p = nbh;
			}
		}
		for (int r = 0; r < 5; r++)
		{
			point.prim[r][k] = point.prim[r][p];
		}
	}
	//
	for (i = 0; i < supersonic_inlet_points; i++)
	{
		min_dist = 100000.00;
		k = supersonic_inlet_points_index[i];
		for (r = 0; r < point.nbhs[k]; r++)
		{
			nbh = point.conn[r][k];
			dx = point.x[nbh] - point.x[k];
			dy = point.y[nbh] - point.y[k];
			dz = point.z[nbh] - point.z[k];
			ds = sqrt(dx * dx + dy * dy + dz * dz);
			if (ds < min_dist && point.status[nbh] != 1)
			{
				min_dist = ds;
				p = nbh;
			}
		}
		for (int r = 0; r < 5; r++)
		{
			point.prim[r][k] = point.prim[r][p];
		}
	}
	res_new = sqrt(sum_res_sqr) / max_points;
	//
	//
}

__global__ void state_update_cuda(points &point, int wall_points, int outer_points, int interior_points, int supersonic_outlet_points, int supersonic_inlet_points, int *wall_points_index, int *outer_points_index, int *interior_points_index, int *supersonic_outlet_points_index, int *supersonic_inlet_points_index, double *sum_res_sqr, double u1_inf, double u2_inf, double u3_inf, double rho_inf, double pi, double pr_inf)
//
//
{
	//
	int k, nbh, p, r;
	double U[5], temp, min_dist;
	double res_sqr;
	double dx, dy, dz, ds;
	//
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < 0 || i >= max_points)
	{
		return;
	}
	//
	if (i < wall_points)
	{
		//
		k = wall_points_index[i];
		primitive_to_conserved_cuda(point, k, U);
		temp = U[0];
		for (int r = 0; r < 5; r++)
		{
			U[r] = U[r] - point.flux_res[r][k];
		}
		U[3] = 0.00;
		//
		res_sqr = (U[0] - temp) * (U[0] - temp);
		sum_res_sqr[k] = res_sqr;
		//                                        print*, i, k, point.flux_res[r][k]
		conserved_to_primitive_cuda(point, k, U);
		//
	}
	//
	if (i < outer_points)
	{
		//
		k = outer_points_index[i];
		conserved_vector_Ubar_cuda(point, k, U, u1_inf, u2_inf, u3_inf, rho_inf, pi, pr_inf);
		temp = U[0];
		for (int r = 0; r < 5; r++)
		{
			U[r] = U[r] - point.flux_res[r][k];
		}
		//

		conserved_to_primitive_cuda(point, k, U);
		//
	}
	//
	if (i < interior_points)
	{
		//
		k = interior_points_index[i];
		primitive_to_conserved_cuda(point, k, U);
		temp = U[0];
		for (int r = 0; r < 5; r++)
		{
			U[r] = U[r] - point.flux_res[r][k];
		}
		//
		res_sqr = (U[0] - temp) * (U[0] - temp);
		sum_res_sqr[k] = res_sqr;
		//
		conserved_to_primitive_cuda(point, k, U);
		//
	}
	//
	if (i < supersonic_outlet_points)
	{
		min_dist = 100000.00;
		k = supersonic_outlet_points_index[i];
		for (int r = 0; r < point.nbhs[k]; r++)
		{
			nbh = point.conn[r][k];
			dx = point.x[nbh] - point.x[k];
			dy = point.y[nbh] - point.y[k];
			dz = point.z[nbh] - point.z[k];
			ds = sqrt(dx * dx + dy * dy + dz * dz);
			if (ds < min_dist && point.status[nbh] != 1)
			{
				min_dist = ds;
				p = nbh;
			}
		}
		for (int r = 0; r < 5; r++)
		{
			point.prim[r][k] = point.prim[r][p];
		}
	}
	//
	if (i < supersonic_inlet_points)
	{
		min_dist = 100000.00;
		k = supersonic_inlet_points_index[i];
		for (r = 0; r < point.nbhs[k]; r++)
		{
			nbh = point.conn[r][k];
			dx = point.x[nbh] - point.x[k];
			dy = point.y[nbh] - point.y[k];
			dz = point.z[nbh] - point.z[k];
			ds = sqrt(dx * dx + dy * dy + dz * dz);
			if (ds < min_dist && point.status[nbh] != 1)
			{
				min_dist = ds;
				p = nbh;
			}
		}
		for (int r = 0; r < 5; r++)
		{
			point.prim[r][k] = point.prim[r][p];
		}
	}

	//
	//
}
//

__global__ void state_update_wall(points &point, int wall_points, int *wall_points_index, double *sum_res_sqr)
{
	int k;
	double U[5], temp;
	double res_sqr;
	//
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < 0 || i >= wall_points)
	{
		return;
	}
	k = wall_points_index[i];
	primitive_to_conserved_cuda(point, k, U);
	temp = U[0];
	for (int r = 0; r < 5; r++)
	{
		U[r] = U[r] - point.flux_res[r][k];
	}
	U[3] = 0.00;
	//
	res_sqr = (U[0] - temp) * (U[0] - temp);
	sum_res_sqr[k] = res_sqr;
	//                                        print*, i, k, point.flux_res[r][k]
	conserved_to_primitive_cuda(point, k, U);
}
//
//
__global__ void state_update_outer(points &point, int outer_points, int *outer_points_index, double u1_inf, double u2_inf, double u3_inf, double rho_inf, double pi, double pr_inf)
{
	int k;
	double U[5];
	//
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < 0 || i >= outer_points)
	{
		return;
	}
	k = outer_points_index[i];
	conserved_vector_Ubar_cuda(point, k, U, u1_inf, u2_inf, u3_inf, rho_inf, pi, pr_inf);
	// temp = U[0];
	for (int r = 0; r < 5; r++)
	{
		U[r] = U[r] - point.flux_res[r][k];
	}
	//

	conserved_to_primitive_cuda(point, k, U);
}
//
//
__global__ void state_update_interior(points &point, int interior_points, int *interior_points_index, double *sum_res_sqr)
{
	int k;
	double U[5], temp;
	double res_sqr;
	//
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < 0 || i >= interior_points)
	{
		return;
	}
	k = interior_points_index[i];
	primitive_to_conserved_cuda(point, k, U);
	temp = U[0];
	for (int r = 0; r < 5; r++)
	{
		U[r] = U[r] - point.flux_res[r][k];
	}
	//
	res_sqr = (U[0] - temp) * (U[0] - temp);
	sum_res_sqr[k] = res_sqr;
	//
	conserved_to_primitive_cuda(point, k, U);
}
//
//
void state_update_supersonic_outlet()
{
	int k, nbh, p ,i;
	double min_dist;
	double dx, dy, dz, ds;
	//
	for (i = 0; i < supersonic_outlet_points; i++)
	{
		min_dist = 100000.00;
		k = supersonic_outlet_points_index[i];
		for (int r = 0; r < point.nbhs[k]; r++)
		{
			nbh = point.conn[r][k];
			dx = point.x[nbh] - point.x[k];
			dy = point.y[nbh] - point.y[k];
			dz = point.z[nbh] - point.z[k];
			ds = sqrt(dx * dx + dy * dy + dz * dz);
			if (ds < min_dist && point.status[nbh] != 1)
			{
				min_dist = ds;
				p = nbh;
			}
		}
		for (int r = 0; r < 5; r++)
		{
			point.prim[r][k] = point.prim[r][p];
		}
	}
}
//
//
void state_update_supersonic_inlet()
{
	int k, nbh, p, r,i;
	double min_dist;
	double dx, dy, dz, ds;
	//
	for (i = 0; i < supersonic_inlet_points; i++)
	{
		min_dist = 100000.00;
		k = supersonic_inlet_points_index[i];
		for (r = 0; r < point.nbhs[k]; r++)
		{
			nbh = point.conn[r][k];
			dx = point.x[nbh] - point.x[k];
			dy = point.y[nbh] - point.y[k];
			dz = point.z[nbh] - point.z[k];
			ds = sqrt(dx * dx + dy * dy + dz * dz);
			if (ds < min_dist && point.status[nbh] != 1)
			{
				min_dist = ds;
				p = nbh;
			}
		}
		for (int r = 0; r < 5; r++)
		{
			point.prim[r][k] = point.prim[r][p];
		}
	}
}


// MULTI GPU FUNCTIONS

__global__ void state_update_wall_multi_nccl(int myRank,splitPoints *splitPoint, int wallPointsLocal, int *wallPointsLocalIndex, double *sum_res_sqr)
{
	int k;
	double U[5], temp;
	double res_sqr;
	//
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < 0 || i >= wallPointsLocal)
	{
		return;
	}
	k = wallPointsLocalIndex[i];
	primitive_to_conserved_multi_nccl(splitPoint, k, U);
	temp = U[0];
	for (int r = 0; r < 5; r++)
	{
		// if(U[r]!=U[r]){
		// 	printf("k is %d\n",k);
		// 	printf("r is %d\n",r);
		// }
		U[r] = U[r] - splitPoint[k].flux_res[r];
		// if(splitPoint[k].prim[0]<0){
		// 	printf("Wall splitPoint[%d].prim[%d] %.15f",splitPoint[k].globalIndex,r,splitPoint[k].prim[0]);
		// }
		// if(k==50001 && ((splitPoint[k].flux_res[0]<0 && 0-splitPoint[k].flux_res[0]>10e-15) ||(splitPoint[k].flux_res[4]<0 && 0-splitPoint[k].flux_res[4]>10e-15)) ){
		// 	printf("Wall Here %d\n",k);
		// }
		// if(splitPoint[k].flux_res[r]!=splitPoint[k].flux_res[r]){
		// 	printf("k is %d\n",k);
		// 	printf("r is %d\n",r);
		// }
	}
	U[3] = 0.00;
	//
	res_sqr = (U[0] - temp) * (U[0] - temp);

	sum_res_sqr[k] = res_sqr;
	//                                        print*, i, k, point.flux_res[r][k]
	conserved_to_primitive_multi_nccl(splitPoint, k, U);
}

__global__ void state_update_interior_multi_nccl(splitPoints *splitPoint, int interiorPointsLocal, int *interiorPointsLocalIndex, double *sum_res_sqr)
{
	int k;
	double U[5], temp;
	double res_sqr;
	//
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < 0 || i >= interiorPointsLocal)
	{
		return;
	}
	k = interiorPointsLocalIndex[i];
	primitive_to_conserved_multi_nccl(splitPoint, k, U);
	temp = U[0];
	for (int r = 0; r < 5; r++)
	{
		U[r] = U[r] - splitPoint[k].flux_res[r];
		// if(splitPoint[k].prim[0]<0){
		// 	printf("Interior splitPoint[%d].prim[%d] %.15f",splitPoint[k].globalIndex,r,splitPoint[k].prim[0]);
		// }
		// if(k==50001 && ((splitPoint[k].flux_res[0]<0 && 0-splitPoint[k].flux_res[0]>10e-15) ||(splitPoint[k].flux_res[4]<0 && 0-splitPoint[k].flux_res[4]>10e-15)) ){
		// 	printf("Interior Here %d\n",k);
		// }
		// if(splitPoint[k].flux_res[r]!=splitPoint[k].flux_res[r]){
		// 	printf("Interior Here %d\n",splitPoint[k].globalIndex);
		// }
	}
	//
	res_sqr = (U[0] - temp) * (U[0] - temp);
	sum_res_sqr[k] = res_sqr;
	//
	conserved_to_primitive_multi_nccl(splitPoint, k, U);
}

__global__ void state_update_outer_multi_nccl(int myRank,splitPoints *splitPoint, int outerPointsLocal, int *outerPointsLocalIndex, double u1_inf, double u2_inf, double u3_inf, double rho_inf, double pi, double pr_inf)
{
	int k;
	double U[5];
	//
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < 0 || i >= outerPointsLocal)
	{
		return;
	}
	k = outerPointsLocalIndex[i];
	conserved_vector_Ubar_multi_cuda(splitPoint, k, U, u1_inf, u2_inf, u3_inf, rho_inf, pi, pr_inf);
	// temp = U[0];
	for (int r = 0; r < 5; r++)
	{
		U[r] = U[r] - splitPoint[k].flux_res[r];
		// if(splitPoint[k].prim[0]<0){
		// 	printf("Outer splitPoint[%d].prim[%d] %.15f",splitPoint[k].globalIndex,r,splitPoint[k].prim[0]);
		// }
		// if(k==249755 && myRank==1){
		// 	printf("splitPoint[k].flux_res[r] %.15f\n",splitPoint[k].flux_res[r]);
		// }
		// if(splitPoint[k].flux_res[r]!=splitPoint[k].flux_res[r]){
		// 	printf("Outer Here %d\n",k);
		// }
	}
	//

	conserved_to_primitive_multi_nccl(splitPoint, k, U);
}

__global__ void state_update_symmetric_multi_nccl(splitPoints *splitPoint, double power, double VL_CONST, double pi, int symmetryPointsLocal, int *symmetryPointsLocalIndex,int* globalToLocalIndex,int **globalToGhostIndex,transferPoints **receiveBuffer,int *partVector){
	int k;
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	double delx, delz;
	if (i < 0 || i >= symmetryPointsLocal)
	{
		return;
	}
	k=symmetryPointsLocalIndex[i];

	for (int j=0;j<splitPoint[k].numberOfLocalNbhs;j++){
		int nbh = splitPoint[k].localNbhs[j];
		nbh = globalToLocalIndex[nbh];
		delx = splitPoint[nbh].x-splitPoint[k].x;
		delz = splitPoint[nbh].z-splitPoint[k].z;
		if(abs(delx) <10e-9 && abs(delz) <10e-9){
			for(int r=0;r<5;r++){
				splitPoint[k].prim[r]=splitPoint[nbh].prim[r];
			}
		}
	}
	for (int j = 0; j < splitPoint[k].numberOfGhostNbhs; j++)
    {
        int nbh = splitPoint[k].ghostNbhs[j];
        int device=partVector[nbh];
        int ghostIndex=globalToGhostIndex[device][nbh];
		delx = receiveBuffer[device][ghostIndex].x-splitPoint[k].x;
		delz = receiveBuffer[device][ghostIndex].z-splitPoint[k].z;
		if(abs(delx) <10e-9 && abs(delz) <10e-9){
			for(int r=0;r<5;r++){
				splitPoint[k].prim[r]=receiveBuffer[device][ghostIndex].prim[r];
			}
		}
	}
}