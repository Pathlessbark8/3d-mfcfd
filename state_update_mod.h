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
