#pragma once
//
//
//
//
#include "data_structure_mod.h"
#include "primitive_to_conserved_vector_mod.h"
#include "conserved_to_primitive_vector_mod.h"
#include "conserved_vector_Ubar_mod.h"
//
//
void state_update()
//
//
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
	for(int i=0;i<symmetry_points;i++){
		k=symmetry_points_index[i];
		double dely,delx;
		for (int j=0;j<point.nbhs[k];j++){
			int nbh = point.conn[j][k];
			delx = point.x[nbh]-point.x[k];
			dely = point.y[nbh]-point.y[k];
			if(delx <10e-9 && dely <10e-9){
				for(int r=0;r<5;r++){
					point.prim[r][k]=point.prim[r][k];
				}
			}
		}
	}
	res_new = sqrt(sum_res_sqr) / max_points;
	//
	//
}
//
//
//
