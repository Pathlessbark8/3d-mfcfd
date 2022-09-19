
//
#pragma once
#include "parameter_mod.h"
//
//
struct points
{
    //
    //		geometry based attributes ..
    //
    double x[max_points], y[max_points], z[max_points];
    //		double  nor_neigh
    double tan1[max_points][3], tan2[max_points][3], nor[max_points][3];
    int status[max_points];
    int nbhs[max_points];
    int conn[max_points][27];
    double min_dist[max_points];
    //
    int xpos_nbhs[max_points], xneg_nbhs[max_points];
    int ypos_nbhs[max_points], yneg_nbhs[max_points];
    int zpos_nbhs[max_points], zneg_nbhs[max_points];
    int xpos_conn[max_points][27], xneg_conn[max_points][27];
    int ypos_conn[max_points][27], yneg_conn[max_points][27];
    int zpos_conn[max_points][27], zneg_conn[max_points][27];
    //
    //		flow field based attributes ..
    //
    double U[max_points][5];
    double delUp[max_points][5];
    double delUn[max_points][5];
    double prim[max_points][5];
    double flux_res[max_points][5];
    double q[max_points][5];
    double qm[max_points][2][5];
    double dq[max_points][3][5];
    double temp[max_points][3][5];
    double delt[max_points];
    //
    //		for the implicit solver ..
    //
    int alias[max_points], point_with_alias[max_points];
    //
    //
};
//
points point;
//
//
int interior_points, wall_points, outer_points,symmetry_points;
int supersonic_inlet_points, supersonic_outlet_points;
//
int interior_points_index[max_points];
int wall_points_index[max_points];
int outer_points_index[max_points];
int supersonic_outlet_points_index[max_points];
int supersonic_inlet_points_index[max_points];
int symmetry_points_index[max_points];
//
double res_old, res_new, max_res, residue;
int max_res_point;

int threads_per_block = 128;
double sum_res_sqr=0;
//
//

//
