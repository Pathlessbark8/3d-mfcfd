
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
    double tan1[3][max_points], tan2[3][max_points], nor[3][max_points];
    int status[max_points];
    int nbhs[max_points];
    int conn[27][max_points];
    double min_dist[max_points];
    //
    int xpos_nbhs[max_points], xneg_nbhs[max_points];
    int ypos_nbhs[max_points], yneg_nbhs[max_points];
    int zpos_nbhs[max_points], zneg_nbhs[max_points];
    int xpos_conn[22][max_points], xneg_conn[22][max_points];
    int ypos_conn[22][max_points], yneg_conn[22][max_points];
    int zpos_conn[22][max_points], zneg_conn[22][max_points];
    //
    //		flow field based attributes ..
    //
    double U[5][max_points];
    double delUp[5][max_points];
    double delUn[5][max_points];
    double prim[5][max_points];
    double flux_res[5][max_points];
    double q[5][max_points];
    double qm[2][5][max_points];
    double dq[3][5][max_points];
    double temp[3][5][max_points];
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
int interior_points, wall_points, outer_points;
int supersonic_inlet_points, supersonic_outlet_points;
//
int interior_points_index[max_points];
int wall_points_index[max_points];
int outer_points_index[max_points];
int supersonic_outlet_points_index[max_points];
int supersonic_inlet_points_index[max_points];
//
double res_old, res_new, max_res, residue;
int max_res_point;

int threads_per_block = 128;
//
//

//
