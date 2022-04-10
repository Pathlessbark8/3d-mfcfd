
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
    int counter[max_points];
    //
    //
};
//
points point;


struct splitPoints{
    double x,y,z;
    double tan1[3],tan2[3],nor[3];
    int status;
    int nbhs;
    int conn[27];
    double min_dist;
    int xpos_nbhs,xneg_nbhs;
    int ypos_nbhs,yneg_nbhs;
    int zpos_nbhs,zneg_nbhs;
    int xpos_conn[22],xneg_conn[22];
    int ypos_conn[22],yneg_conn[22];
    int zpos_conn[22],zneg_conn[22];
    double U[5];
    double delUp[5];
    double delUn[5];
    double prim[5];
    double flux_res[5];
    double q[5];
    double qm[2][5];
    double dq[3][5];
    double temp[3][5];
    double delt;
    int alias;
    int point_with_alias;
    int counter;
};

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
double sum_res_sqr=0;
//
splitPoints *splitPoint[max_devices];

void assign(splitPoints &splitPoint,int i){
    splitPoint.x=point.x[i];
    splitPoint.y=point.y[i];
    splitPoint.z=point.z[i];
    splitPoint.status=point.status[i];
    splitPoint.nbhs=point.nbhs[i];
    splitPoint.min_dist=point.min_dist[i];  
    splitPoint.alias=point.alias[i];
    splitPoint.point_with_alias=point.point_with_alias[i];
    splitPoint.counter=point.counter[i];
    for(int j=0;j<27;j++){
        splitPoint.conn[j]=point.conn[j][i];
    }
    for(int j=0;j<3;j++){
        splitPoint.tan1[j]=point.tan1[j][i];
        splitPoint.tan2[j]=point.tan2[j][i];
        splitPoint.nor[j]=point.nor[j][i];
    }
    for(int j=0;j<5;j++){
        splitPoint.U[j]=point.U[j][i];
        splitPoint.delUp[j]=point.delUp[j][i];
        splitPoint.delUn[j]=point.delUn[j][i];
        splitPoint.prim[j]=point.prim[j][i];
        splitPoint.flux_res[j]=point.flux_res[j][i];
        splitPoint.q[j]=point.q[j][i];
        splitPoint.qm[0][j]=point.qm[0][j][i];
        splitPoint.qm[1][j]=point.qm[1][j][i];
        splitPoint.dq[0][j]=point.dq[0][j][i];
        splitPoint.dq[1][j]=point.dq[1][j][i];
        splitPoint.dq[2][j]=point.dq[2][j][i];
        splitPoint.temp[0][j]=point.temp[0][j][i];
        splitPoint.temp[1][j]=point.temp[1][j][i];
        splitPoint.temp[2][j]=point.temp[2][j][i];
    }
    splitPoint.delt=point.delt[i];
}
//

//
