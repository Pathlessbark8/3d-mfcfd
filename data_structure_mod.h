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
#include "parameter_mod.h"

//
//
struct points
{
    //
    //		geometry based attributes ..
    //
    double x[max_points], y[max_points], z[max_points];
    //	
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
    double flux_res[5]={0};
    double q[5];
    double qm[2][5];
    double dq[3][5];
    double temp[3][5];
    double delt;
    int alias;
    int point_with_alias;
    int numberOfGhostNbhs;
    int numberOfLocalNbhs;
    int ghostNbhs[27];
    int localNbhs[27];
    int numberOfLocalxposNbhs;
    int numberOfLocalxnegNbhs;
    int numberOfLocalyposNbhs;
    int numberOfLocalynegNbhs;
    int numberOfLocalzposNbhs;
    int numberOfLocalznegNbhs;
    int numberOfGhostxposNbhs;
    int numberOfGhostxnegNbhs;
    int numberOfGhostyposNbhs;
    int numberOfGhostynegNbhs;
    int numberOfGhostzposNbhs;
    int numberOfGhostznegNbhs;
    int localxpos_conn[22],localxneg_conn[22];
    int localypos_conn[22],localyneg_conn[22];
    int localzpos_conn[22],localzneg_conn[22];
    int ghostxpos_conn[22],ghostxneg_conn[22];
    int ghostypos_conn[22],ghostyneg_conn[22];
    int ghostzpos_conn[22],ghostzneg_conn[22];
    bool isGhost;
    int globalIndex;
    int ghostIndex[max_partitions];
    int numberOfPartitionsToSendTo;
    int partitions[max_partitions];
};


struct transferPoints{
    int globalIndex;
    double x,y,z;
    double q[5];
    double dq[3][5];
    double prim[5];
    double delt;
    double qm[2][5];
    double min_dist;
};

transferPoints **sendBuffer;
transferPoints **receiveBuffer;

int globalToLocalIndex[max_points]={0};
int **globalToGhostIndex;
int *localToGlobalIndex;
int **ghostToGlobalIndex;
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

// For the parallel version
int interiorPointsLocal=0, wallPointsLocal=0, outerPointsLocal=0, symmetryPointsLocal=0;
int supersonicInletPointsLocal=0, supersonicOutletPointsLocal=0;
int *interiorPointsLocalIndex;
int *wallPointsLocalIndex;
int *outerPointsLocalIndex;
int *supersonicOutletPointsLocalIndex;
int *supersonicInletPointsLocalIndex;
int *symmetryPointsLocalIndex;
//
//
double res_old, res_new, max_res, residue;
int max_res_point;

int threads_per_block = 128;
double sum_res_sqr=0;
//
splitPoints *splitPoint;

void findNatureOfLocalPoints(splitPoints &splitPoint){
    if (splitPoint.status == 0)
            interiorPointsLocal = interiorPointsLocal + 1;
        else if (splitPoint.status == 1)
            wallPointsLocal = wallPointsLocal + 1;
        else if (splitPoint.status == 2)
            outerPointsLocal = outerPointsLocal + 1;
        else if (splitPoint.status == 3)
            symmetryPointsLocal = symmetryPointsLocal + 1;
        else if (splitPoint.status== 6)
        {
            supersonicOutletPointsLocal = supersonicOutletPointsLocal + 1;
        }
        else if (splitPoint.status== 5)
        {
            supersonicInletPointsLocal = supersonicInletPointsLocal + 1;
        }
}

void allocateSizeForNatureOfLocalPoints(){
    interiorPointsLocalIndex=new int[interiorPointsLocal];
    wallPointsLocalIndex=new int[wallPointsLocal];
    outerPointsLocalIndex=new int[outerPointsLocal];
    supersonicOutletPointsLocalIndex=new int[supersonicOutletPointsLocal];
    supersonicInletPointsLocalIndex=new int[supersonicInletPointsLocal];
    symmetryPointsLocalIndex=new int[symmetryPointsLocal];
    interiorPointsLocal=0;
    wallPointsLocal=0;
    outerPointsLocal=0;
    supersonicOutletPointsLocal=0;
    supersonicInletPointsLocal=0;
    symmetryPointsLocal=0;
}

void assignNatureOfLocalPoints(splitPoints &splitPoint,int k){
    if (splitPoint.status == 0)
        {
            interiorPointsLocalIndex[interiorPointsLocal] = k;
            // cout<<"Verify "<<k<<endl;
            interiorPointsLocal = interiorPointsLocal + 1;
        }
        else if (splitPoint.status == 1)
        { 
            wallPointsLocalIndex[wallPointsLocal] = k;
            wallPointsLocal = wallPointsLocal + 1;
        }
        else if (splitPoint.status == 2)
        {
            outerPointsLocalIndex[outerPointsLocal] = k;
            outerPointsLocal = outerPointsLocal + 1;
        }
        else if (splitPoint.status == 3)
        {
            symmetryPointsLocalIndex[symmetryPointsLocal] = k;
            symmetryPointsLocal = symmetryPointsLocal + 1;
        }
        else if (splitPoint.status == 6)
        {
            supersonicOutletPointsLocalIndex[supersonicOutletPointsLocal] = k;
            supersonicOutletPointsLocal = supersonicOutletPointsLocal + 1;
        }
        else if (splitPoint.status == 5)
        {
            supersonicInletPointsLocalIndex[supersonicInletPointsLocal] = k;
            supersonicInletPointsLocal = supersonicInletPointsLocal + 1;
        }
}

void assign(splitPoints &splitPoint,int i,int myRank){
    splitPoint.globalIndex=i;
    splitPoint.x=point.x[i];
    splitPoint.y=point.y[i];
    splitPoint.z=point.z[i];
    splitPoint.status=point.status[i];
    splitPoint.nbhs=point.nbhs[i];
    splitPoint.min_dist=point.min_dist[i];  
    splitPoint.alias=point.alias[i];
    splitPoint.point_with_alias=point.point_with_alias[i];
    //
    splitPoint.numberOfGhostNbhs=0;
    for(int j=0;j<point.nbhs[i];j++){
        splitPoint.conn[j]=point.conn[j][i];
        if(partVector[splitPoint.conn[j]]!=partVector[i]){
            splitPoint.numberOfGhostNbhs++;
        }
    }
    splitPoint.numberOfLocalNbhs=splitPoint.nbhs-splitPoint.numberOfGhostNbhs;
    if(splitPoint.numberOfGhostNbhs!=0){
        splitPoint.numberOfGhostNbhs=0;
    }

    splitPoint.numberOfLocalNbhs=0;
    for(int j=0;j<point.nbhs[i];j++){
        if(partVector[splitPoint.conn[j]]!=partVector[i]){
            splitPoint.ghostNbhs[splitPoint.numberOfGhostNbhs++]=splitPoint.conn[j];
        }
        else{
            splitPoint.localNbhs[splitPoint.numberOfLocalNbhs++]=splitPoint.conn[j];
        }
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
    splitPoint.isGhost=false;
    splitPoint.numberOfPartitionsToSendTo=0;
}

//

//