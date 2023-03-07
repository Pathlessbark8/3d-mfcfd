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
//	First written on 18.06.2021 ..
//
//
#pragma once
//
//
#include "data_structure_mod.h"
#include "fstream"
#include <string.h>
//
//
//
void read_input_point_data()
//
//
{
    //
    int k, r, counter;
    int interior_count, wall_count, outer_count, symmetry_count;
    int supersonic_inlet_count, supersonic_outlet_count;
    //
    //
    //
    //		
    std::fstream fin;
    fin.open("/home/nsm/3d-mfcfd/inputFiles/"+to_string(max_points)+"/partGrid-"+to_string(max_points)+".dat", std::ios::in);
    cout<<"/home/nsm/3d-mfcfd/inputFiles/"+to_string(max_points)+"/partGrid-"+to_string(max_points)+".dat"<<endl;
    //
    for (k = 0; k < max_points; k++)
    {
        fin >> counter  >> point.x[k] >> point.y[k] >> point.z[k]>> point.status[k]>>point.min_dist[k];
        fin >> point.tan1[0][k] >> point.tan1[1][k] >> point.tan1[2][k];
        fin >> point.tan2[0][k] >> point.tan2[1][k] >> point.tan2[2][k];
        fin >> point.nor[0][k] >> point.nor[1][k] >> point.nor[2][k];
        if(point.status[k]==0){
            point.nor[2][k]=1;
        }
        fin >> point.nbhs[k];
        for (r = 0; r < point.nbhs[k]; r++)
        {
            fin >> point.conn[r][k];
            point.conn[r][k]-=1;
            // if(k==9){
            //     cout<<point.conn[r][k]<<endl;
            // }
        }
        // if(k==1){
        //     cout<<counter  <<" "<< point.x[k] <<" "<< point.y[k] <<" "<< point.z[k]<<" "<< point.status[k]<<" "<<point.min_dist[k];
        //     cout<<" "<< point.tan1[0][k] <<" "<< point.tan1[1][k] <<" "<< point.tan1[2][k];
        // cout<<" "<< point.tan2[0][k] <<" "<< point.tan2[1][k] <<" "<< point.tan2[2][k];
        // cout<<" "<< point.nor[0][k]<<" "<< point.nor[1][k] <<" "<< point.nor[2][k];
        // cout<<" "<< point.nbhs[k];
        // }
        // fin >> point.min_dist[k];
    }
    fin.close();
    //
    //		Finding the number of interior, wall, outer and other boundary points ..
    //
    //
    interior_points = 0;
    wall_points = 0;
    outer_points = 0;
    supersonic_outlet_points = 0;
    supersonic_inlet_points = 0;
    symmetry_points = 0;

    //
    for (k = 0; k < max_points; k++)
    {
        if (point.status[k] == 0)
            interior_points = interior_points + 1;
        else if (point.status[k] == 1)
            wall_points = wall_points + 1;
        else if (point.status[k] == 2)
            outer_points = outer_points + 1;
        else if(point.status[k] == 3)
            symmetry_points=symmetry_points+1;
        else if (point.status[k] == 6)
        {
            supersonic_outlet_points = supersonic_outlet_points + 1;
            // cout<<k<<endl;
        }
        else if (point.status[k] == 5)
        {
            supersonic_inlet_points = supersonic_inlet_points + 1;
            // cout<<k<<"Here"<<endl;
        }

    }
    //
    //		Allocating the size of the respective points ..
    //
    //allocate(interior_points_index(interior_points))
    // allocate(wall_points_index(wall_points))
    // allocate(outer_points_index(outer_points))
    // allocate(supersonic_outlet_points_index(supersonic_outlet_points))
    // allocate(supersonic_inlet_points_index(supersonic_inlet_points))
    //
    //
    //		Finding the indices of the interior, wall, outer and other boundary points ..
    //
    interior_count = -1;
    wall_count = -1;
    outer_count = -1;
    supersonic_inlet_count = -1;
    supersonic_outlet_count = -1;
    symmetry_count=-1;
    //
    for (k = 0; k < max_points; k++)
    {
        if (point.status[k] == 0)
        {
            interior_count = interior_count + 1;
            interior_points_index[interior_count] = k;
        }
        else if (point.status[k] == 1)
        {
            wall_count = wall_count + 1;
            wall_points_index[wall_count] = k;
        }
        else if (point.status[k] == 2)
        {
            outer_count = outer_count + 1;
            outer_points_index[outer_count] = k;
        }
        else if (point.status[k] == 3){
            symmetry_count=symmetry_count+1;
            symmetry_points_index[symmetry_count] = k;
        }
        else if (point.status[k] == 6)
        {
            supersonic_outlet_count = supersonic_outlet_count + 1;
            supersonic_outlet_points_index[supersonic_outlet_count] = k;
        }
        else if (point.status[k] == 5)
        {
            supersonic_inlet_count = supersonic_inlet_count + 1;
            supersonic_inlet_points_index[supersonic_inlet_count] = k;
        }
    }
    //
    //
}

void read_input_point_data_multi()
//
//
{
    //		
    std::fstream fin;
    fin.open("/home/nsm/3d-mfcfd/inputFiles/"+to_string(max_points)+"/partGrid-"+to_string(max_points)+".dat", std::ios::in);
    cout<<"/home/nsm/3d-mfcfd/inputFiles/"+to_string(max_points)+"/partGrid-"+to_string(max_points)+".dat"<<endl;
    //
    int r,k,i=0;
    double x,y,z,min_dist,tan1[3],tan2[3],nor[3];
    int counter,status,nbhs,conn[27];
    for (k = 0; k < max_points && i<local_points; k++)
    {
        fin>> counter >> x >> y >> z >> status >> min_dist;
        fin >> tan1[0] >> tan1[1] >> tan1[2];
        fin >> tan2[0] >> tan2[1] >> tan2[2];
        fin >> nor[0] >> nor[1] >> nor[2];
        fin >> nbhs;
        for (r = 0; r < nbhs; r++)
        {
            fin >> conn[r];
            conn[r]-=1;
        }
        if(localToGlobalIndex[i]==k){
            splitPoint[i].globalIndex=k;
            splitPoint[i].x = x;
            splitPoint[i].y = y;
            splitPoint[i].z = z;
            splitPoint[i].status = status;
            splitPoint[i].min_dist = min_dist;
            for(int r=0;r<3;r++){
                splitPoint[i].tan1[r]=tan1[r];
                splitPoint[i].tan2[r]=tan2[r];
                splitPoint[i].nor[r]=nor[r];
            }
            splitPoint[i].nbhs=nbhs;
            for(int r=0;r<nbhs;r++){
                splitPoint[i].conn[r]=conn[r];
            }


            splitPoint[i].numberOfGhostNbhs=0;
            for(int j=0;j<splitPoint[i].nbhs;j++){
                if(partVector[splitPoint[i].conn[j]]!=partVector[k]){
                    splitPoint[i].numberOfGhostNbhs++;
                }
            }
            splitPoint[i].numberOfLocalNbhs=splitPoint[i].nbhs-splitPoint[i].numberOfGhostNbhs;
            if(splitPoint[i].numberOfGhostNbhs!=0){
                splitPoint[i].numberOfGhostNbhs=0;
            }

            splitPoint[i].numberOfLocalNbhs=0;
            for(int j=0;j<splitPoint[i].nbhs;j++){
                if(partVector[splitPoint[i].conn[j]]!=partVector[k]){
                    splitPoint[i].ghostNbhs[splitPoint[i].numberOfGhostNbhs++]=splitPoint[i].conn[j];
                }
                else{
                    splitPoint[i].localNbhs[splitPoint[i].numberOfLocalNbhs++]=splitPoint[i].conn[j];
                }
            }
            // splitPoint[i].delt=point.delt[i];
            splitPoint[i].isGhost=false;
            splitPoint[i].numberOfPartitionsToSendTo=0;
            i++;
        }
    }
    fin.close();
    //
    //		Finding the number of interior, wall, outer and other boundary points ..
    //
    //
    interiorPointsLocal=0;
    wallPointsLocal=0;
    outerPointsLocal=0;
    symmetryPointsLocal=0;
    supersonicInletPointsLocal=0;
    supersonicOutletPointsLocal=0;

    //
    for (k = 0; k < local_points; k++)
    {
        if (splitPoint[k].status == 0)
            interiorPointsLocal = interiorPointsLocal + 1;
        else if (splitPoint[k].status == 1)
            wallPointsLocal = wallPointsLocal + 1;
        else if (splitPoint[k].status == 2)
            outerPointsLocal = outerPointsLocal + 1;
        else if (splitPoint[k].status == 3)
            symmetryPointsLocal = symmetryPointsLocal + 1;
        else if (splitPoint[k].status== 6)
        {
            supersonicOutletPointsLocal = supersonicOutletPointsLocal + 1;
        }
        else if (splitPoint[k].status== 5)
        {
            supersonicInletPointsLocal = supersonicInletPointsLocal + 1;
        }

    }
    //
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
    //
    for (k = 0; k < local_points; k++)
    {
        if (splitPoint[k].status == 0)
        {
            interiorPointsLocalIndex[interiorPointsLocal] = k;
            // cout<<"Verify "<<k<<endl;
            interiorPointsLocal = interiorPointsLocal + 1;
        }
        else if (splitPoint[k].status == 1)
        { 
            wallPointsLocalIndex[wallPointsLocal] = k;
            wallPointsLocal = wallPointsLocal + 1;
        }
        else if (splitPoint[k].status == 2)
        {
            outerPointsLocalIndex[outerPointsLocal] = k;
            outerPointsLocal = outerPointsLocal + 1;
        }
        else if (splitPoint[k].status == 3)
        {
            symmetryPointsLocalIndex[symmetryPointsLocal] = k;
            symmetryPointsLocal = symmetryPointsLocal + 1;
        }
        else if (splitPoint[k].status == 6)
        {
            supersonicOutletPointsLocalIndex[supersonicOutletPointsLocal] = k;
            supersonicOutletPointsLocal = supersonicOutletPointsLocal + 1;
        }
        else if (splitPoint[k].status == 5)
        {
            supersonicInletPointsLocalIndex[supersonicInletPointsLocal] = k;
            supersonicInletPointsLocal = supersonicInletPointsLocal + 1;
        }
    }
    
    
}
//
//
