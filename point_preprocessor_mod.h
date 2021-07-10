//
//	First written on 18.06.2021 ..
//
//
#pragma once
//
//
#include "data_structure_mod.h"
#include "fstream"
//
//
//
void read_input_point_data()
//
//
{
    //
    int k, r, counter;
    int interior_count, wall_count, outer_count;
    int supersonic_inlet_count, supersonic_outlet_count;
    //
    //
    //
    //		OPEN(UNIT=101,FILE="hemisphere-testcase/3d_input_data",FORM="FORMATTED",STATUS="OLD",ACTION="READ")
    //OPEN(UNIT=101,FILE="3d_input_data",FORM="FORMATTED",STATUS="OLD",ACTION="READ")
    std::fstream fin;
    fin.open("/home/dhruv/3d_new_data", std::ios::in);
    //
    for (k = 0; k < max_points; k++)
    {
        fin >> counter >> point.status[k] >> point.x[k] >> point.y[k] >> point.z[k];
        fin >> point.tan1[0][k] >> point.tan1[1][k] >> point.tan1[2][k];
        fin >> point.tan2[0][k] >> point.tan2[1][k] >> point.tan2[2][k];
        fin >> point.nor[0][k] >> point.nor[1][k] >> point.nor[2][k];
        fin >> point.nbhs[k];
        for (r = 0; r < point.nbhs[k]; r++)
        {
            fin >> point.conn[r][k];
            point.conn[r][k] -= 1;
        }
        fin >> point.min_dist[k];
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
    //
    for (k = 0; k < max_points; k++)
    {
        if (point.status[k] == 0)
            interior_points = interior_points + 1;
        else if (point.status[k] == 1)
            wall_points = wall_points + 1;
        else if (point.status[k] == 2)
            outer_points = outer_points + 1;
        else if (point.status[k] == 6)
            supersonic_outlet_points = supersonic_outlet_points + 1;
        else if (point.status[k] == 5)
            supersonic_inlet_points = supersonic_inlet_points + 1;
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
//
//
