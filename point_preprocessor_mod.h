/*  
	MFCFD is a 3D Computational Fluid Dynamics Solver based off q-LSKUM
    Copyright (C) 2022 Dhruv Saxena
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
    //		
    std::fstream fin;
    fin.open("/home/anil/new_3d_code/3d-mfcfd/inputFiles/3d-grid-580485.dat", std::ios::in);
    //
    for (k = 0; k < max_points; k++)
    {
        fin >> counter  >> point.x[k] >> point.y[k] >> point.z[k]>> point.status[k]>>point.min_dist[k];
        fin >> point.tan1[0][k] >> point.tan1[1][k] >> point.tan1[2][k];
        fin >> point.tan2[0][k] >> point.tan2[1][k] >> point.tan2[2][k];
        fin >> point.nor[0][k] >> point.nor[1][k] >> point.nor[2][k];
        fin >> point.nbhs[k];
        for (r = 0; r < point.nbhs[k]; r++)
        {
            fin >> point.conn[r][k];
            point.conn[r][k]-=1;
        }
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
