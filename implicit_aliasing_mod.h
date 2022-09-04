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
//	First written on 23.06.2021.
//
#include "data_structure_mod.h"
//
//
//
//	This subroutine sets the alias indexing for all points ..
//
//
void aliasing()
//
//
{
    //
    int i, j, p, k, count;
    //
    //
    for (i = 0; i < max_points; i++)
    {
        point.alias[i] = 0;
        point.point_with_alias[i] = 0;
    }
    //
    //
    point.alias[0] = 1;
    point.point_with_alias[0] = 1;
    count = 0;
    //
    //
    for (i = 0; i < max_points; i++)
    {
        p = point.point_with_alias[i];
        for (j = 0; j < point.nbhs[p]; j++)
        {
            k = point.conn[j][p];
            if (point.alias[k] == 0)
            {
                count = count + 1;
                point.alias[k] = count;
                point.point_with_alias[count] = k;
            }
        }
    }
    //
}
//
