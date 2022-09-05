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
            k = point.conn[p][j];
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
