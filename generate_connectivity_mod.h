//
//	First written on 23.06.2021
//
//
#pragma once

#include "data_structure_mod.h"
//
//

//
void get_interior_point_neighbours(int i)
//
{
    //
    int k, r, count;
    double xi, yi, zi, xk, yk, zk;
    double delx, dely, delz;
    double dels, delt, deln;
    //
    xi = point.x[i];
    yi = point.y[i];
    zi = point.z[i];
    //
    point.xpos_nbhs[i] = 0;
    point.xneg_nbhs[i] = 0;
    point.ypos_nbhs[i] = 0;
    point.yneg_nbhs[i] = 0;
    point.zpos_nbhs[i] = 0;
    point.zneg_nbhs[i] = 0;
    //
    for (r = 0; r < point.nbhs[i]; r++)
    {
        k = point.conn[r][i];
        xk = point.x[k];
        yk = point.y[k];
        zk = point.z[k];
        //
        delx = xk - xi;
        dely = yk - yi;
        delz = zk - zi;
        //
        dels = delx * point.tan1[0][i] + dely * point.tan1[1][i] + delz * point.tan1[2][i];
        delt = delx * point.tan2[0][i] + dely * point.tan2[1][i] + delz * point.tan2[2][i];
        deln = delx * point.nor[0][i] + dely * point.nor[1][i] + delz * point.nor[2][i];
        //
        if (dels <= 0.00)
        {
            point.xpos_nbhs[i] = point.xpos_nbhs[i] + 1;
            count = point.xpos_nbhs[i];
            point.xpos_conn[count-1][i] = k;
        }
        //
        if (dels > 0.00)
        {
            point.xneg_nbhs[i] = point.xneg_nbhs[i] + 1;
            count = point.xneg_nbhs[i];
            point.xneg_conn[count-1][i] = k;
        }
        //
        if (delt <= 0.00)
        {
            point.ypos_nbhs[i] = point.ypos_nbhs[i] + 1;
            count = point.ypos_nbhs[i];
            point.ypos_conn[count-1][i] = k;
        }
        //
        if (delt > 0.00)
        {
            point.yneg_nbhs[i] = point.yneg_nbhs[i] + 1;
            count = point.yneg_nbhs[i];
            point.yneg_conn[count-1][i] = k;
        }
        //
        if (deln <= 0.00)
        {
            point.zpos_nbhs[i] = point.zpos_nbhs[i] + 1;
            count = point.zpos_nbhs[i];
            point.zpos_conn[count-1][i] = k;
        }
        //
        if (deln > 0.00)
        {
            point.zneg_nbhs[i] = point.zneg_nbhs[i] + 1;
            count = point.zneg_nbhs[i];
            point.zneg_conn[count-1][i] = k;
        }
        //
    }
    //
    //
}
//
//
//
void get_wall_point_neighbours(int i)
//
{
    //
    int k, r, count;
    double xi, yi, zi, xk, yk, zk;
    double delx, dely, delz;
    double dels, delt, deln;
    //
    xi = point.x[i];
    yi = point.y[i];
    zi = point.z[i];
    //
    point.xpos_nbhs[i] = 0;
    point.xneg_nbhs[i] = 0;
    point.ypos_nbhs[i] = 0;
    point.yneg_nbhs[i] = 0;
    point.zneg_nbhs[i] = 0;
    //
    for (r = 0; r < point.nbhs[i]; r++)
    {
        k = point.conn[r][i];
        xk = point.x[k];
        yk = point.y[k];
        zk = point.z[k];
        //
        delx = xk - xi;
        dely = yk - yi;
        delz = zk - zi;
        //
        dels = delx * point.tan1[0][i] + dely * point.tan1[1][i] + delz * point.tan1[2][i];
        delt = delx * point.tan2[0][i] + dely * point.tan2[1][i] + delz * point.tan2[2][i];
        deln = delx * point.nor[0][i] + dely * point.nor[1][i] + delz * point.nor[2][i];
        //
        //
        if (dels <= 0.00)
        {
            point.xpos_nbhs[i] = point.xpos_nbhs[i] + 1;
            count = point.xpos_nbhs[i];
            point.xpos_conn[count-1][i] = k;
        }
        //
        if (dels >= 0.00)
        {
            point.xneg_nbhs[i] = point.xneg_nbhs[i] + 1;
            count = point.xneg_nbhs[i];
            point.xneg_conn[count-1][i] = k;
        }
        //
        if (delt <= 0.00)
        {
            point.ypos_nbhs[i] = point.ypos_nbhs[i] + 1;
            count = point.ypos_nbhs[i];
            point.ypos_conn[count-1][i] = k;
        }
        //
        if (delt >= 0.00)
        {
            point.yneg_nbhs[i] = point.yneg_nbhs[i] + 1;
            count = point.yneg_nbhs[i];
            point.yneg_conn[count-1][i] = k;
        }
        //
        point.zneg_nbhs[i] = point.zneg_nbhs[i] + 1;
        count = point.zneg_nbhs[i];
        point.zneg_conn[count-1][i] = k;
    }
    //
    //
}
//
//
//
void get_outer_point_neighbours(int i)
//
//
{
    //
    int k, r, count;
    double xi, yi, zi, xk, yk, zk;
    double delx, dely, delz;
    double dels, delt, deln;
    //
    xi = point.x[i];
    yi = point.y[i];
    zi = point.z[i];
    //
    point.xpos_nbhs[i] = 0;
    point.xneg_nbhs[i] = 0;
    point.ypos_nbhs[i] = 0;
    point.yneg_nbhs[i] = 0;
    point.zneg_nbhs[i] = 0;
    //
    for (r = 0; r < point.nbhs[i]; r++)
    {
        k = point.conn[r][i];
        xk = point.x[k];
        yk = point.y[k];
        zk = point.z[k];
        //
        delx = xk - xi;
        dely = yk - yi;
        delz = zk - zi;
        //
        dels = delx * point.tan1[0][i] + dely * point.tan1[1][i] + delz * point.tan1[2][i];
        delt = delx * point.tan2[0][i] + dely * point.tan2[1][i] + delz * point.tan2[2][i];
        deln = delx * point.nor[0][i] + dely * point.nor[1][i] + delz * point.nor[2][i];
        //
        if (dels <= 0.00)
        {
            point.xpos_nbhs[i] = point.xpos_nbhs[i] + 1;
            count = point.xpos_nbhs[i];
            point.xpos_conn[count-1][i] = k;
        }
        //
        if (dels >= 0.00)
        {
            point.xneg_nbhs[i] = point.xneg_nbhs[i] + 1;
            count = point.xneg_nbhs[i];
            point.xneg_conn[count-1][i] = k;
        }
        //
        if (delt <= 0.00)
        {
            point.ypos_nbhs[i] = point.ypos_nbhs[i] + 1;
            count = point.ypos_nbhs[i];
            point.ypos_conn[count-1][i] = k;
        }
        //
        if (delt >= 0.00)
        {
            point.yneg_nbhs[i] = point.yneg_nbhs[i] + 1;
            count = point.yneg_nbhs[i];
            point.yneg_conn[count-1][i] = k;
        }
        //
        point.zpos_nbhs[i] = point.zpos_nbhs[i] + 1;
        count = point.zpos_nbhs[i];
        point.zpos_conn[count-1][i] = k;
    }
    //
    //
}
//
//

//
void generate_split_stencils()
//
{
    //
    int i, k;

    for (k = 0; k < interior_points; k++)
    {
        i = interior_points_index[k];
        get_interior_point_neighbours(i);
    }
    //
    for (k = 0; k < wall_points; k++)
    {
        i = wall_points_index[k];
        get_wall_point_neighbours(i);
    }
    //
    for (k = 0; k < outer_points; k++)
    {
        i = outer_points_index[k];
        get_outer_point_neighbours(i);
    }
    //
}
//
//
