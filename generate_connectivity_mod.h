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
        k = point.conn[i][r];
        // if(i==21){
        //         printf("%d %d \n",i,k);
        // }
        xk = point.x[k];
        yk = point.y[k];
        zk = point.z[k];
        //
        delx = xk - xi;
        dely = yk - yi;
        delz = zk - zi;
        //
        dels = delx * point.tan1[i][0] + dely * point.tan1[i][1] + delz * point.tan1[i][2];
        delt = delx * point.tan2[i][0] + dely * point.tan2[i][1] + delz * point.tan2[i][2];
        deln = delx * point.nor[i][0] + dely * point.nor[i][1] + delz * point.nor[i][2];
        //
        if (dels <= 0.00)
        {
            point.xpos_nbhs[i] = point.xpos_nbhs[i] + 1;
            count = point.xpos_nbhs[i];
            point.xpos_conn[i][count-1] = k;
        }
        //
        if (dels > 0.00)
        {
            point.xneg_nbhs[i] = point.xneg_nbhs[i] + 1;
            count = point.xneg_nbhs[i];
            point.xneg_conn[i][count-1] = k;
        }
        //
        if (delt <= 0.00)
        {
            point.ypos_nbhs[i] = point.ypos_nbhs[i] + 1;
            count = point.ypos_nbhs[i];
            point.ypos_conn[i][count-1] = k;
        }
        //
        if (delt > 0.00)
        {
            point.yneg_nbhs[i] = point.yneg_nbhs[i] + 1;
            count = point.yneg_nbhs[i];
            point.yneg_conn[i][count-1] = k;
        }
        //
        if (deln <= 0.00)
        {
            point.zpos_nbhs[i] = point.zpos_nbhs[i] + 1;
            count = point.zpos_nbhs[i];
            point.zpos_conn[i][count-1] = k;
            // if(i==21){
            //     printf("%d %d %d\n",i,count-1,k);
            // }
        }
        //
        if (deln > 0.00)
        {
            point.zneg_nbhs[i] = point.zneg_nbhs[i] + 1;
            count = point.zneg_nbhs[i];
            point.zneg_conn[i][count-1] = k;
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
    double dels, delt;
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
        k = point.conn[i][r];
        xk = point.x[k];
        yk = point.y[k];
        zk = point.z[k];
        //
        delx = xk - xi;
        dely = yk - yi;
        delz = zk - zi;
        //
        dels = delx * point.tan1[i][0] + dely * point.tan1[i][1] + delz * point.tan1[i][2];
        delt = delx * point.tan2[i][0] + dely * point.tan2[i][1] + delz * point.tan2[i][2];
        // deln = delx * point.nor[i][0] + dely * point.nor[i][1] + delz * point.nor[i][2];
        //
        //
        if (dels <= 0.00)
        {
            point.xpos_nbhs[i] = point.xpos_nbhs[i] + 1;
            count = point.xpos_nbhs[i];
            point.xpos_conn[i][count-1] = k;
        }
        //
        if (dels >= 0.00)
        {
            point.xneg_nbhs[i] = point.xneg_nbhs[i] + 1;
            count = point.xneg_nbhs[i];
            point.xneg_conn[i][count-1] = k;
        }
        //
        if (delt <= 0.00)
        {
            point.ypos_nbhs[i] = point.ypos_nbhs[i] + 1;
            count = point.ypos_nbhs[i];
            point.ypos_conn[i][count-1] = k;
        }
        //
        if (delt >= 0.00)
        {
            point.yneg_nbhs[i] = point.yneg_nbhs[i] + 1;
            count = point.yneg_nbhs[i];
            point.yneg_conn[i][count-1] = k;
        }
        //
        point.zneg_nbhs[i] = point.zneg_nbhs[i] + 1;
        count = point.zneg_nbhs[i];
        point.zneg_conn[i][count-1] = k;
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
    double dels, delt;
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
        k = point.conn[i][r];
        xk = point.x[k];
        yk = point.y[k];
        zk = point.z[k];
        //
        delx = xk - xi;
        dely = yk - yi;
        delz = zk - zi;
        //
        dels = delx * point.tan1[i][0] + dely * point.tan1[i][1] + delz * point.tan1[i][2];
        delt = delx * point.tan2[i][0] + dely * point.tan2[i][1] + delz * point.tan2[i][2];
        // deln = delx * point.nor[i][0] + dely * point.nor[i][1] + delz * point.nor[i][2];
        //
        if (dels <= 0.00)
        {
            point.xpos_nbhs[i] = point.xpos_nbhs[i] + 1;
            count = point.xpos_nbhs[i];
            point.xpos_conn[i][count-1] = k;
        }
        //
        if (dels >= 0.00)
        {
            point.xneg_nbhs[i] = point.xneg_nbhs[i] + 1;
            count = point.xneg_nbhs[i];
            point.xneg_conn[i][count-1] = k;
        }
        //
        if (delt <= 0.00)
        {
            point.ypos_nbhs[i] = point.ypos_nbhs[i] + 1;
            count = point.ypos_nbhs[i];
            point.ypos_conn[i][count-1] = k;
        }
        //
        if (delt >= 0.00)
        {
            point.yneg_nbhs[i] = point.yneg_nbhs[i] + 1;
            count = point.yneg_nbhs[i];
            point.yneg_conn[i][count-1] = k;
        }
        //

        point.zpos_nbhs[i] = point.zpos_nbhs[i] + 1;
        count = point.zpos_nbhs[i];
        // if(i==20){
        //     cout<<count<<endl;
        // }
        // if(&point.zpos_conn[21][0]==&point.zpos_conn[i][count-1]){
        //     cout<<"Issue\n";
        //     cout<<i<<" "<<count-1<<endl;
        // }
        // auto ptr1 = &point.zpos_conn[i][count-1];
        // auto ptr2 = point.zpos_conn[i]
        point.zpos_conn[i][count-1] = k;
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
    // cout<<"After Interior\n";
    // for (int j = 0; j < point.zpos_nbhs[21]; j++)
	// //
	// {
    //     int k = point.zpos_conn[21][j];
    //     cout<<21<<" "<<j<<" "<<k<<endl;
    // }
    //
    for (k = 0; k < wall_points; k++)
    {
        i = wall_points_index[k];
        get_wall_point_neighbours(i);
    }
    // cout<<"After Wall\n";
    // for (int j = 0; j < point.zpos_nbhs[21]; j++)
	// //
	// {
    //     int k = point.zpos_conn[21][j];
    //     cout<<21<<" "<<j<<" "<<k<<endl;
    // }
    //
    for (k = 0; k < outer_points; k++)
    {
        i = outer_points_index[k];
        get_outer_point_neighbours(i);
    }
    //
    // cout<<"After Outer\n";
    // for (int j = 0; j < point.zpos_nbhs[21]; j++)
	// //
	// {
    //     int k = point.zpos_conn[21][j];
    //     cout<<21<<" "<<j<<" "<<k<<endl;
    // }
}
//
//
