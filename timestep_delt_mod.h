#pragma once
//
//
#include "data_structure_mod.h"
//

void timestep_delt()
//
//
{
    //
    int i, k, r;
    double rho, u1, u2, u3, pr;
    double delx, dely, delz, dist;
    double min_delt;
    double mod_u, delt;
    //
    //
    for (i = 0; i < max_points; i++)
    {
        //
        min_delt = 1.00;
        //
        for (r = 0; r < point.nbhs[i]; r++)
        //
        {
            k = point.conn[i][r];
            //
            rho = point.prim[k][0];
            u1 = point.prim[k][1];
            u2 = point.prim[k][2];
            u3 = point.prim[k][3];
            pr = point.prim[k][4];
            //
            delx = point.x[k] - point.x[i];
            dely = point.y[k] - point.y[i];
            delz = point.z[k] - point.z[i];
            //
            dist = sqrt(delx * delx + dely * dely + delz * delz);
            //
            mod_u = sqrt(u1 * u1 + u2 * u2 + u3 * u3);

            delt = CFL * dist / (mod_u + 3.00 * sqrt(pr / rho));
            //
            if (delt < min_delt)
            {
                min_delt = delt;
            }
            //
        }
        //
        point.delt[i] = min_delt;
        //
    }
    //
}
//
__global__ void timestep_delt_cuda(points &point,double CFL)
//
//
{
    //
    int k, r;
    double rho, u1, u2, u3, pr;
    double delx, dely, delz, dist;
    double min_delt;
    double mod_u, delt;
    //
    int bx = blockIdx.x;
    int tx = threadIdx.x;
    int i = bx * blockDim.x + tx;
    if (i < 0 || i >= max_points)
    {
        return;
    }
    //
    //
    min_delt = 1.00;
    //
    for (r = 0; r < point.nbhs[i]; r++)
    //
    {
        k = point.conn[i][r];
        //
        rho = point.prim[k][0];
        u1 = point.prim[k][1];
        u2 = point.prim[k][2];
        u3 = point.prim[k][3];
        pr = point.prim[k][4];
        //
        delx = point.x[k] - point.x[i];
        dely = point.y[k] - point.y[i];
        delz = point.z[k] - point.z[i];
        //
        dist = sqrt(delx * delx + dely * dely + delz * delz);
        //
        mod_u = sqrt(u1 * u1 + u2 * u2 + u3 * u3);

        delt = CFL * dist / (mod_u + 3.00 * sqrt(pr / rho));
        //
        if (delt < min_delt)
        {
            min_delt = delt;
        }
        //
    }
    //
    point.delt[i] = min_delt;
    //
    //
}
//
