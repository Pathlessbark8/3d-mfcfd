//
//
#pragma once
#include "data_structure_mod.h"
//
//
//
void qtilde_to_primitive(double *q, double *prim)
//
{
    //
    double beta, temp, temp1, temp2;
    //
    //
    beta = -q[4] * 0.50;
    //
    temp = 0.50 / beta;
    //
    prim[1] = q[1] * temp;
    prim[2] = q[2] * temp;
    prim[3] = q[3] * temp;
    //
    temp1 = q[0] + beta * (prim[1] * prim[1] + prim[2] * prim[2] + prim[3] * prim[3]);
    temp2 = temp1 - (log(beta) * 2.50);
    //
    prim[0] = exp(temp2);
    prim[4] = prim[0] * temp;
}
//
//
//
