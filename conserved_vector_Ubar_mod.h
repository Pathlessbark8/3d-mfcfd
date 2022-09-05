#pragma once
//
//
#include "data_structure_mod.h"
//

//
void conserved_vector_Ubar(int k, double *Ubar)
//
{
    //
    //
    double u1_inf_rot, u2_inf_rot, u3_inf_rot, rho_e_inf;
    double u1, u2, u3, pr, rho, u1_rot, u2_rot, u3_rot, rho_e;
    double beta, S3, B3_inf, A3n_inf;
    double B3, A3p, temp1, temp2;
    //
    //
    //
    u1_inf_rot = u1_inf * point.tan1[k][0] + u2_inf * point.tan1[k][1] + u3_inf * point.tan1[k][2];
    u2_inf_rot = u1_inf * point.tan2[k][0] + u2_inf * point.tan2[k][1] + u3_inf * point.tan2[k][2];
    u3_inf_rot = u1_inf * point.nor[k][0] + u2_inf * point.nor[k][1] + u3_inf * point.nor[k][2];
    //
    rho_e_inf = pr_inf * 2.50 + 0.5 * rho_inf * (u1_inf_rot * u1_inf_rot + u2_inf_rot * u2_inf_rot + u3_inf_rot * u3_inf_rot);
    //
    beta = (0.5 * rho_inf) / pr_inf;
    S3 = u3_inf_rot * sqrt(beta);
    B3_inf = exp(-S3 * S3) / (2 * sqrt(pi * beta));
    A3n_inf = 0.5 * (1 - erf(S3));
    //
    rho = point.prim[k][0];
    u1 = point.prim[k][1];
    u2 = point.prim[k][2];
    u3 = point.prim[k][3];
    pr = point.prim[k][4];
    //
    u1_rot = u1 * point.tan1[k][0] + u2 * point.tan1[k][1] + u3 * point.tan1[k][2];
    u2_rot = u1 * point.tan2[k][0] + u2 * point.tan2[k][1] + u3 * point.tan2[k][2];
    u3_rot = u1 * point.nor[k][0] + u2 * point.nor[k][1] + u3 * point.nor[k][2];
    //
    rho_e = pr * 2.50 + 0.5 * rho * (u1_rot * u1_rot + u2_rot * u2_rot + u3_rot * u3_rot);
    //
    beta = (0.5 * rho) / pr;
    S3 = u3_rot * sqrt(beta);
    B3 = exp(-S3 * S3) / (2 * sqrt(pi * beta));
    A3p = 0.5 * (1 + erf(S3));
    //
    Ubar[0] = (rho_inf * A3n_inf) + (rho * A3p);
    Ubar[1] = (rho_inf * u1_inf_rot * A3n_inf) + (rho * u1_rot * A3p);
    Ubar[2] = (rho_inf * u2_inf_rot * A3n_inf) + (rho * u2_rot * A3p);
    //
    temp1 = rho_inf * (u3_inf_rot * A3n_inf - B3_inf);
    temp2 = rho * (u3_rot * A3p + B3);
    Ubar[3] = temp1 + temp2;
    //
    temp1 = (rho_e_inf * A3n_inf - 0.5 * rho_inf * u3_inf_rot * B3_inf);
    temp2 = (rho_e * A3p + 0.5 * rho * u3_rot * B3);
    //
    Ubar[4] = temp1 + temp2;
    //
    //
}

//
__device__ void conserved_vector_Ubar_cuda(points &point,int k, double *Ubar,double u1_inf,double u2_inf,double u3_inf,double rho_inf,double pi,double pr_inf)
//
{
    //
    //
    double u1_inf_rot, u2_inf_rot, u3_inf_rot, rho_e_inf;
    double u1, u2, u3, pr, rho, u1_rot, u2_rot, u3_rot, rho_e;
    double beta, S3, B3_inf, A3n_inf;
    double B3, A3p, temp1, temp2;
    //
    //
    //
    u1_inf_rot = u1_inf * point.tan1[k][0] + u2_inf * point.tan1[k][1] + u3_inf * point.tan1[k][2];
    u2_inf_rot = u1_inf * point.tan2[k][0] + u2_inf * point.tan2[k][1] + u3_inf * point.tan2[k][2];
    u3_inf_rot = u1_inf * point.nor[k][0] + u2_inf * point.nor[k][1] + u3_inf * point.nor[k][2];
    //
    rho_e_inf = pr_inf * 2.50 + 0.5 * rho_inf * (u1_inf_rot * u1_inf_rot + u2_inf_rot * u2_inf_rot + u3_inf_rot * u3_inf_rot);
    //
    beta = (0.5 * rho_inf) / pr_inf;
    S3 = u3_inf_rot * sqrt(beta);
    B3_inf = exp(-S3 * S3) / (2 * sqrt(pi * beta));
    A3n_inf = 0.5 * (1 - erf(S3));
    //
    rho = point.prim[k][0];
    u1 = point.prim[k][1];
    u2 = point.prim[k][2];
    u3 = point.prim[k][3];
    pr = point.prim[k][4];
    //
    u1_rot = u1 * point.tan1[k][0] + u2 * point.tan1[k][1] + u3 * point.tan1[k][2];
    u2_rot = u1 * point.tan2[k][0] + u2 * point.tan2[k][1] + u3 * point.tan2[k][2];
    u3_rot = u1 * point.nor[k][0] + u2 * point.nor[k][1] + u3 * point.nor[k][2];
    //
    rho_e = pr * 2.50 + 0.5 * rho * (u1_rot * u1_rot + u2_rot * u2_rot + u3_rot * u3_rot);
    //
    beta = (0.5 * rho) / pr;
    S3 = u3_rot * sqrt(beta);
    B3 = exp(-S3 * S3) / (2 * sqrt(pi * beta));
    A3p = 0.5 * (1 + erf(S3));
    //
    Ubar[0] = (rho_inf * A3n_inf) + (rho * A3p);
    Ubar[1] = (rho_inf * u1_inf_rot * A3n_inf) + (rho * u1_rot * A3p);
    Ubar[2] = (rho_inf * u2_inf_rot * A3n_inf) + (rho * u2_rot * A3p);
    //
    temp1 = rho_inf * (u3_inf_rot * A3n_inf - B3_inf);
    temp2 = rho * (u3_rot * A3p + B3);
    Ubar[3] = temp1 + temp2;
    //
    temp1 = (rho_e_inf * A3n_inf - 0.5 * rho_inf * u3_inf_rot * B3_inf);
    temp2 = (rho_e * A3p + 0.5 * rho * u3_rot * B3);
    //
    Ubar[4] = temp1 + temp2;
    //
    //
}