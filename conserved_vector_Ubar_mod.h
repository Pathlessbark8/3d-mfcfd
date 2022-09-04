
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
    u1_inf_rot = u1_inf * point.tan1[0][k] + u2_inf * point.tan1[1][k] + u3_inf * point.tan1[2][k];
    u2_inf_rot = u1_inf * point.tan2[0][k] + u2_inf * point.tan2[1][k] + u3_inf * point.tan2[2][k];
    u3_inf_rot = u1_inf * point.nor[0][k] + u2_inf * point.nor[1][k] + u3_inf * point.nor[2][k];
    //
    rho_e_inf = pr_inf * 2.50 + 0.5 * rho_inf * (u1_inf_rot * u1_inf_rot + u2_inf_rot * u2_inf_rot + u3_inf_rot * u3_inf_rot);
    //
    beta = (0.5 * rho_inf) / pr_inf;
    S3 = u3_inf_rot * sqrt(beta);
    B3_inf = exp(-S3 * S3) / (2 * sqrt(pi * beta));
    A3n_inf = 0.5 * (1 - erf(S3));
    //
    rho = point.prim[0][k];
    u1 = point.prim[1][k];
    u2 = point.prim[2][k];
    u3 = point.prim[3][k];
    pr = point.prim[4][k];
    //
    u1_rot = u1 * point.tan1[0][k] + u2 * point.tan1[1][k] + u3 * point.tan1[2][k];
    u2_rot = u1 * point.tan2[0][k] + u2 * point.tan2[1][k] + u3 * point.tan2[2][k];
    u3_rot = u1 * point.nor[0][k] + u2 * point.nor[1][k] + u3 * point.nor[2][k];
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
    u1_inf_rot = u1_inf * point.tan1[0][k] + u2_inf * point.tan1[1][k] + u3_inf * point.tan1[2][k];
    u2_inf_rot = u1_inf * point.tan2[0][k] + u2_inf * point.tan2[1][k] + u3_inf * point.tan2[2][k];
    u3_inf_rot = u1_inf * point.nor[0][k] + u2_inf * point.nor[1][k] + u3_inf * point.nor[2][k];
    //
    rho_e_inf = pr_inf * 2.50 + 0.5 * rho_inf * (u1_inf_rot * u1_inf_rot + u2_inf_rot * u2_inf_rot + u3_inf_rot * u3_inf_rot);
    //
    beta = (0.5 * rho_inf) / pr_inf;
    S3 = u3_inf_rot * sqrt(beta);
    B3_inf = exp(-S3 * S3) / (2 * sqrt(pi * beta));
    A3n_inf = 0.5 * (1 - erf(S3));
    //
    rho = point.prim[0][k];
    u1 = point.prim[1][k];
    u2 = point.prim[2][k];
    u3 = point.prim[3][k];
    pr = point.prim[4][k];
    //
    u1_rot = u1 * point.tan1[0][k] + u2 * point.tan1[1][k] + u3 * point.tan1[2][k];
    u2_rot = u1 * point.tan2[0][k] + u2 * point.tan2[1][k] + u3 * point.tan2[2][k];
    u3_rot = u1 * point.nor[0][k] + u2 * point.nor[1][k] + u3 * point.nor[2][k];
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

__device__ void conserved_vector_Ubar_multi_cuda(splitPoints *splitPoint,int k, double *Ubar,double u1_inf,double u2_inf,double u3_inf,double rho_inf,double pi,double pr_inf)
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
    u1_inf_rot = u1_inf * splitPoint[k].tan1[0] + u2_inf * splitPoint[k].tan1[1] + u3_inf * splitPoint[k].tan1[2];
    u2_inf_rot = u1_inf * splitPoint[k].tan2[0] + u2_inf * splitPoint[k].tan2[1] + u3_inf * splitPoint[k].tan2[2];
    u3_inf_rot = u1_inf * splitPoint[k].nor[0] + u2_inf * splitPoint[k].nor[1] + u3_inf * splitPoint[k].nor[2];
    //
    rho_e_inf = pr_inf * 2.50 + 0.5 * rho_inf * (u1_inf_rot * u1_inf_rot + u2_inf_rot * u2_inf_rot + u3_inf_rot * u3_inf_rot);
    //
    beta = (0.5 * rho_inf) / pr_inf;
    S3 = u3_inf_rot * sqrt(beta);
    B3_inf = exp(-S3 * S3) / (2 * sqrt(pi * beta));
    A3n_inf = 0.5 * (1 - erf(S3));
    //
    rho = splitPoint[k].prim[0];
    u1 = splitPoint[k].prim[1];
    u2 = splitPoint[k].prim[2];
    u3 = splitPoint[k].prim[3];
    pr = splitPoint[k].prim[4];
    //
    u1_rot = u1 * splitPoint[k].tan1[0] + u2 * splitPoint[k].tan1[1] + u3 * splitPoint[k].tan1[2];
    u2_rot = u1 * splitPoint[k].tan2[0] + u2 * splitPoint[k].tan2[1] + u3 * splitPoint[k].tan2[2];
    u3_rot = u1 * splitPoint[k].nor[0] + u2 * splitPoint[k].nor[1] + u3 * splitPoint[k].nor[2];
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