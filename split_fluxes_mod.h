#pragma once
#include "parameter_mod.h"
#include<iostream>
using namespace std;
//
//
//
void flux_Gxp(double *Gxp, double *t1, double *t2, double *n, double *prim)
//
//
{
	//
	double u1, u2, u3, rho, pr;
	double ut1, ut2, un;
	double beta, S1, B1, A1pos;
	double temp1, temp2, u_sqr;
	//
	//
	rho = prim[0];
	u1 = prim[1];
	u2 = prim[2];
	u3 = prim[3];
	pr = prim[4];
	//
	ut1 = t1[0] * u1 + t1[1] * u2 + t1[2] * u3;
	ut2 = t2[0] * u1 + t2[1] * u2 + t2[2] * u3;
	un = n[0] * u1 + n[1] * u2 + n[2] * u3;
	//
	//
	beta = 0.50 * rho / pr;
	S1 = ut1 * sqrt(beta);
	B1 = 0.50 * exp(-S1 * S1) / sqrt(pi * beta);
	A1pos = 0.50 * (1.00 + erf(S1));
	//
	u_sqr = ut1 * ut1 + ut2 * ut2 + un * un;
	//
	//     Expressions for the split fluxes ..
	//
	temp1 = (ut1 * A1pos + B1);
	Gxp[0] = rho * temp1;
	//
	Gxp[1] = pr * A1pos + rho * ut1 * temp1;
	//
	Gxp[2] = rho * ut2 * temp1;
	//
	Gxp[3] = rho * un * temp1;
	//
	temp1 = 2.50 * pr + 0.50 * rho * u_sqr;
	//		temp1 = 2.50*pr
	Gxp[4] = (temp1 + pr) * ut1 * A1pos + (temp1 + 0.50 * pr) * B1;
	//
	//		Gxp[4] = (3.00*pr + 0.50*rho*u_sqr)*B1
	//
}
//
//
//
void flux_Gxn(double *Gxn, double *t1, double *t2, double *n, double *prim)
//
//
{
	//
	double u1, u2, u3, rho, pr;
	double ut1, ut2, un;
	double beta, S1, B1, A1neg;
	double temp1, temp2, u_sqr;
	//
	//
	rho = prim[0];
	u1 = prim[1];
	u2 = prim[2];
	u3 = prim[3];
	pr = prim[4];
	//
	//
	ut1 = t1[0] * u1 + t1[1] * u2 + t1[2] * u3;
	ut2 = t2[0] * u1 + t2[1] * u2 + t2[2] * u3;
	un = n[0] * u1 + n[1] * u2 + n[2] * u3;
	//
	beta = 0.5 * rho / pr;
	S1 = ut1 * sqrt(beta);
	B1 = 0.5 * exp(-S1 * S1) / sqrt(pi * beta);
	A1neg = 0.5 * (1 - erf(S1));
	//
	u_sqr = ut1 * ut1 + ut2 * ut2 + un * un;
	//
	//     Expressions for the split fluxes ..
	//
	temp1 = (ut1 * A1neg - B1);
	Gxn[0] = rho * temp1;
	//
	Gxn[1] = pr * A1neg + rho * ut1 * temp1;
	//
	Gxn[2] = rho * ut2 * temp1;
	//
	Gxn[3] = rho * un * temp1;
	//
	temp1 = 2.50 * pr + 0.5 * rho * u_sqr;
	Gxn[4] = (temp1 + pr) * ut1 * A1neg - (temp1 + 0.5 * pr) * B1;
	//
	//
}
//
//
void flux_Gyp(double *Gyp, double *t1, double *t2, double *n, double *prim)
//
//
{
	//
	double u1, u2, u3, rho, pr;
	double ut1, ut2, un;
	double beta, S2, B2, A2pos;
	double temp1, temp2, u_sqr;
	//
	//
	rho = prim[0];
	u1 = prim[1];
	u2 = prim[2];
	u3 = prim[3];
	pr = prim[4];
	//
	ut1 = t1[0] * u1 + t1[1] * u2 + t1[2] * u3;
	ut2 = t2[0] * u1 + t2[1] * u2 + t2[2] * u3;
	un = n[0] * u1 + n[1] * u2 + n[2] * u3;
	//
	//
	beta = 0.5 * rho / pr;
	S2 = ut2 * sqrt(beta);
	B2 = 0.5 * exp(-S2 * S2) / sqrt(pi * beta);
	A2pos = 0.5 * (1 + erf(S2));
	//
	u_sqr = ut1 * ut1 + ut2 * ut2 + un * un;
	//
	//     Expressions for the split fluxes ..
	//
	temp1 = (ut2 * A2pos + B2);
	Gyp[0] = rho * temp1;
	//
	Gyp[1] = rho * ut1 * temp1;
	//
	Gyp[2] = pr * A2pos + rho * ut2 * temp1;
	//
	Gyp[3] = rho * un * temp1;
	//
	temp1 = 2.50 * pr + 0.5 * rho * u_sqr;
	Gyp[4] = (temp1 + pr) * ut2 * A2pos + (temp1 + 0.5 * pr) * B2;
	//
	//
}
//
//
void flux_Gyn(double *Gyn, double *t1, double *t2, double *n, double *prim)
//
//
{
	//
	double u1, u2, u3, rho, pr;
	double ut1, ut2, un;
	double beta, S2, B2, A2neg;
	double temp1, temp2, u_sqr;
	//
	//
	rho = prim[0];
	u1 = prim[1];
	u2 = prim[2];
	u3 = prim[3];
	pr = prim[4];
	//
	ut1 = t1[0] * u1 + t1[1] * u2 + t1[2] * u3;
	ut2 = t2[0] * u1 + t2[1] * u2 + t2[2] * u3;
	;
	un = n[0] * u1 + n[1] * u2 + n[2] * u3;
	//
	//
	beta = 0.5 * rho / pr;
	S2 = ut2 * sqrt(beta);
	B2 = 0.5 * exp(-S2 * S2) / sqrt(pi * beta);
	A2neg = 0.5 * (1 - erf(S2));
	//
	u_sqr = ut1 * ut1 + ut2 * ut2 + un * un;
	//
	//     Expressions for the split fluxes ..
	//
	temp1 = (ut2 * A2neg - B2);
	Gyn[0] = rho * temp1;
	//
	Gyn[1] = rho * ut1 * temp1;
	//
	Gyn[2] = pr * A2neg + rho * ut2 * temp1;
	//
	Gyn[3] = rho * un * temp1;
	//
	temp1 = 2.50 * pr + 0.5 * rho * u_sqr;
	Gyn[4] = (temp1 + pr) * ut2 * A2neg - (temp1 + 0.5 * pr) * B2;
	//
	//
}
//
//
void flux_Gzp(double *Gzp, double *t1, double *t2, double *n, double *prim)
//
//
{
	//
	double u1, u2, u3, rho, pr;
	double ut1, ut2, un;
	double beta, S3, B3, A3pos;
	double temp1, temp2, u_sqr;
	//
	//
	rho = prim[0];
	u1 = prim[1];
	u2 = prim[2];
	u3 = prim[3];
	pr = prim[4];
	//
	ut1 = t1[0] * u1 + t1[1] * u2 + t1[2] * u3;
	ut2 = t2[0] * u1 + t2[1] * u2 + t2[2] * u3;
	un = n[0] * u1 + n[1] * u2 + n[2] * u3;
	//
	//
	beta = 0.5 * rho / pr;
	S3 = un * sqrt(beta);
	B3 = 0.5 * exp(-S3 * S3) / sqrt(pi * beta);
	A3pos = 0.5 * (1 + erf(S3));
	//
	u_sqr = ut1 * ut1 + ut2 * ut2 + un * un;
	//
	//     Expressions for the split fluxes ..
	//
	temp1 = (un * A3pos + B3);
	Gzp[0] = rho * temp1;
	//
	Gzp[1] = rho * ut1 * temp1;
	//
	Gzp[2] = rho * ut2 * temp1;
	//
	Gzp[3] = pr * A3pos + rho * un * temp1;
	//
	//
	temp1 = 2.50 * pr + 0.5 * rho * u_sqr;
	Gzp[4] = (temp1 + pr) * un * A3pos + (temp1 + 0.5 * pr) * B3;
	//
	//
}
//
//
void flux_Gzn(double *Gzn, double *t1, double *t2, double *n, double *prim)
//
//
{
	//
	double u1, u2, u3, rho, pr;
	double ut1, ut2, un;
	double beta, S3, B3, A3neg;
	double temp1, temp2, u_sqr;
	//
	//
	rho = prim[0];
	u1 = prim[1];
	u2 = prim[2];
	u3 = prim[3];
	pr = prim[4];
	//
	ut1 = t1[0] * u1 + t1[1] * u2 + t1[2] * u3;
	ut2 = t2[0] * u1 + t2[1] * u2 + t2[2] * u3;
	un = n[0] * u1 + n[1] * u2 + n[2] * u3;

	//
	//
	beta = 0.5 * rho / pr;
	S3 = un * sqrt(beta);
	B3 = 0.5 * exp(-S3 * S3) / sqrt(pi * beta);
	A3neg = 0.5 * (1 - erf(S3));
	//
	u_sqr = ut1 * ut1 + ut2 * ut2 + un * un;
	//
	//     Expressions for the split fluxes ..
	//
	temp1 = (un * A3neg - B3);
	Gzn[0] = rho * temp1;
	//
	Gzn[1] = rho * ut1 * temp1;
	//
	Gzn[2] = rho * ut2 * temp1;
	//
	Gzn[3] = pr * A3neg + rho * un * temp1;
	//
	//
	temp1 = 2.50 * pr + 0.5 * rho * u_sqr;
	Gzn[4] = (temp1 + pr) * un * A3neg - (temp1 + 0.5 * pr) * B3;
	//
	//
}
//
//
//
