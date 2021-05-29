#pragma once
#include "parameter_mod.h"
//
//
void flux_Gwxn(double *G, double *t1, double *t2, double *n, double *prim)
//
//
{
	//

	double u1, u2, u3, rho, pr;

	double ut1, ut2, un;
	double beta, S1, B1, S3, B3, A1neg, A3neg;
	double temp1, temp2, temp3, u_sqr;
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
	S1 = ut1 * sqrt(beta);
	S3 = un * sqrt(beta);
	B1 = 0.5 * exp(-S1 * S1) / sqrt(pi * beta);
	B3 = 0.5 * exp(-S3 * S3) / sqrt(pi * beta);
	A1neg = 0.5 * (1 - erf(S1));
	A3neg = 0.5 * (1 - erf(S3));
	//
	u_sqr = ut1 * ut1 + ut2 * ut2 + un * un;
	//
	//     Expressions for the split fluxes ..
	//
	temp1 = (ut1 * A1neg - B1);
	G[0] = rho * temp1 * A3neg;
	//
	G[1] = (pr * A1neg + rho * ut1 * temp1) * A3neg;
	//
	G[2] = rho * ut2 * temp1 * A3neg;
	//
	G[3] = rho * temp1 * (un * A3neg - B3);
	//
	temp2 = 2.50 * pr + 0.5 * rho * u_sqr;

	temp3 = (temp2 + pr) * ut1 * A1neg - (temp2 + 0.5 * pr) * B1;

	G[4] = temp3 * A3neg - 0.5 * rho * un * B3 * temp1;
	//
	//
}
//
//
void flux_Gwxp(double *G, double *t1, double *t2, double *n, double *prim)
//
//
{
	//

	double u1, u2, u3, rho, pr;

	double ut1, ut2, un;
	double beta, S1, B1, S3, B3, A1pos, A3neg;
	double temp1, temp2, temp3, u_sqr;
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
	S1 = ut1 * sqrt(beta);
	S3 = un * sqrt(beta);
	B1 = 0.5 * exp(-S1 * S1) / sqrt(pi * beta);
	B3 = 0.5 * exp(-S3 * S3) / sqrt(pi * beta);
	A1pos = 0.5 * (1 + erf(S1));
	A3neg = 0.5 * (1 - erf(S3));
	//
	u_sqr = ut1 * ut1 + ut2 * ut2 + un * un;
	//
	//     Expressions for the split fluxes ..
	//
	temp1 = (ut1 * A1pos + B1);
	G[0] = rho * temp1 * A3neg;
	//
	G[1] = (pr * A1pos + rho * ut1 * temp1) * A3neg;
	//
	G[2] = rho * ut2 * temp1 * A3neg;
	//
	G[3] = rho * temp1 * (un * A3neg - B3);
	//
	temp2 = 2.50 * pr + 0.5 * rho * u_sqr;

	temp3 = (temp2 + pr) * ut1 * A1pos + (temp2 + 0.5 * pr) * B1;

	G[4] = temp3 * A3neg - 0.5 * rho * un * B3 * temp1;
	//
	//
}
//
//
void flux_Goxp(double *G, double *t1, double *t2, double *n, double *prim)
//
//
{
	//

	double u1, u2, u3, rho, pr;

	double ut1, ut2, un;
	double beta, S1, B1, S3, B3, A1pos, A3pos;
	double temp1, temp2, temp3, u_sqr;
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
	beta = 0.5 * rho / pr;
	S1 = ut1 * sqrt(beta);
	S3 = un * sqrt(beta);
	B1 = 0.5 * exp(-S1 * S1) / sqrt(pi * beta);
	B3 = 0.5 * exp(-S3 * S3) / sqrt(pi * beta);
	A1pos = 0.5 * (1 + erf(S1));
	A3pos = 0.5 * (1 + erf(S3));
	//
	u_sqr = ut1 * ut1 + ut2 * ut2 + un * un;
	//
	//     Expressions for the split fluxes ..
	//
	temp1 = (ut1 * A1pos + B1);
	G[0] = rho * temp1 * A3pos;
	//
	G[1] = (pr * A1pos + rho * ut1 * temp1) * A3pos;
	//
	G[2] = rho * ut2 * temp1 * A3pos;
	//
	G[3] = rho * temp1 * (un * A3pos + B3);
	//
	temp2 = 2.50 * pr + 0.5 * rho * u_sqr;

	temp3 = (temp2 + pr) * ut1 * A1pos + (temp2 + 0.5 * pr) * B1;

	G[4] = temp3 * A3pos + 0.5 * rho * un * B3 * temp1;
	//
	//
}
//
//
void flux_Goxn(double *G, double *t1, double *t2, double *n, double *prim)
//
//
{
	//

	double u1, u2, u3, rho, pr;

	double ut1, ut2, un;
	double beta, S1, B1, S3, B3, A1neg, A3pos;
	double temp1, temp2, temp3, u_sqr;
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
	beta = 0.5 * rho / pr;
	S1 = ut1 * sqrt(beta);
	S3 = un * sqrt(beta);
	B1 = 0.5 * exp(-S1 * S1) / sqrt(pi * beta);
	B3 = 0.5 * exp(-S3 * S3) / sqrt(pi * beta);
	A1neg = 0.5 * (1 - erf(S1));
	A3pos = 0.5 * (1 + erf(S3));
	//
	u_sqr = ut1 * ut1 + ut2 * ut2 + un * un;
	//
	//     Expressions for the split fluxes ..
	//
	temp1 = (ut1 * A1neg - B1);
	G[0] = rho * temp1 * A3pos;
	//
	G[1] = (pr * A1neg + rho * ut1 * temp1) * A3pos;
	//
	G[2] = rho * ut2 * temp1 * A3pos;
	//
	G[3] = rho * temp1 * (un * A3pos + B3);
	//
	temp2 = 2.50 * pr + 0.5 * rho * u_sqr;

	temp3 = (temp2 + pr) * ut1 * A1neg - (temp2 + 0.5 * pr) * B1;

	G[4] = temp3 * A3pos + 0.5 * rho * un * B3 * temp1;
	//
	//
}
//
//
void flux_Gwyn(double *G, double *t1, double *t2, double *n, double *prim)
//
//
{
	//

	double u1, u2, u3, rho, pr;

	double ut1, ut2, un;
	double beta, S2, B2, S3, B3, A2neg, A3neg;
	double temp1, temp2, temp3, u_sqr;
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
	beta = 0.5 * rho / pr;
	S2 = ut2 * sqrt(beta);
	S3 = un * sqrt(beta);
	B2 = 0.5 * exp(-S2 * S2) / sqrt(pi * beta);
	B3 = 0.5 * exp(-S3 * S3) / sqrt(pi * beta);
	A2neg = 0.5 * (1 - erf(S2));
	A3neg = 0.5 * (1 - erf(S3));
	//
	u_sqr = ut1 * ut1 + ut2 * ut2 + un * un;
	//
	//     Expressions for the split fluxes ..
	//
	temp1 = (ut2 * A2neg - B2);
	G[0] = rho * temp1 * A3neg;
	//
	G[1] = rho * ut1 * temp1 * A3neg;

	G[2] = (pr * A2neg + rho * ut2 * temp1) * A3neg;
	//
	G[3] = rho * temp1 * (un * A3neg - B3);
	//
	temp2 = 2.50 * pr + 0.5 * rho * u_sqr;

	temp3 = (temp2 + pr) * ut2 * A2neg - (temp2 + 0.5 * pr) * B2;

	G[4] = temp3 * A3neg - 0.5 * rho * un * B3 * temp1;
	//
	//
}
//
//
void flux_Gwyp(double *G, double *t1, double *t2, double *n, double *prim)
//
//
{
	//

	double u1, u2, u3, rho, pr;

	double ut1, ut2, un;
	double beta, S2, B2, S3, B3, A2pos, A3neg;
	double temp1, temp2, temp3, u_sqr;
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
	beta = 0.5 * rho / pr;
	S2 = ut2 * sqrt(beta);
	S3 = un * sqrt(beta);
	B2 = 0.5 * exp(-S2 * S2) / sqrt(pi * beta);
	B3 = 0.5 * exp(-S3 * S3) / sqrt(pi * beta);
	A2pos = 0.5 * (1 + erf(S2));
	A3neg = 0.5 * (1 - erf(S3));
	//
	u_sqr = ut1 * ut1 + ut2 * ut2 + un * un;
	//
	//     Expressions for the split fluxes ..
	//
	temp1 = (ut2 * A2pos + B2);
	G[0] = rho * temp1 * A3neg;
	//
	G[1] = rho * ut1 * temp1 * A3neg;

	G[2] = (pr * A2pos + rho * ut2 * temp1) * A3neg;
	//
	G[3] = rho * temp1 * (un * A3neg - B3);
	//
	temp2 = 2.50 * pr + 0.5 * rho * u_sqr;

	temp3 = (temp2 + pr) * ut2 * A2pos + (temp2 + 0.5 * pr) * B2;

	G[4] = temp3 * A3neg - 0.5 * rho * un * B3 * temp1;
	//
	//
}
//
//
void flux_Goyp(double *G, double *t1, double *t2, double *n, double *prim)
//
//
{
	//

	double u1, u2, u3, rho, pr;

	double ut1, ut2, un;
	double beta, S2, B2, S3, B3, A2pos, A3pos;
	double temp1, temp2, temp3, u_sqr;
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
	beta = 0.5 * rho / pr;
	S2 = ut2 * sqrt(beta);
	S3 = un * sqrt(beta);
	B2 = 0.5 * exp(-S2 * S2) / sqrt(pi * beta);
	B3 = 0.5 * exp(-S3 * S3) / sqrt(pi * beta);
	A2pos = 0.5 * (1 + erf(S2));
	A3pos = 0.5 * (1 + erf(S3));
	//
	u_sqr = ut1 * ut1 + ut2 * ut2 + un * un;
	//
	//     Expressions for the split fluxes ..
	//
	temp1 = (ut2 * A2pos + B2);
	G[0] = rho * temp1 * A3pos;
	//
	G[1] = rho * ut1 * temp1 * A3pos;

	G[2] = (pr * A2pos + rho * ut2 * temp1) * A3pos;
	//
	G[3] = rho * temp1 * (un * A3pos + B3);
	//
	temp2 = 2.50 * pr + 0.5 * rho * u_sqr;

	temp3 = (temp2 + pr) * ut2 * A2pos + (temp2 + 0.5 * pr) * B2;

	G[4] = temp3 * A3pos + 0.5 * rho * un * B3 * temp1;
	//
	//
}
//
//
void flux_Goyn(double *G, double *t1, double *t2, double *n, double *prim)
//
//
{
	//

	double u1, u2, u3, rho, pr;

	double ut1, ut2, un;
	double beta, S2, B2, S3, B3, A2neg, A3pos;
	double temp1, temp2, temp3, u_sqr;
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
	beta = 0.5 * rho / pr;
	S2 = ut2 * sqrt(beta);
	S3 = un * sqrt(beta);
	B2 = 0.5 * exp(-S2 * S2) / sqrt(pi * beta);
	B3 = 0.5 * exp(-S3 * S3) / sqrt(pi * beta);
	A2neg = 0.5 * (1 - erf(S2));
	A3pos = 0.5 * (1 + erf(S3));
	//
	u_sqr = ut1 * ut1 + ut2 * ut2 + un * un;
	//
	//     Expressions for the split fluxes ..
	//
	temp1 = (ut2 * A2neg - B2);
	G[0] = rho * temp1 * A3pos;
	//
	G[1] = rho * ut1 * temp1 * A3pos;

	G[2] = (pr * A2neg + rho * ut2 * temp1) * A3pos;
	//
	G[3] = rho * temp1 * (un * A3pos + B3);
	//
	temp2 = 2.50 * pr + 0.5 * rho * u_sqr;

	temp3 = (temp2 + pr) * ut2 * A2neg - (temp2 + 0.5 * pr) * B2;

	G[4] = temp3 * A3pos + 0.5 * rho * un * B3 * temp1;
	//
	//
}
//
//
