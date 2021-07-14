
#pragma once
#include "data_structure_mod.h"

//	The following subroutine computes the Venkatakrishnan limiter ..
//
//
void venkat_limiter(double *qtilde, double *phi, int k)
//
//
{
	//
	int r;
	double q, del_neg, del_pos;
	double epsi, num, den, temp;
	//
	for (r = 0; r < 5; r++)
	{
		q = point.q[r][k];
		del_neg = qtilde[r] - q;
		//
		if (abs(del_neg) <= 10e-6)
			phi[r] = 1.00;
		else if (abs(del_neg) > 10e-6)
		{
			if (del_neg > 0.00)
				del_pos = point.qm[0][r][k] - q;
			else if (del_neg < 0.00)
				del_pos = point.qm[1][r][k] - q;
			//
			epsi = VL_CONST * point.min_dist[k];
			epsi = pow(epsi, 3.00);
			num = (del_pos * del_pos) + (epsi * epsi); // Numerator ..
			num = num * del_neg + 2.00 * del_neg * del_neg * del_pos;
			den = del_pos * del_pos + 2.00 * del_neg * del_neg; // Denominator ..
			den = den + del_neg * del_pos + epsi * epsi;
			den = den * del_neg;
			temp = num / den;
			//
			if (temp < 1.0)
				phi[r] = temp;
			else
				phi[r] = 1.0;
		}
	}
	//
}
//



__device__ void venkat_limiter_cuda(points &point,double *qtilde, double *phi, int k,int VL_CONST)
//
//
{
	//
	int r;
	double q, del_neg, del_pos;
	double epsi, num, den, temp;
	//
	for (r = 0; r < 5; r++)
	{
		q = point.q[r][k];
		del_neg = qtilde[r] - q;
		//
		if (abs(del_neg) <= 10e-6)
			phi[r] = 1.00;
		else if (abs(del_neg) > 10e-6)
		{
			if (del_neg > 0.00)
				del_pos = point.qm[0][r][k] - q;
			else if (del_neg < 0.00)
				del_pos = point.qm[1][r][k] - q;
			//
			epsi = VL_CONST * point.min_dist[k];
			epsi = pow(epsi, 3.00);
			num = (del_pos * del_pos) + (epsi * epsi); // Numerator ..
			num = num * del_neg + 2.00 * del_neg * del_neg * del_pos;
			den = del_pos * del_pos + 2.00 * del_neg * del_neg; // Denominator ..
			den = den + del_neg * del_pos + epsi * epsi;
			den = den * del_neg;
			temp = num / den;
			//
			if (temp < 1.0)
				phi[r] = temp;
			else
				phi[r] = 1.0;
		}
	}
	//
}
