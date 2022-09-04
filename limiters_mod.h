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


__device__ void venkat_limiter_multi_nccl(splitPoints *splitPoint,double *qtilde, double *phi, int k,int VL_CONST)
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
		q = splitPoint[k].q[r];
		del_neg = qtilde[r] - q;
		//
		if (abs(del_neg) <= 10e-6)
			phi[r] = 1.00;
		else if (abs(del_neg) > 10e-6)
		{
			if (del_neg > 0.00)
				del_pos = splitPoint[k].qm[0][r] - q;
			else if (del_neg < 0.00)
				del_pos = splitPoint[k].qm[1][r] - q;
			//
			epsi = VL_CONST * splitPoint[k].min_dist;
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

__device__ void venkat_limiter_multi_nccl_ghost(transferPoints **receiveBuffer,double *qtilde, double *phi, int device,int ghostIndex, int VL_CONST)
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
		q = receiveBuffer[device][ghostIndex].q[r];
		del_neg = qtilde[r] - q;
		//
		if (abs(del_neg) <= 10e-6)
			phi[r] = 1.00;
		else if (abs(del_neg) > 10e-6)
		{
			if (del_neg > 0.00)
				del_pos = receiveBuffer[device][ghostIndex].qm[0][r] - q;
			else if (del_neg < 0.00)
				del_pos = receiveBuffer[device][ghostIndex].qm[1][r] - q;
			//
			epsi = VL_CONST * receiveBuffer[device][ghostIndex].min_dist;
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