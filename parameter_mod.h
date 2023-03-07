/*  
	MFCFD is a 3D Computational Fluid Dynamics Solver based off q-LSKUM
    Copyright (C) 2022 
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
#include "math.h"
#include <nccl.h>


		const int max_points = 256000000;
		const int max_partitions=8;
		int max_iters = 1000;//1000
//
//		Flow conditions ..
//
		double Mach = 0.63;
		double aoa = 2.00;
//
//		Other parameters ..
//
		double gamma_new = 1.40;
		double pi=4.00*atan(1.00);
		double theta = aoa*pi/180.00;
		double power = 2.00;
		double VL_CONST = 50.00;
		double CFL = 0.2;
		int inner_iterations = 0;
//
//
//		Freestream values of the primitive variables ..
//
		double rho_inf = 1.00;
		double u1_inf = Mach*cos(theta);
		double u2_inf = 0.00;
		double u3_inf = Mach*sin(theta);
		double pr_inf = 1.00/1.40;

		
		int numDevices = 1;
		int partVector[max_points];
		long long numberOfPointsPerDevice;
		int local_points=0;

//      Variables Required for NCCL

		ncclUniqueId id;
//
//