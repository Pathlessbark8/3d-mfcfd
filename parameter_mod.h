#pragma once
#include "math.h"
#include <nccl.h>


		const int max_points =  580485;
		const int max_devices=	32;
		int max_iters = 1;//1000
//
//		Flow conditions ..
//
		double Mach = 0.7;
		double aoa = 0.00;
//
//		Other parameters ..
//
		double gamma_new = 1.40;
		double pi=4.00*atan(1.00);
		double theta = aoa*pi/180.00;
		double power = 0.00;
		double VL_CONST = 2.00;
		double CFL = 0.5;
		int inner_iterations = 0;
//
//
//		Freestream values of the primitive variables ..
//
		double rho_inf = 1.00;
		double u1_inf = Mach*cos(theta);
		double u2_inf = Mach*sin(theta);
		double u3_inf = 0.00;
		double pr_inf = 1.00/1.40;

		int numDevices = 1;
		int partVector[max_points];
		int numberOfPointsPerDevice[max_devices];

//      Variables Required for NCCL

		ncclComm_t *comms;
		int *devs;
		ncclUniqueId id;

//
//