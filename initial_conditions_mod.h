#pragma once
//
#include "data_structure_mod.h"
//

//
//
void initial_conditions()
//
{
	//
	int k, i;
	double U[5], x1, y1, z1;
	//
	// OPEN(UNIT=101,FILE="stored-solution.dat",FORM="FORMATTED",STATUS="OLD",ACTION="READ")
	fstream fin;
	fin.open("/home/dhruv/primal-solution.dat", ios::in);
	//
	for (k = 0; k < max_points; k++)
	{
		// point.prim[0][k] = rho_inf;
		// point.prim[1][k] = u1_inf;
		// point.prim[2][k] = u2_inf;
		// point.prim[3][k] = u3_inf;
		// point.prim[4][k] = pr_inf;
		// for(int r=0;r<5;r++)
		// {
		// point.delUp[r][k] = 0.00;
		// point.delUn[r][k] = 0.00;
		// }
		//				read(101,*) point.prim(1,k), point.prim(2,k), point.prim(3,k), point.prim(4,k), point.prim(5,k)
		fin >> x1 >> y1 >> z1 >> point.prim[0][k] >> point.prim[1][k] >> point.prim[2][k] >> point.prim[3][k] >> point.prim[4][k];
	}
	fin.close();
	//
}
//
//
