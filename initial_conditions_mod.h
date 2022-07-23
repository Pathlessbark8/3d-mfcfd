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
	int k;
	//
	// fstream fin;
	// fin.open("primal-solution.dat", ios::in);
	//
	for (k = 0; k < max_points; k++)
	{
		point.prim[k][0] = rho_inf;
		point.prim[k][1] = u1_inf;
		point.prim[k][2] = u2_inf;
		point.prim[k][3] = u3_inf;
		point.prim[k][4] = pr_inf;
		// for (int r = 0; r < 5; r++)
		// {
		// 	point.delUp[k][r] = 0.00;
		// 	point.delUn[k][r] = 0.00;
		// }
		// fin >> x1 >> y1 >> z1 >> point.prim[k][0] >> point.prim[k][1] >> point.prim[k][2] >> point.prim[k][3] >> point.prim[k][4];
	}
	// fin.close();
	//
}
//
//
