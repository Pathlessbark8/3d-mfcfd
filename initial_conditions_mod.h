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
#include "data_structure_mod.h"
//

//
//
void initial_conditions()
//
{
	//
	int k;
	double x1, y1, z1;
	//
	// fstream fin;
	// fin.open("primal-solution.dat", ios::in);
	//
	for (k = 0; k < max_points; k++)
	{
		point.prim[0][k] = rho_inf;
		point.prim[1][k] = u1_inf;
		point.prim[2][k] = u2_inf;
		point.prim[3][k] = u3_inf;
		point.prim[4][k] = pr_inf;
		// for (int r = 0; r < 5; r++)
		// {
		// 	point.delUp[r][k] = 0.00;
		// 	point.delUn[r][k] = 0.00;
		// }
		// fin >> x1 >> y1 >> z1 >> point.prim[0][k] >> point.prim[1][k] >> point.prim[2][k] >> point.prim[3][k] >> point.prim[4][k];
	}
	// fin.close();
	//
}
//
//
