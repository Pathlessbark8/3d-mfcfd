#pragma once
//
//	First written on 05.02.2021
//
//
#include "data_structure_mod.h"
#include "generate_connectivity_mod.h"
#include "fpi_solver_mod.h"
//
//
//
//
void q_lskum()
//
{
	//
	int t;
	//
	//
	fstream fout;
	fout.open("residue_cpp.dat", ios::out);
	//
	//
	generate_split_stencils();
	//
	for (t = 1; t <= max_iters; t++)
	{
		fpi_solver(t);
		cout << setprecision(13) << t << " " << res_new << " " << residue << endl;
		fout << setprecision(13) << t << " " << res_new << " " << residue << endl;
	}
	fout.close();
	//
	//
}
//
//
//
