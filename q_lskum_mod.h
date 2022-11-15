#pragma once
//
//	First written on 05.02.2021
//
//
#include "data_structure_mod.h"
#include "generate_connectivity_mod.h"
#include "fpi_solver_mod.h"
#include <iomanip>
#include <chrono>

//
//
//
//
using namespace std::chrono;
void q_lskum()
//
{
	//
	int t;
	//
	// fstream fout;
	// fout.open("/home/anil/new_3d_code/3d-mfcfd/temp/inner_3.dat", ios::out);
	//
	//
	generate_split_stencils();
	//
	auto start = high_resolution_clock::now();
	for (t = 1; t <= max_iters; t++)
	{
		fpi_solver(t);
		cout <<scientific<< setprecision(13) << t << " " << res_new << " " << residue << endl;
		// fout << scientific<<setprecision(13) << t << " " << res_new << " " << residue << endl;
	}
	auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    cout << "Time Taken :" << duration.count() / 1000000.0 << endl;
	// fout << "Time Taken :" << duration.count() / 1000000.0 << endl;
	// fout.close();
	//
	//
}
//
//
//
