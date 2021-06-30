#include <iostream>
#include <fstream>
#include <random>
#include "split_fluxes_mod.h"
#include <iomanip>
#include "octant_fluxes_mod.h"
#include "wall_flux_dGxneg_mod.h"
#include "wall_flux_dGyneg_mod.h"
#include "wall_flux_dGxpos_mod.h"
#include "wall_flux_dGypos_mod.h"
#include "wall_flux_dGzneg_mod.h"
#include "point_preprocessor_mod.h"
#include "compute_conserved_vector_mod.h"
#include "timestep_delt_mod.h"
#include "generate_connectivity_mod.h"
#include "implicit_aliasing_mod.h"
#include "flux_residual_mod.h"
#include "initial_conditions_mod.h"
#include "interior_flux_dGxneg_mod.h"
#include "interior_flux_dGxpos_mod.h"
#include "interior_flux_dGyneg_mod.h"
#include "interior_flux_dGypos_mod.h"
using namespace std;
//  int main(){
// srand(time(0));

// fstream fout;
// fout.open("test_data.dat",ios::out);
// for(int i=0;i<1000;i++)
// {
//     for(int j=0;j<5;j++)
//     {
//         fout<<rand()%100 + 1<<" ";
//     }
//     fout<<endl;
// }
// fout.close();
// }

int main()
{
    // fstream fin,fout;
    // double G[5];
    // fin.open("test_data.dat",ios::in);
    // fout.open("cpp_out.dat",ios::out);
    // fout<<fixed<<scientific;
    // fout<<setprecision(15);
    // for(int i=0;i<1000;i++)
    // {
    //     for(int j=0;j<5;j++)
    //     {
    //         fin>>G[j];
    //     }
    //     wall_dGx_neg(G,i);
    //      wall_dGx_pos(G,i);
    //      wall_dGy_neg(G,i);
    //      wall_dGy_pos(G,i);
    //      wall_dGz_neg(G,i);
    //     for(int j=0;j<5;j++)
    //     {
    //         fout<<G[j]<<" ";
    //     }
    //     fout<<endl;
    // }
    // fin.close();
    // fout.close();
    double G[5] = {0};
    read_input_point_data();
    // cout<<outer_points<<endl;
    initial_conditions();
    // cout<< wall_points<<endl;
    generate_split_stencils();
    aliasing();
    compute_conserved_vector();
    eval_q_variables();
    //  for(int r=0;r<5;r++)
    // {
    //      cout<<point.dq[0][r][0]<<" ";
    // }
    // cout<<endl;
    eval_q_derivatives();
    //  for(int r=0;r<5;r++)
    // {
    //     cout<<point.dq[0][r][0]<<" ";
    // }
    // cout<<endl;
 
    timestep_delt();
    
    eval_flux_residual();
    // fstream fout;
    // fout.open("flux_res_var_cpp",ios::out);
    // for(int i=0;i<wall_points;i++)
    // {
    //     int k = wall_points_index[i];
    //     for(int j=0;j<5;j++)
    //     {
    //         fout<<std::scientific<<point.flux_res[j][k]<<" ";
    //     }
    //     fout<<endl;
    // }
    // fout.close();
    // for (int i = 0; i < wall_points; i++)
    // {
    //     int k = wall_points_index[i];
        // for (int r = 0; r < 5; r++)
        // {
        //     cout << point.dq[0][r][k] << " ";
        // }
        // cout << endl;
        // wall_dGx_neg(G,k);
        // wall_dGx_pos(G,k);
        // wall_dGy_neg(G,k);
        // wall_dGy_pos(G,k);
        // wall_dGz_neg(G,k);
        // for (int j = 0; j < 5; j++)
        // {
        //     cout << G[j] << " ";
        // }
        // cout << endl;
    // }
    // for (int i = 0; i < interior_points; i++)
    // {
    //     int k = interior_points_index[i];
    //     // // cout<<"k "<<k<<endl;
    //     interior_dGx_neg(G, k);
    //     for (int j = 0; j < 5; j++)
    //     {
    //         cout << G[j] << " ";
    //     }
    //     cout << endl;
    // }
}