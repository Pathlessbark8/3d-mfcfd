#include<iostream>
#include<fstream>
#include<random>
#include"split_fluxes_mod.h"
#include<iomanip>
#include "octant_fluxes_mod.h"
using namespace std;
// int main(){
// srand(time(0));

// fstream fout;
// fout.open("test_data.dat",ios::out);
// for(int i=0;i<1000;i++)
// {
//     for(int j=0;j<19;j++)
//     {
//         fout<<rand()<<" ";
//     }
//     fout<<endl;
// }
// fout.close();
// }

int main(){
    fstream fin,fout;
    double G[5],t1[3],prim[5],n[3],t2[3];
    fin.open("test_data.dat",ios::in);
    fout.open("cpp_out.dat",ios::out);
    fout<<fixed<<scientific;
    fout<<setprecision(15);
    for(int i=0;i<1000;i++)
    {
        for(int j=0;j<5;j++)
        {
            fin>>G[j];
        }
        for(int j=0;j<3;j++)
        {
            fin>>t1[j];
        }
        for(int j=0;j<3;j++)
        {
            fin>>t2[j];
        }
        for(int j=0;j<3;j++)
        {
            fin>>n[j];
        }
        for(int j=0;j<5;j++)
        {
            fin>>prim[j];
        }
        flux_Gxp(G,t1,t2,n,prim);
        flux_Gxn(G,t1,t2,n,prim);
        flux_Gyn(G,t1,t2,n,prim);
        flux_Gyp(G,t1,t2,n,prim);
        flux_Gzn(G,t1,t2,n,prim);
        flux_Gzp(G,t1,t2,n,prim);
        flux_Gwxn(G,t1,t2,n,prim);
        flux_Gwxp(G,t1,t2,n,prim);
        flux_Gwyn(G,t1,t2,n,prim);
        flux_Gwyp(G,t1,t2,n,prim);
        flux_Goxn(G,t1,t2,n,prim);
        flux_Goxp(G,t1,t2,n,prim);
        flux_Goyn(G,t1,t2,n,prim);
        flux_Goyp(G,t1,t2,n,prim);
        for(int j=0;j<5;j++)
        {
            fout<<G[j]<<" ";
        }
        for(int j=0;j<3;j++)
        {
            fout<<t1[j]<<" ";
        }
        for(int j=0;j<3;j++)
        {
            fout<<t2[j]<<" ";
        }
        for(int j=0;j<3;j++)
        {
            fout<<n[j]<<" ";
        }
        for(int j=0;j<5;j++)
        {
            fout<<prim[j]<<" ";
        }
        fout<<endl;
    }
    fin.close();
    fout.close();
}