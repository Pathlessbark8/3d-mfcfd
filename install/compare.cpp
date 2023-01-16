#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <iomanip>

using namespace std;

int main(){
    int max_points=1000000;
    fstream fin1,fin2;
    fin1>>setprecision(16);
    fin2>>setprecision(16);
    fin1.open("/home/nsm/3d-mfcfd/install/output_fluxres.dat",ios::in);
    fin2.open("/home/nsm/3d-mfcfd/install/output_fluxres_serial.dat",ios::in);
    cout<<setprecision(16);
    double prim1[5] = {};
    double prim2[5] = {};
    for(int i=0;i<max_points;++i){
        for(int r=0;r<5;++r){
            fin1>>prim1[r];
            fin2>>prim2[r];
            // cout<<prim1[r]<<" "<<prim2[r]<<endl;
            if(abs(prim1[r]-prim2[r])<1e-11){
                
            }
            else{
                cout<<"Values don't match for: "<<endl;
                cout<<i<<" "<<r<<" "<<prim1[r]<<" "<<prim2[r]<<endl;
                // exit(0);
            }
        }
    }
    fin1.close();
    fin2.close();
    return 0;
}