#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>

double L = 1;
double v = pow(10,-4);
double T_i = 50;
double T_ic = 100;
double T_f = 50;
double a = 0.2;
double tf = 2500;

int main(){

    double dx = 0.01;
    double dt = 0.1*pow(dx,2)/v;

    int nx = 1/dx + 1;
    int ny = 1/dx + 1;
    double nt = tf/dt + 1;

    std::vector<double> x(nx);
    for (int i = 0; i<nx; ++i){
    x[i] = dx*i;
    }

    std::vector<double> y(ny);
    for (int i = 0; i<ny; ++i){
    y[i] = dx*i;
    }
    
    std::vector<std::vector<std::vector<double>>> T(nx, std::vector<std::vector<double>>(ny, std::vector<double>(nt)));

    //Condiciones de frontera

    for (int k = 0; k<nt; ++k){
        for (int i = 0; i < nx; ++i){
            T[i][0][k] = T_f;
            T[i][ny-1][k] = T_f;
        }
        for (int j = 0; j < ny; ++j){
            T[0][j][k] = T_f;
            T[nx-1][j][k] = T_f;
        }
    }

    //Condiciones iniciales

    for (int i = 0; i < nx; ++i){
        for (int j = 0; j < ny; ++j){
            T[i][j][0] = T_i;
        }
    }

    double ind_02 = 0.2/dx;
    double ind_04 = 0.4/dx;
    double ind_06 = 0.6/dx;

    for (int i = ind_02; i < ind_04+1; ++i){
        for (int j = ind_04; j < ind_06+1; ++j){
            T[i][j][0] = T_ic;
        }
    }

    // Diferencias finitas
    for (int k = 0; k < nt-1; ++k){
        for (int i = 1; i < nx-1; ++i){
            for (int j = 1; j < ny-1; ++j){
                T[i][j][k+1] = T[i][j][k] + v*dt/pow(dx,2)*(T[i+1][j][k]+T[i-1][j][k]+T[i][j+1][k]+T[i][j-1][k]-4*T[i][j][k]);
            }
        }
    }

    std::ofstream outfile0;
    std::ofstream outfile100;
    std::ofstream outfile1000;
    std::ofstream outfile2500;
    std::ofstream outfile_xy;

    double ind_100 = 100/dt;
    double ind_1000 = 1000/dt;
    double ind_2500 = 2500/dt;

    outfile0.open("heat_0.dat");
        for (int i = 0; i < nx; ++i){
            for (int j = 0; j < ny; ++j){
                outfile0 << T[i][j][0] << " ";
            }
            outfile0 << "\n";
        }
    outfile0.close();

    outfile100.open("heat_100.dat");
        for (int i = 0; i < nx; ++i){
            for (int j = 0; j < ny; ++j){
                outfile100 << T[i][j][ind_100] << " ";
            }
            outfile100 << "\n";
        }
    outfile100.close();

    outfile1000.open("heat_1000.dat");
        for (int i = 0; i < nx; ++i){
            for (int j = 0; j < ny; ++j){
                outfile1000 << T[i][j][ind_1000] << " ";
            }
            outfile1000 << "\n";
        }
    outfile1000.close();

    outfile2500.open("heat_2500.dat");
        for (int i = 0; i < nx; ++i){
            for (int j = 0; j < ny; ++j){
                outfile2500 << T[i][j][ind_2500] << " ";
            }
            outfile2500 << "\n";
        }
    outfile2500.close();

    outfile_xy.open("xy.dat");
    for (int i = 0; i < nx; ++i){
        outfile_xy << x[i] << " " << y[i] << "\n";
    }
    outfile_xy.close();

    return 0;

}