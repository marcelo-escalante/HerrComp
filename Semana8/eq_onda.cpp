#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>

double fx_inicial(double x);

double L = 2.0;
double c = 300.0;
double h0 = 0.1;
double tf = 0.1;

int main(){


    double a = 0.0;
    double b = L;
    double t0 = 0.0;
    int nx = 100;
    double dx = (b-a)/(nx-1);
    double dt = 0.5*dx/c;
    int nt = (tf-t0)/(dt)+1;
    double cp = dx/dt;

    std::vector<double> x(nx);
    for (int i = 0; i<nx; ++i){
    x[i] = dx*i;
    }

    std::vector<std::vector<double>> u(nx, std::vector<double>(nt));

    // Condiciones de frontera
    for (int j = 0; j < nt; ++j){
        u[0][j] = 0;
        u[nx-1][j] = 0;
    }

    // Condiciones iniciales
    for (int i = 0; i<nx; ++i){
        u[i][0] = fx_inicial(x[i]);
    }

    for (int i = 1; i<nx-1; ++i){
        u[i][1] = u[i][0]+0.5*pow(c/cp,2)*(u[i+1][0]+u[i-1][0]-2*u[i][0]); 
    }

    // Diferencias finitas
    for (int j = 1; j<nt-1; ++j){
        for (int i = 1; i < nx-1; ++i){
          u[i][j+1] = 2*u[i][j]-u[i][j-1]+pow(c/cp,2)*(u[i+1][j]-2*u[i][j]+u[i-1][j]);
        }
    }

    std::ofstream outfile;
    std::ofstream outfile_x;
    outfile.open("wave.dat");
    outfile_x.open("x.dat");
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < nt;j++){
            outfile << u[i][j] << " ";
        }
        outfile << "\n";
        outfile_x << x[i] << " ";
    }
    outfile.close();
    outfile_x.close();

   
    return 0;
}


double fx_inicial(double x){
    if (x<=0.5*L){
        return 2*h0/L*x;
    } else{
        return -2*h0/L*x+2*h0;
    }
}