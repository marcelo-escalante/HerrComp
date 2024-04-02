#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>

void euler(double M, double G, double e, std::vector<double> &t, double x0, double y0, double vx0, double vy0);
void leapfrog(double M, double G, double e, std::vector<double> &t, double x0, double y0, double vx0, double vy0);
void rk4(double M, double G, double e, std::vector<double> &t, double x0, double y0, double vx0, double vy0);

int main(){

  double M = 1.0; 
  double G = 2.95900*pow(10,-4);
  double e = 0.01671; 
  double a = 1.000003;

  double x0 = -0.983;
  double y0 = 0.0;
  double vx0 = 0.0;
  double vy0 = -sqrt(G*M/a*(1+e)/(1-e));

  double start = 0.0;
  double end = 2200.0; //2200 días (seis años y diez días)
  double h = 0.001;
  double num_points = num_points = (end - start) / h + 1;

  std::vector<double> t(num_points);
  for (int i = 0; i<num_points; ++i){
    t[i] = start + h*i;
  }

  euler(M,G,e,t,x0,y0,vx0,vy0);
  leapfrog(M,G,e,t,x0,y0,vx0,vy0);
  rk4(M,G,e,t,x0,y0,vx0,vy0);

  return 0;

}

//Implementación Euler

void euler(double M, double G, double e, std::vector<double> &t, double x0, double y0, double vx0, double vy0){

  const int num_points = t.size();
  std::vector<double> x(num_points);
  std::vector<double> vx(num_points);
  std::vector<double> y(num_points);
  std::vector<double> vy(num_points);
  std::vector<double> R3(num_points);
  double h = (t.back() - t.front()) / (num_points - 1);
  x[0] = x0;
  y[0] = y0;
  vx[0] = vx0;
  vy[0] = vy0;

  for (int i = 0; i < num_points - 1; i++) {
    R3[i] = pow(sqrt(pow(x[i],2)+pow(y[i],2)),3);
    x[i + 1] = x[i] + h * vx[i];
    y[i + 1] = y[i] + h * vy[i];
    vx[i + 1] = vx[i] + h * (-G*M/R3[i])*x[i];
    vy[i + 1] = vy[i] + h * (-G*M/R3[i])*y[i];
  }

  std::ofstream outfile;
  outfile.open("data1_euler.dat");
  for (int i = 0; i < num_points; i++) {
    outfile << x[i] << " " << y[i] << std::endl;
  }
  outfile.close();
}

//Implementación LeapFrog

void leapfrog(double M, double G, double e, std::vector<double> &t, double x0, double y0, double vx0, double vy0){

  const int num_points = t.size();
  std::vector<double> x(num_points);
  std::vector<double> vx(num_points);
  std::vector<double> y(num_points);
  std::vector<double> vy(num_points);
  std::vector<double> R3(num_points);
  double h = (t.back() - t.front()) / (num_points - 1);
  x[0] = x0;
  y[0] = y0;

  R3[0] = pow(sqrt(pow(x0,2)+pow(y0,2)),3);
  vx[0] = vx0 + h * 0.5 * (-G*M/R3[0] * x[0]);
  vy[0] = vy0 + h * 0.5 * (-G*M/R3[0] * y[0]);

  for (int i = 0; i < num_points - 1; i++) {
    x[i + 1] = x[i] + h * vx[i];
    y[i + 1] = y[i] + h * vy[i];
    R3[i+1] = pow(sqrt(pow(x[i+1],2)+pow(y[i+1],2)),3);
    vx[i + 1] = vx[i] + h * (-G*M/R3[i+1])*x[i+1];
    vy[i + 1] = vy[i] + h * (-G*M/R3[i+1])*y[i+1];
  }

  std::ofstream outfile;
  outfile.open("data2_leapfrog.dat");
  for (int i = 0; i < num_points; i++) {
    outfile << x[i] << " " << y[i] << std::endl;
  }
  outfile.close();
}

//Implementación Runge Kutta orden 4

//Funciones auxiliares
double R3(double x, double y);
double dvx(double G, double m, double x, double y);
double dvy(double G, double m, double x, double y);
//

void rk4(double M, double G, double e, std::vector<double> &t, double x0, double y0, double vx0, double vy0) {

  const int num_points = t.size();
  std::vector<double> x(num_points);
  std::vector<double> vx(num_points);
  std::vector<double> y(num_points);
  std::vector<double> vy(num_points);
  std::vector<double> R3(num_points);
  double h = (t.back() - t.front()) / (num_points - 1);
  x[0] = x0;
  y[0] = y0;
  vx[0] = vx0;
  vy[0] = vy0;

  for (int i = 0; i < num_points - 1; ++i) {

    double k1x = h * vx[i];
    double k1y = h * vy[i];
    
    double k1vx = h * dvx(G,M,x[i],y[i]);
    double k1vy = h * dvy(G,M,x[i],y[i]);

    double k2x = h * (vx[i]+k1x/2);
    double k2y = h * (vy[i]+k1y/2);
    
    double k2vx = h * dvx(G,M,x[i]+k1vx/2,y[i]+k1vy/2);
    double k2vy = h * dvy(G,M,x[i]+k1vx/2,y[i]+k1vy/2);

    double k3x = h * (vx[i]+k2x/2);
    double k3y = h * (vy[i]+k2y/2);
    
    double k3vx = h * dvx(G,M,x[i]+k2vx/2,y[i]+k2vy/2);
    double k3vy = h * dvy(G,M,x[i]+k2vx/2,y[i]+k2vy/2);

    double k4x = h * (vx[i]+k3x/2);
    double k4y = h * (vy[i]+k3y/2);
    
    double k4vx = h * dvx(G,M,x[i]+k3vx/2,y[i]+k3vy/2);
    double k4vy = h * dvy(G,M,x[i]+k3vx/2,y[i]+k3vy/2);

    x[i + 1] = x[i] + (k1x + 2 * k2x + 2 * k3x + k4x) / 6.0;
    y[i + 1] = y[i] + (k1y + 2 * k2y + 2 * k3y + k4y) / 6.0;

    vx[i + 1] = vx[i] + (k1vx + 2 * k2vx + 2 * k3vx + k4vx) / 6.0;
    vy[i + 1] = vy[i] + (k1vy + 2 * k2vy + 2 * k3vy + k4vy) / 6.0;
  }

  std::ofstream outfile;
  outfile.open("data3_rk4.dat");
  for (int i = 0; i < num_points; i++) {
    outfile << x[i] << "\t" << y[i] << std::endl;
  }
  outfile.close();

}

double R3(double x, double y) {
  return pow(hypot(x,y),3);
}

double dvx(double G, double m, double x, double y) {
  return -G*m/R3(x,y)*(x);
}

double dvy(double G, double m, double x, double y) {
  return -G*m/R3(x,y)*(y);
}
