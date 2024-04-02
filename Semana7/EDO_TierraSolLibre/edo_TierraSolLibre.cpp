#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>

void euler(double M, double G, double e, double m, std::vector<double> &t, 
  double x0t, double y0t, double vx0t, double vy0t, 
  double x0s, double y0s, double vx0s, double vy0s);

void leapfrog(double M, double G, double e, double m, std::vector<double> &t, 
  double x0t, double y0t, double vx0t, double vy0t, 
  double x0s, double y0s, double vx0s, double vy0s);

void rk4(double M, double G, double e, double m, std::vector<double> &t, 
  double x0t, double y0t, double vx0t, double vy0t, 
  double x0s, double y0s, double vx0s, double vy0s);

int main(){

  double M = 1.0;
  double G = 2.95900e-4;
  double e = 0.01671; 
  double a = 1.000003;
  double m = 3.0027e-6;

  double x0t = -0.983;
  double y0t = 0.0;
  double vx0t = 0.0;
  double vy0t = -sqrt(G*M/a*(1+e)/(1-e));

  double x0s = 0.0;
  double y0s = 0.0;
  double vx0s = 0.0;
  double vy0s = 0.0; 

  double start = 0.0;
  double end = 2200.0; // 2200 días (seis años y diez días)
  double h = 0.001;
  double num_points = num_points = (end - start) / h + 1;

  std::vector<double> t(num_points);
  for (int i = 0; i<num_points; ++i){
    t[i] = start + h*i;
  }

  euler(M,G,e,m,t,x0t,y0t,vx0t,vy0t,x0s,y0s,vx0s,vy0s);
  leapfrog(M,G,e,m,t,x0t,y0t,vx0t,vy0t,x0s,y0s,vx0s,vy0s);
  rk4(M,G,e,m,t,x0t,y0t,vx0t,vy0t,x0s,y0s,vx0s,vy0s);

  return 0;

}

//Implementación Euler

void euler(double M, double G, double e, double m, std::vector<double> &t, 
  double x0t, double y0t, double vx0t, double vy0t, 
  double x0s, double y0s, double vx0s, double vy0s){

  const int num_points = t.size();

  std::vector<double> xt(num_points);
  std::vector<double> vxt(num_points);
  std::vector<double> yt(num_points);
  std::vector<double> vyt(num_points);

  std::vector<double> xs(num_points);
  std::vector<double> vxs(num_points);
  std::vector<double> ys(num_points);
  std::vector<double> vys(num_points);

  std::vector<double> R3(num_points);

  double h = (t.back() - t.front()) / (num_points - 1);

  xt[0] = x0t;
  yt[0] = y0t;
  vxt[0] = vx0t;
  vyt[0] = vy0t;

  xs[0] = x0s;
  ys[0] = y0s;
  vxs[0] = vx0s;
  vys[0] = vy0s;

  for (int i = 0; i < num_points - 1; ++i) {

    R3[i] = pow(hypot((xt[i]-xs[i]),yt[i]-ys[i]),3);

    // Iteración Tierra
    xt[i + 1] = xt[i] + h * vxt[i];
    yt[i + 1] = yt[i] + h * vyt[i];
    vxt[i + 1] = vxt[i] + h * (G*M/R3[i])*(xs[i]-xt[i]);
    vyt[i + 1] = vyt[i] + h * (G*M/R3[i])*(ys[i]-yt[i]);

    //Iteración Sol
    xs[i + 1] = xs[i] + h * vxs[i];
    ys[i + 1] = ys[i] + h * vys[i];
    vxs[i + 1] = vxs[i] + h * (G*m/R3[i])*(xt[i]-xs[i]);
    vys[i + 1] = vys[i] + h * (G*m/R3[i])*(yt[i]-ys[i]);
  }

  std::ofstream outfile;
  outfile.open("data1_euler.dat");
  for (int i = 0; i < num_points; i++) {
    outfile << xt[i] << "\t" << yt[i] << "\t" << xs[i] << "\t" << ys[i] << std::endl;
  }
  outfile.close();
}

//Implementación LeapFrog

void leapfrog(double M, double G, double e, double m, std::vector<double> &t, 
  double x0t, double y0t, double vx0t, double vy0t, 
  double x0s, double y0s, double vx0s, double vy0s){

  const int num_points = t.size();

  std::vector<double> xt(num_points);
  std::vector<double> vxt(num_points);
  std::vector<double> yt(num_points);
  std::vector<double> vyt(num_points);

  std::vector<double> xs(num_points);
  std::vector<double> vxs(num_points);
  std::vector<double> ys(num_points);
  std::vector<double> vys(num_points);

  std::vector<double> R3(num_points);

  double h = (t.back() - t.front()) / (num_points - 1);

  xt[0] = x0t;
  yt[0] = y0t;
  vxt[0] = vx0t;
  vyt[0] = vy0t;

  xs[0] = x0s;
  ys[0] = y0s;
  vxs[0] = vx0s;
  vys[0] = vy0s;

  for (int i = 0; i < num_points - 1; i++) {

    xt[i + 1] = xt[i] + h * vxt[i];
    yt[i + 1] = yt[i] + h * vyt[i];

    xs[i + 1] = xs[i] + h * vxs[i];
    ys[i + 1] = ys[i] + h * vys[i];

    R3[i+1] = pow(hypot((xt[i+1]-xs[i+1]),yt[i+1]-ys[i+1]),3);

    vxt[i + 1] = vxt[i] + h * (G*M/R3[i+1])*(xs[i+1]-xt[i+1]);
    vyt[i + 1] = vyt[i] + h * (G*M/R3[i+1])*(ys[i+1]-yt[i+1]);

    vxs[i + 1] = vxs[i] + h * (G*m/R3[i+1])*(xt[i+1]-xs[i+1]);
    vys[i + 1] = vys[i] + h * (G*m/R3[i+1])*(yt[i+1]-ys[i+1]);
  }

  std::ofstream outfile;
  outfile.open("data2_leapfrog.dat");
  for (int i = 0; i < num_points; i++) {
    outfile << xt[i] << "\t" << yt[i] << "\t" << xs[i] << "\t" << ys[i] << std::endl;
  }
  outfile.close();
}

//Implementación Runge Kutta orden 4

//Funciones auxiliares
double R3(double xt, double yt, double xs, double ys);
double dvxt(double G, double m, double xt, double yt, double xs, double ys);
double dvyt(double G, double m, double xt, double yt, double xs, double ys);
double dvxs(double G, double m, double xt, double yt, double xs, double ys);
double dvys(double G, double m, double xt, double yt, double xs, double ys);
//

void rk4(double M, double G, double e, double m, std::vector<double> &t, 
  double x0t, double y0t, double vx0t, double vy0t, 
  double x0s, double y0s, double vx0s, double vy0s){

  const int num_points = t.size();

  std::vector<double> xt(num_points);
  std::vector<double> vxt(num_points);
  std::vector<double> yt(num_points);
  std::vector<double> vyt(num_points);

  std::vector<double> xs(num_points);
  std::vector<double> vxs(num_points);
  std::vector<double> ys(num_points);
  std::vector<double> vys(num_points);

  double h = (t.back() - t.front()) / (num_points - 1);

  xt[0] = x0t;
  yt[0] = y0t;
  vxt[0] = vx0t;
  vyt[0] = vy0t;

  xs[0] = x0s;
  ys[0] = y0s;
  vxs[0] = vx0s;
  vys[0] = vy0s;

  for (int i = 0; i < num_points - 1; ++i) {

    double k1xt = h * vxt[i];
    double k1yt = h * vyt[i];
    double k1xs = h * vxs[i];
    double k1ys = h * vys[i];

    double k1vxt = h * dvxt(G,M,xt[i],yt[i],xs[i],ys[i]);
    double k1vyt = h * dvyt(G,M,xt[i],yt[i],xs[i],ys[i]);
    double k1vxs = h * dvxs(G,m,xt[i],yt[i],xs[i],ys[i]);
    double k1vys = h * dvys(G,m,xt[i],yt[i],xs[i],ys[i]);

    double k2xt = h * (vxt[i] + k1vxt / 2);
    double k2yt = h * (vyt[i] + k1vyt / 2);
    double k2xs = h * (vxs[i] + k1vxs / 2);
    double k2ys = h * (vys[i] + k1vys / 2);

    double k2vxt = h * dvxt(G,M,xt[i]+k1vxt/2,yt[i]+k1vyt/2,xs[i]+k1vxs/2,ys[i]+k1vys/2);
    double k2vyt = h * dvyt(G,M,xt[i]+k1vxt/2,yt[i]+k1vyt/2,xs[i]+k1vxs/2,ys[i]+k1vys/2);
    double k2vxs = h * dvxs(G,m,xt[i]+k1vxt/2,yt[i]+k1vyt/2,xs[i]+k1vxs/2,ys[i]+k1vys/2);
    double k2vys = h * dvys(G,m,xt[i]+k1vxt/2,yt[i]+k1vyt/2,xs[i]+k1vxs/2,ys[i]+k1vys/2);

    double k3xt = h * (vxt[i] + k2vxt / 2);
    double k3yt = h * (vyt[i] + k2vyt / 2);
    double k3xs = h * (vxs[i] + k2vxs / 2);
    double k3ys = h * (vys[i] + k2vys / 2);

    double k3vxt = h * dvxt(G,M,xt[i]+k2vxt/2,yt[i]+k2vyt/2,xs[i]+k2vxs/2,ys[i]+k2vys/2);
    double k3vyt = h * dvyt(G,M,xt[i]+k2vxt/2,yt[i]+k2vyt/2,xs[i]+k2vxs/2,ys[i]+k2vys/2);
    double k3vxs = h * dvxs(G,m,xt[i]+k2vxt/2,yt[i]+k2vyt/2,xs[i]+k2vxs/2,ys[i]+k2vys/2);
    double k3vys = h * dvys(G,m,xt[i]+k2vxt/2,yt[i]+k2vyt/2,xs[i]+k2vxs/2,ys[i]+k2vys/2);

    double k4xt = h * (vxt[i] + k3vxt / 2);
    double k4yt = h * (vyt[i] + k3vyt / 2);
    double k4xs = h * (vxs[i] + k3vxs / 2);
    double k4ys = h * (vys[i] + k3vys / 2);

    double k4vxt = h * dvxt(G,M,xt[i]+k3vxt/2,yt[i]+k3vyt/2,xs[i]+k3vxs/2,ys[i]+k3vys/2);
    double k4vyt = h * dvyt(G,M,xt[i]+k3vxt/2,yt[i]+k3vyt/2,xs[i]+k3vxs/2,ys[i]+k3vys/2);
    double k4vxs = h * dvxs(G,m,xt[i]+k3vxt/2,yt[i]+k3vyt/2,xs[i]+k3vxs/2,ys[i]+k3vys/2);
    double k4vys = h * dvys(G,m,xt[i]+k3vxt/2,yt[i]+k3vyt/2,xs[i]+k3vxs/2,ys[i]+k3vys/2);

    xt[i + 1] = xt[i] + (k1xt + 2 * k2xt + 2 * k3xt + k4xt) / 6.0;
    yt[i + 1] = yt[i] + (k1yt + 2 * k2yt + 2 * k3yt + k4yt) / 6.0;
    xs[i + 1] = xs[i] + (k1xs + 2 * k2xs + 2 * k3xs + k4xs) / 6.0;
    ys[i + 1] = ys[i] + (k1ys + 2 * k2ys + 2 * k3ys + k4ys) / 6.0;
    
    vxt[i + 1] = vxt[i] + (k1vxt + 2 * k2vxt + 2 * k3vxt + k4vxt) / 6.0;
    vyt[i + 1] = vyt[i] + (k1vyt + 2 * k2vyt + 2 * k3vyt + k4vyt) / 6.0;
    vxs[i + 1] = vxs[i] + (k1vxs + 2 * k2vxs + 2 * k3vxs + k4vxs) / 6.0;
    vys[i + 1] = vys[i] + (k1vys + 2 * k2vys + 2 * k3vys + k4vys) / 6.0;

  }

  std::ofstream outfile;
  outfile.open("data3_rk4.dat");
  for (int i = 0; i < num_points; i++) {
    outfile << xt[i] << "\t" << yt[i] << "\t" << xs[i] << "\t" << ys[i] << std::endl;
  }
  outfile.close();
}

double R3(double xt, double yt, double xs, double ys) {
  return pow(hypot(xt-xs,yt-ys),3);
}

double dvxt(double G, double m, double xt, double yt, double xs, double ys) {
  return G*m/R3(xt,yt,xs,ys)*(xs-xt);
}

double dvyt(double G, double m, double xt, double yt, double xs, double ys) {
  return G*m/R3(xt,yt,xs,ys)*(ys-yt);
}

double dvxs(double G, double m, double xt, double yt, double xs, double ys) {
  return G*m/R3(xt,yt,xs,ys)*(xt-xs);
}

double dvys(double G, double m, double xt, double yt, double xs, double ys) {
  return G*m/R3(xt,yt,xs,ys)*(yt-ys);
}
