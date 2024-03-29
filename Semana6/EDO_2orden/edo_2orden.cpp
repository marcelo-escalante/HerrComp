#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <math.h>

std::vector<double> euler(double y0, double v0, std::vector<double> &t, double k, double m, double b);
std::vector<double> leap_frog(double y0, double v0, std::vector<double> &t, double k, double m);
std::vector<double> rk4(double y0, double v0, std::vector<double> &t, double k, double m, double b);
void export_data(std::vector<double> &t, std::vector<double> &y, std::string s);

int main(void) {

  // y'' = -k/m*y - b*y'

  double y0 = 0.1;
  double v0 = 0.0;

  double k = 50.0;
  double m = 0.2;
  double b = 0.08;

  double start = 0.0;
  double end = 3.0;
  double const h = 0.0001;
  double num_points = num_points = (end - start) / h + 1;

  std::vector<double> t(num_points);

  for (int i = 0; i < num_points; ++i) {
    t[i] = start + i * h;
  }

  // Sin fricción

  std::vector<double> y1 = euler(y0, v0, t, k, m, 0.0);
  export_data(t, y1, "data1_euler.dat");

  std::vector<double> y2 = leap_frog(y0, v0, t, k, m);
  export_data(t, y2, "data2_leapfrog.dat");

  // Solución analítica sin fricción
  std::vector<double> ysol(num_points);
  double w = sqrt(k/m);
  for (int i = 0; i < num_points; ++i){
    ysol[i] = v0/w*sin(w*t[i]) + y0*cos(w*t[i]);
  }

  export_data(t, ysol, "data3_analitica.dat");

  std::vector<double> error_euler(num_points);
  std::vector<double> error_leapfrog(num_points);

  for (int i = 0; i < num_points; ++i){
    error_euler[i] = std::abs(ysol[i]-y1[i]);
    error_leapfrog[i] = std::abs(ysol[i]-y2[i]);
  }

  export_data(t, error_euler, "data4_error_euler.dat");
  export_data(t, error_leapfrog, "data5_error_leapfrog.dat");

  // Con fricción

  std::vector<double> y3 = euler(y0, v0, t, k, m, b);
  export_data(t, y3, "data6_euler_friccion.dat");

  std::vector<double> y4 = rk4(y0, v0, t, k, m, b);
  export_data(t, y4, "data7_rk4_friccion.dat");

  // Solución analítica con fricción
  // Caso subamortiguado
  double zeta = b/2;
  double wd = sqrt(k/m-pow(zeta,2));
  double c1 = y0, c2 = (v0+zeta*y0)/wd;
  std::vector<double> ysol2(num_points);
  for (int i = 0; i < num_points; ++i){
    ysol2[i] = exp(-zeta*t[i])*(c1*cos(wd*t[i]) + c2*sin(wd*t[i]));
  }

  export_data(t, ysol2, "data8_analitica_friccion.dat");

  std::vector<double> error_euler_friccion(num_points);
  std::vector<double> error_rk4_friccion(num_points);

  for (int i = 0; i < num_points; ++i){
    error_euler_friccion[i] = std::abs(ysol2[i]-y3[i]);
    error_rk4_friccion[i] = std::abs(ysol2[i]-y4[i]);
  }

  export_data(t, error_euler_friccion, "data9_error_euler_friccion.dat");
  export_data(t, error_rk4_friccion, "data10_error_rk4_friccion.dat");

  return 0;

}

std::vector<double> euler(double y0, double v0, std::vector<double> &t, double k, double m, double b) {
  const int num_points = t.size();
  std::vector<double> y(num_points);
  std::vector<double> v(num_points);
  double h = (t.back() - t.front()) / (num_points - 1);
  y[0] = y0;
  v[0] = v0;
  for (int i = 0; i < num_points - 1; i++) {
    y[i + 1] = y[i] + h * v[i];
    v[i + 1] = v[i] + h * (-k / m * y[i] - b * v[i]);
  }

  return y;

}

std::vector<double> leap_frog(double y0, double v0, std::vector<double> &t, double k, double m) {
  const int num_points = t.size();
  std::vector<double> y(num_points);
  std::vector<double> v(num_points);
  double h = (t.back() - t.front()) / (num_points - 1);
  y[0] = y0;
  v[0] = v0 + h * 0.5 * (-k / m * y[0]);
  for (int i = 0; i < num_points - 1; i++) {
    y[i + 1] = y[i] + h * v[i];
    v[i + 1] = v[i] + h * (-k / m * y[i+1]);
  }

  return y;

}

std::vector<double> rk4(double y0, double v0, std::vector<double> &t, double k, double m, double b) {
  const int num_points = t.size();
  std::vector<double> y(num_points);
  std::vector<double> v(num_points);
  double h = (t.back() - t.front()) / (num_points - 1);
  y[0] = y0;
  v[0] = v0;

  for (int i = 0; i < num_points - 1; i++) {

    double k1y = h * v[i];
    double k1v = h * (-k / m * y[i] - b * v[i]);

    double k2y = h * (v[i] + k1v / 2);
    double k2v = h * (-k / m * (y[i] + k1y / 2) - b * (v[i] + k1v / 2));

    double k3y = h * (v[i] + k2v / 2);
    double k3v = h * (-k / m * (y[i] + k2y / 2) - b * (v[i] + k2v / 2));

    double k4y = h * (v[i] + k3v);
    double k4v = h * (-k / m * (y[i] + k3y) - b * (v[i] + k3v));

    y[i + 1] = y[i] + (k1y + 2 * k2y + 2 * k3y + k4y) / 6.0;
    v[i + 1] = v[i] + (k1v + 2 * k2v + 2 * k3v + k4v) / 6.0;
  }

  return y;

}

void export_data(std::vector<double> &t, std::vector<double> &y, std::string s){
  const int num_points = t.size();
  std::ofstream outfile;
  outfile.open(s);
  for (int i = 0; i < num_points; i++) {
    outfile << t[i] << "\t" << y[i] << std::endl;
  }
  outfile.close();

}