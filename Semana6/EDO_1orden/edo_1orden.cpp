#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

std::vector<double> euler(double y0, std::vector<double> &t);
std::vector<double> rk4(double y0, std::vector<double> &t);
void export_data(std::vector<double> &t, std::vector<double> &y, std::string s);


int main() {

  // y' = -y

  double start = 0.0;
  double end = 5.0;
  const double h = 0.0001;
  const double num_points = (end - start) / h + 1;

  std::vector<double> t(num_points);

  for (int i = 0; i < num_points; i++) {
    t[i] = start + i * h;
  }
  
  std::vector<double> y1 = euler(1.0, t);
  export_data(t, y1, "data1_euler.dat");

  std::vector<double> y2 = rk4(1.0, t);
  export_data(t, y2, "data2_rk4.dat");

  // Solución analítica 
  std::vector<double> ysol(num_points);
  for (int i = 0; i < num_points; i++){
    ysol[i] = exp(-t[i]);
  }

  export_data(t, ysol, "data3_analitica.dat");

  std::vector<double> error_euler(num_points);
  std::vector<double> error_rk4(num_points);

  for (int i=0; i < num_points; i++) {
    error_euler[i] = std::abs(ysol[i]-y1[i]);
    error_rk4[i] = std::abs(ysol[i]-y2[i]);
  }

  export_data(t, error_euler, "data4_error_euler.dat");
  export_data(t, error_rk4, "data5_error_rk4.dat");

  return 0;
}

std::vector<double> euler(double y0, std::vector<double> &t) {
  const int num_points = t.size();
  std::vector<double> y(num_points);
  y[0] = y0;
  double h = (t.back() - t.front()) / (num_points - 1);
  for (int i = 0; i < num_points - 1; i++) {
    y[i + 1] = y[i] - h * (y[i]);
  }

  return y;

}

std::vector<double> rk4(double y0, std::vector<double> &t) {
  const int num_points = t.size();
  std::vector<double> y(num_points);
  y[0] = y0;
  double h = (t.back() - t.front()) / (num_points - 1);

  for (int i = 0; i < num_points - 1; i++) {
    double k1 = h * (-y[i]);
    double k2 = h * (-(y[i] + k1 / 2));
    double k3 = h * (-(y[i] + k2 / 2));
    double k4 = h * (-(y[i] + k3));
    y[i + 1] = y[i] + 1.0 / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
  }

  return y;

}

void export_data(std::vector<double> &t, std::vector<double> &y, std::string s) {
  const int num_points = t.size();
  std::ofstream outfile;
  outfile.open(s);
  for (int i = 0; i < num_points; i++) {
    outfile << t[i] << "\t" << y[i] << std::endl;
  }
  outfile.close();

}