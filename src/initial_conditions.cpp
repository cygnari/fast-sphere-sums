#include <vector>
#include <iostream>
#include <queue>
#include <cmath>

#include "general_utils.hpp"
#include "structs.hpp"

void sh_43(const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos, std::vector<double>& potential) {
  // real spherical harmonic, degree 4, order 3
  double constant = 3.0/4.0*sqrt(35.0/(2.0*M_PI));
  double x, y, z;
  for (int i = 0; i < xcos.size(); i++) {
    x = xcos[i], y = ycos[i], z = zcos[i];
    potential[i] = x*(x*x-3.0*y*y)*z;
  }
}

void initialize_condition(const RunConfig& run_information, const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos, std::vector<double>& potential) {
  // initial condition
  if (run_information.initial_condition == "SH43") {
    // 4 3 spherical harmonic
    sh_43(xcos, ycos, zcos, potential);
  }
}
