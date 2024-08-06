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

void sh_10(const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos, std::vector<double>& potential) {
  double constant = sqrt(3.0/(4.0*M_PI));
  for (int i = 0; i < xcos.size(); i++) {
    potential[i] = constant * zcos[i];
  }
}

void one(const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos, std::vector<double>& potential) {
  for (int i = 0; i < potential.size(); i++) {
    potential[i] = 1.0;
  }
}

void z_co(const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos, std::vector<double>& potential) {
  for (int i = 0; i < potential.size(); i++) {
    potential[i] = zcos[i];
  }
}

void test_func_one(const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos, std::vector<double>& potential) {
  double x, y, z;
  for (int i = 0; i < potential.size(); i++) {
    x = xcos[i], y = ycos[i], z = zcos[i];
    potential[i] = 1.0 / (10+2.0*x+3.0*y*y + 4.0*z*z*z+0.5*x*y+0.25*y*z+0.33*x*z+x*y*z);
  }
}

void x_cubed(const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos, std::vector<double>& potential) {
  for (int i = 0; i < potential.size(); i++) {
    potential[i] = pow(xcos[i], 3);
  }
}

void initialize_condition(const RunConfig& run_information, const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos, std::vector<double>& potential) {
  // initial condition
  if (run_information.initial_condition == "SH43") {
    // 4 3 spherical harmonic
    sh_43(xcos, ycos, zcos, potential);
  } else if (run_information.initial_condition == "SH10") {
    sh_10(xcos, ycos, zcos, potential);
  } else if (run_information.initial_condition == "One") {
    one(xcos, ycos, zcos, potential);
  } else if (run_information.initial_condition == "Z") {
    z_co(xcos, ycos, zcos, potential);
  } else if (run_information.initial_condition == "tf1") {
    test_func_one(xcos, ycos, zcos, potential);
  } else if (run_information.initial_condition == "x3") {
    x_cubed(xcos, ycos, zcos, potential);
  }
}

void balance_conditions(std::vector<double>& potential, const std::vector<double>& area) {
  // insure that integral of potential is 0
  double int_pot = 0;
  for (int i = 0; i < potential.size(); i++) {
    int_pot += potential[i] * area[i];
  }
  if (std::abs(int_pot) > 1e-15) {
    for (int i = 0; i < potential.size(); i++) {
      potential[i] -= int_pot / (4 * M_PI);
    }
  }
}
