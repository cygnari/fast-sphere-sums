#include <vector>
#include <cmath>
#include "structs.hpp"
#include "general_utils.hpp"

void direct_sum_invert_laplacian(const RunConfig& run_information, const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos, const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral) {
  // perform direct summation to convolve to invert the laplacian
  double tx, ty, tz, sx, sy, sz;
  // for (int i = 0; i < xcos.size(); i++) { // loop over targets
  for (int i = run_information.two_d_one_lb; i < run_information.two_d_one_ub; i++) {
    tx = xcos[i], ty = ycos[i], tz = zcos[i];
    // for (int j = 0; j < xcos.size(); j++) { // loop over sources
    for (int j = run_information.two_d_two_lb; j < run_information.two_d_two_ub; j++) {
      if (i != j) {
      // skip singularity
        sx = xcos[j], sy = ycos[j], sz = zcos[j];
        integral[i] += -1.0/(4.0*M_PI) * log(1-tx*sx-ty*sy-tz*sz)*potential[j]*area[j];
      }
    }
  }
}

void direct_sum_inverse_biharmonic(const RunConfig& run_information, const std::vector<double>& xcos_t, const std::vector<double>& ycos_t, const std::vector<double>& zcos_t,
                                  const std::vector<double>& xcos_s, const std::vector<double>& ycos_s, const std::vector<double>& zcos_s,
                                  const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral) {
  // perform direct summation to convolve to invert the laplacian
  double tx, ty, tz, sx, sy, sz;
  const double gfc = 2.0 - pow(M_PI, 2)/6.0;
  for (int i = run_information.two_d_one_lb; i < run_information.two_d_one_ub; i++) {
    tx = xcos_t[i], ty = ycos_t[i], tz = zcos_t[i];
    for (int j = run_information.two_d_two_lb; j < run_information.two_d_two_ub; j++) {
      sx = xcos_s[j], sy = ycos_s[j], sz = zcos_s[j];
      // if ((abs(tx - sx) < 1e-15) and (abs(ty - sy) < 1e-15) and (abs(tz - sz) < 1e-15)) {
        // no singularity
      // integral[i] += (gfc+dilog(0.5*(1-tx*sx-ty*sy-tz*sz)))*potential[j]*area[j];
      integral[i] += -1.0/(4.0*M_PI)*dilog(0.5*(1-tx*sx-ty*sy-tz*sz))*potential[j]*area[j];
      // integral[i] += 1.0/(1-tx*sx-ty*sy-tz*sz)*potential[j]*area[j];
      // }
    }
  }
}
