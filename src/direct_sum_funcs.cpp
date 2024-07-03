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

void direct_sum_invert_laplacian_reg(const RunConfig& run_information, const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos, const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral) {
  // perform direct summation to convolve to invert the laplacian
  double tx, ty, tz, sx, sy, sz;
  double eps=run_information.kernel_eps;
  // for (int i = 0; i < xcos.size(); i++) { // loop over targets
  for (int i = run_information.two_d_one_lb; i < run_information.two_d_one_ub; i++) {
    tx = xcos[i], ty = ycos[i], tz = zcos[i];
    // for (int j = 0; j < xcos.size(); j++) { // loop over sources
    for (int j = run_information.two_d_two_lb; j < run_information.two_d_two_ub; j++) {
      // if (i != j) {
      // skip singularity
      sx = xcos[j], sy = ycos[j], sz = zcos[j];
      integral[i] += -1.0/(4.0*M_PI) * log(1-tx*sx-ty*sy-tz*sz+eps)*potential[j]*area[j];
      // }
    }
  }
}

void direct_sum_laplacian(const RunConfig& run_information, const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos, const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral) {
  // perform direct summation to convolve to invert the laplacian
  double tx, ty, tz, sx, sy, sz;
  double eps=run_information.kernel_eps;
  double h1, h2, gamm, gamm2, gamm4;
  for (int i = run_information.two_d_one_lb; i < run_information.two_d_one_ub; i++) {
    tx = xcos[i], ty = ycos[i], tz = zcos[i];
    for (int j = run_information.two_d_two_lb; j < run_information.two_d_two_ub; j++) {
      sx = xcos[j], sy = ycos[j], sz = zcos[j];
      gamm = tx * sx + ty * sy + tz * sz;
      gamm2 = gamm * gamm;
      gamm4 = gamm2 * gamm2;
      h1 = exp(-1/eps)*pow(M_PI,-1.5)*(2*gamm*(-1.0+3*eps+gamm2));
      h2 = exp((gamm2-1)/eps)/(sqrt(eps)*M_PI)*(2*eps*eps+2*gamm2*(gamm2-1)+eps*(7*gamm2-1))*(1+std::erf(gamm/eps));
      integral[i] += (h1+h2)/pow(eps,2.5) * (potential[i]-potential[j]) * area[j];
    }
  }
}

void direct_sum_inverse_biharmonic(const RunConfig& run_information, const std::vector<double>& xcos_t, const std::vector<double>& ycos_t, const std::vector<double>& zcos_t,
                                  const std::vector<double>& xcos_s, const std::vector<double>& ycos_s, const std::vector<double>& zcos_s,
                                  const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral) {
  // perform direct summation to convolve to invert the bilaplacian
  double tx, ty, tz, sx, sy, sz;
  const double gfc = 2.0 - pow(M_PI, 2)/6.0;
  for (int i = run_information.two_d_one_lb; i < run_information.two_d_one_ub; i++) {
    tx = xcos_t[i], ty = ycos_t[i], tz = zcos_t[i];
    for (int j = run_information.two_d_two_lb; j < run_information.two_d_two_ub; j++) {
      sx = xcos_s[j], sy = ycos_s[j], sz = zcos_s[j];
      integral[i] += -1.0/(4.0*M_PI)*dilog(0.5*(1-tx*sx-ty*sy-tz*sz))*potential[j]*area[j];
    }
  }
}

void direct_sum_sal(const RunConfig& run_information, const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos, const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral) {
  // perform direct summation to compute the SAL potential
  double tx, ty, tz, sx, sy, sz;
  // for (int i = 0; i < xcos.size(); i++) { // loop over targets
  for (int i = run_information.two_d_one_lb; i < run_information.two_d_one_ub; i++) {
    tx = xcos[i], ty = ycos[i], tz = zcos[i];
    // for (int j = 0; j < xcos.size(); j++) { // loop over sources
    for (int j = run_information.two_d_two_lb; j < run_information.two_d_two_ub; j++) {
      sx = xcos[j], sy = ycos[j], sz = zcos[j];
      integral[i] += sal_gf_interp_40(tx*sx+ty*sy+tz*sz)*potential[j]*area[j];
    }
  }
}

void direct_sum_sal_lat_deriv(const RunConfig& run_information, const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos, const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral) {
  // perform direct summation to compute the SAL potential
  double tx, ty, tz, sx, sy, sz;
  // for (int i = 0; i < xcos.size(); i++) { // loop over targets
  for (int i = run_information.two_d_one_lb; i < run_information.two_d_one_ub; i++) {
    tx = xcos[i], ty = ycos[i], tz = zcos[i];
    // for (int j = 0; j < xcos.size(); j++) { // loop over sources
    for (int j = run_information.two_d_two_lb; j < run_information.two_d_two_ub; j++) {
      sx = xcos[j], sy = ycos[j], sz = zcos[j];
      integral[i] += sal_gf_lat_deriv(tx, ty, tz, sx, sy, sz)*potential[j]*area[j];
    }
  }
}

void direct_sum_sal_lon_deriv(const RunConfig& run_information, const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos, const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral) {
  // perform direct summation to compute the SAL potential
  double tx, ty, tz, sx, sy, sz;
  // for (int i = 0; i < xcos.size(); i++) { // loop over targets
  for (int i = run_information.two_d_one_lb; i < run_information.two_d_one_ub; i++) {
    tx = xcos[i], ty = ycos[i], tz = zcos[i];
    // for (int j = 0; j < xcos.size(); j++) { // loop over sources
    for (int j = run_information.two_d_two_lb; j < run_information.two_d_two_ub; j++) {
      sx = xcos[j], sy = ycos[j], sz = zcos[j];
      integral[i] += sal_gf_lon_deriv(tx, ty, tz, sx, sy, sz)*potential[j]*area[j];
    }
  }
}
