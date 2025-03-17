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

void direct_sum_bve(const RunConfig& run_information, const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos, const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral_1, std::vector<double>& integral_2, std::vector<double>& integral_3) {
  // perform direct summation to convolve to invert the laplacian
  double tx, ty, tz, sx, sy, sz, part;
  std::fill(integral_1.begin(), integral_1.end(), 0);
  std::fill(integral_2.begin(), integral_2.end(), 0);
  std::fill(integral_3.begin(), integral_3.end(), 0);
  // for (int i = 0; i < xcos.size(); i++) { // loop over targets
  for (int i = run_information.two_d_one_lb; i < run_information.two_d_one_ub; i++) {
    tx = xcos[i], ty = ycos[i], tz = zcos[i];
    // for (int j = 0; j < xcos.size(); j++) { // loop over sources
    for (int j = run_information.two_d_two_lb; j < run_information.two_d_two_ub; j++) {
      if (i != j) {
      // skip singularity
        sx = xcos[j], sy = ycos[j], sz = zcos[j];
        part = -1.0/(4.0*M_PI*(1.0-tx*sx-ty*sy-tz*sz))*potential[j]*area[j];
        integral_1[i] += (ty*sz-tz*sy)*part;
        integral_2[i] += (tz*sx-tx*sz)*part;
        integral_3[i] += (tx*sy-ty*sx)*part;
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

void direct_sum_dirac_delta_fish(const RunConfig& run_information, const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos, const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral) {
  // perform direct summation to convolve to invert the laplacian
  double tx, ty, tz, sx, sy, sz;
  double eps=run_information.kernel_eps;
  double gamma, denom, part;
  // for (int i = 0; i < xcos.size(); i++) { // loop over targets
  for (int i = run_information.two_d_one_lb; i < run_information.two_d_one_ub; i++) {
    tx = xcos[i], ty = ycos[i], tz = zcos[i];
    // for (int j = 0; j < xcos.size(); j++) { // loop over sources
    for (int j = run_information.two_d_two_lb; j < run_information.two_d_two_ub; j++) {
      sx = xcos[j], sy = ycos[j], sz = zcos[j];
      gamma = tx*sx+ty*sy+tz*sz;
      denom = 1+eps-gamma;
      part = 1-gamma;
      // integral[i] += ((2*eps+part*part)/(4*M_PI*denom*denom)+part/(M_PI*denom))*potential[j]*area[j]; // cartesian derivatives
      integral[i] += -(part*part-2*eps*gamma)/(4*M_PI*denom*denom)*potential[j]*area[j];
    }
  }
}

void direct_sum_laplacian_fish(const RunConfig& run_information, const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos, const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral) {
  // perform direct summation to convolve to invert the laplacian
  double tx, ty, tz, sx, sy, sz;
  double eps=run_information.kernel_eps;
  double gamma, denom, part;
  // for (int i = 0; i < xcos.size(); i++) { // loop over targets
  for (int i = run_information.two_d_one_lb; i < run_information.two_d_one_ub; i++) {
    tx = xcos[i], ty = ycos[i], tz = zcos[i];
    // for (int j = 0; j < xcos.size(); j++) { // loop over sources
    for (int j = run_information.two_d_two_lb; j < run_information.two_d_two_ub; j++) {
      sx = xcos[j], sy = ycos[j], sz = zcos[j];
      gamma = tx*sx+ty*sy+tz*sz;
      denom = 1+eps-gamma;
      part = 1+gamma;
      integral[i] += -eps*(2+eps)*(-4-2*eps*gamma+part*part)/(2*M_PI*denom*denom*denom*denom)*potential[j]*area[j];
    }
  }
}

void direct_sum_gradient_x_pse(const RunConfig& run_information, const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos, const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral_1, std::vector<double>& integral_2, std::vector<double>& integral_3) {
  // perform direct summation to convolve to invert the laplacian
  double tx, ty, tz, sx, sy, sz;
  double eps=run_information.kernel_eps;
  double h1, h2, gamm, gamm2, gamm4, gamm6, g2m1;
  double eps2, epss, eps3, eps4;
  double piece1, piece2, coeff1, coeff2, coeff3, coeff4;
  for (int i = run_information.two_d_one_lb; i < run_information.two_d_one_ub; i++) {
    tx = xcos[i], ty = ycos[i], tz = zcos[i];
    for (int j = run_information.two_d_two_lb; j < run_information.two_d_two_ub; j++) {
      sx = xcos[j], sy = ycos[j], sz = zcos[j];
      gamm = tx * sx + ty * sy + tz * sz;
      if (gamm > 1.0-30*eps) {
        gamm2 = gamm*gamm;
        gamm4=gamm2*gamm2;
        eps2=eps*eps;
        // fourth order kernel
        // piece1 = (2-6*eps-2*gamm2)*gamm;
        // piece2 = 3*eps2-2*eps-2*gamm2+8*eps*gamm2+2*gamm4;
        // coeff1 = 2*exp(1/eps)*pow(M_PI,1.5);
        // coeff2 = exp((gamm2-1)/eps)/(2*sqrt(eps)*M_PI);
        // coeff3 = 1+std::erf(gamm/sqrt(eps));
        // // coeff4 = 1.0/pow(eps,2.5)*(potential[i]-potential[j])*area[j];
        // coeff4 = -1.0/pow(eps,2.5)*potential[j]*area[j];
        // h1 = (tx*piece1+sx*piece2)/coeff1;
        // h2 = coeff2*(2*eps2*(tx-3*gamm*sx)+2*gamm2*(gamm2-1)*(tx-gamm*sx)+(7*gamm2-1)*tx*eps+3*gamm*sx*eps*(1-3*gamm2))*coeff3;
        // integral_1[i] += (h1+h2)*coeff4;
        // h1 = (ty*piece1+sy*piece2)/coeff1;
        // h2 = coeff2*(2*eps2*(ty-3*gamm*sy)+2*gamm2*(gamm2-1)*(ty-gamm*sy)+(7*gamm2-1)*ty*eps+3*gamm*sy*eps*(1-3*gamm2))*coeff3;
        // integral_2[i] += (h1+h2)*coeff4;
        // h1 = (tz*piece1+sz*piece2)/coeff1;
        // h2 = coeff2*(2*eps2*(tz-3*gamm*sz)+2*gamm2*(gamm2-1)*(tz-gamm*sz)+(7*gamm2-1)*tz*eps+3*gamm*sz*eps*(1-3*gamm2))*coeff3;
        // integral_3[i] += (h1+h2)*coeff4;
        // second order kernel
        piece2 = eps*(eps+gamm2);
        piece1 = -eps*gamm;
        coeff1 = exp(1/eps)*pow(M_PI,1.5);
        coeff2 = -exp((gamm2-1)/eps)*sqrt(eps)/(2*M_PI);
        coeff3 = 1+std::erf(gamm/sqrt(eps));
        coeff4 = 1.0/pow(eps,2.5)*potential[j]*area[j];
        h1 = (tx*piece1+sx*piece2)/coeff1;
        h2 = coeff2*coeff3*(eps*(tx-3*gamm*sx)+2*gamm2*(tx-gamm*sx));
        integral_1[i] += (h1+h2)*coeff4;
        h1 = (ty*piece1+sy*piece2)/coeff1;
        h2 = coeff2*coeff3*(eps*(ty-3*gamm*sy)+2*gamm2*(ty-gamm*sy));
        integral_2[i] += (h1+h2)*coeff4;
        h1 = (tz*piece1+sz*piece2)/coeff1;
        h2 = coeff2*coeff3*(eps*(tz-3*gamm*sz)+2*gamm2*(tz-gamm*sz));
        integral_3[i] += (h1+h2)*coeff4;
      }    
    }
  }
}

void direct_sum_laplacian_pse(const RunConfig& run_information, const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos, const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral) {
  // perform direct summation to convolve to invert the laplacian
  double tx, ty, tz, sx, sy, sz;
  double eps=run_information.kernel_eps;
  double h1, h2, gamm, gamm2, gamm4, gamm6, g2m1;
  double eps2, epss, eps3, eps4;
  for (int i = run_information.two_d_one_lb; i < run_information.two_d_one_ub; i++) {
    tx = xcos[i], ty = ycos[i], tz = zcos[i];
    for (int j = run_information.two_d_two_lb; j < run_information.two_d_two_ub; j++) {
      sx = xcos[j], sy = ycos[j], sz = zcos[j];
      gamm = tx * sx + ty * sy + tz * sz;
      // 4/eps in H_{4,eps} integral
      if (gamm > 1.0 - 30*eps) {
        gamm2 = gamm * gamm;
        gamm4 = gamm2 * gamm2;
        h1 = exp(-1/eps)*pow(M_PI,-1.5)*(2*gamm*(-1.0+3*eps+gamm2));
        h2 = exp((gamm2-1)/eps)/(sqrt(eps)*M_PI)*(2*eps*eps+2*gamm2*(gamm2-1)+eps*(7*gamm2-1))*(1+std::erf(gamm/sqrt(eps)));
        integral[i] += (h1+h2)/pow(eps,2.5) * (potential[i]-potential[j]) * area[j];
      }
      // H_{2,eps}
      // if (gamm > 1.0 - 30*eps) {
      //   gamm2 = gamm * gamm;
      //   eps2 = sqrt(eps);
      //   h1 = exp(-1/eps)*pow(M_PI,-1.5)*(2*gamm*eps);
      //   h2 = exp((gamm2-1)/eps)/M_PI*eps2*(eps+2*gamm2)*(1+std::erf(gamm/eps2));
      //   integral[i] += (h1+h2)/pow(eps,2.5) * (potential[i]-potential[j]) * area[j];
      // }
      // H_{6,eps}
      // if (gamm > 1.0-30*eps) {
      //   gamm2 = gamm*gamm;
      //   gamm4 = gamm2*gamm2;
      //   gamm6 = gamm4*gamm2;
      //   epss = sqrt(eps);
      //   eps2 = eps*eps;
      //   h1 = exp(-1.0/eps)*gamm*(23*eps2+16*eps*(gamm2-1)+2*(gamm2-1)*(gamm2-1))/(2*eps*pow(M_PI,1.5));
      //   h2 = exp((gamm2-1)/eps)/(2*epss*eps2*M_PI)*(eps*(1-6*eps+6*eps2)+2*gamm2*(1-9*eps+15*eps2)+gamm4*(17*eps-4)+gamm6*2)*(1+std::erf(gamm/epss));
      //   integral[i] += (h1+h2)/pow(eps,1.5) * (potential[i]-potential[j]) * area[j];
      // }
      // H_{8,eps}
      // if (gamm > 1.0-30*eps) {
      //   gamm2 = gamm*gamm;
      //   g2m1 = gamm2-1;
      //   gamm4 = gamm2*gamm2;
      //   eps2 = eps*eps;
      //   eps4 = eps2*eps2;
      //   eps3 = eps*eps2;
      //   epss = sqrt(eps);
      //   h1 = exp(-1.0/eps)*gamm*(219*eps3+60*eps*g2m1*g2m1+4*g2m1*g2m1*g2m1+eps2*(236*gamm2-234))/(12*eps2*pow(M_PI,1.5));
      //   h2 = exp(g2m1/eps)/(6*eps2*epss*M_PI)*(24*eps4+2*gamm2*pow(g2m1,3)+12*eps3*(13*gamm2-3)+eps*g2m1*g2m1*(31*gamm2-1)+12*eps2*(1-12*gamm2+11*gamm4))*(1+std::erf(gamm/epss));
      //   integral[i] += (h1+h2)/eps*(potential[i]-potential[j])*area[j];
      // }
    }
  }
}

void direct_sum_inverse_biharmonic(const RunConfig& run_information, const std::vector<double>& xcos_t, const std::vector<double>& ycos_t, const std::vector<double>& zcos_t,
                                  const std::vector<double>& xcos_s, const std::vector<double>& ycos_s, const std::vector<double>& zcos_s,
                                  const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral) {
  // perform direct summation to convolve to invert the bilaplacian
  double tx, ty, tz, sx, sy, sz;
  // const double gfc = 2.0 - pow(M_PI, 2)/6.0;
  for (int i = run_information.two_d_one_lb; i < run_information.two_d_one_ub; i++) {
    tx = xcos_t[i], ty = ycos_t[i], tz = zcos_t[i];
    for (int j = run_information.two_d_two_lb; j < run_information.two_d_two_ub; j++) {
      sx = xcos_s[j], sy = ycos_s[j], sz = zcos_s[j];
      integral[i] += 1.0/(4.0*M_PI)*dilog(0.5*(1+tx*sx+ty*sy+tz*sz))*potential[j]*area[j];
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
