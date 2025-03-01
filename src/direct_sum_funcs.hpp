#ifndef H_DIRECT_SUM_FUNCS_H
#define H_DIRECT_SUM_FUNCS_H

#include <vector>
#include "structs.hpp"

void direct_sum_invert_laplacian(const RunConfig& run_information, const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos, const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral);

void direct_sum_bve(const RunConfig& run_information, const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos, const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral_1, std::vector<double>& integral_2, std::vector<double>& integral_3);

void direct_sum_invert_laplacian_reg(const RunConfig& run_information, const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos, const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral);

void direct_sum_dirac_delta_fish(const RunConfig& run_information, const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos, const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral);

void direct_sum_laplacian_fish(const RunConfig& run_information, const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos, const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral);

void direct_sum_gradient_x_pse(const RunConfig& run_information, const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos, const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral_1, std::vector<double>& integral_2, std::vector<double>& integral_3);

void direct_sum_laplacian_pse(const RunConfig& run_information, const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos, const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral);

void direct_sum_inverse_biharmonic(const RunConfig& run_information, const std::vector<double>& xcos_t, const std::vector<double>& ycos_t, const std::vector<double>& zcos_t,
                                  const std::vector<double>& xcos_s, const std::vector<double>& ycos_s, const std::vector<double>& zcos_s,
                                  const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral);

void direct_sum_sal(const RunConfig& run_information, const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos, const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral);

void direct_sum_sal_lat_deriv(const RunConfig& run_information, const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos, const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral);

void direct_sum_sal_lon_deriv(const RunConfig& run_information, const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos, const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral);

#endif
