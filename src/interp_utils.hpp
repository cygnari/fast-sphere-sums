#ifndef H_INTERP_UTILS_H
#define H_INTERP_UTILS_H

#include "structs.hpp"

void fekete_init(std::vector<std::vector<double>> &points, const int degree);

void interp_mat_init_sbb(std::vector<double> &mat, const std::vector<std::vector<double>> &points,
   const int degree, const int point_count);

std::vector<double> interp_vals_sbb(const double s, const double t, const double u, const int degree);

double interp_eval_sbb(const std::vector<double> &alphas, const double s, const double t, const double u,
   const int degree);

double bli_coeff(const int j, const int degree);

std::vector<std::vector<double>> interp_vals_bli(const double xi, const double eta, const double min_xi, const double max_xi, const double min_eta, const double max_eta, const int degree);

std::vector<double> bli_interp_points_shift(const double min_x, const double max_x, const int degree);

#endif
