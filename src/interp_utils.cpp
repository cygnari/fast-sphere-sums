#include <vector>
#include <cmath>
#include <tuple>
#include "general_utils.hpp"

void fekete_init(std::vector<std::vector<double>> &points, const int degree) {
  double delta_x = 1.0 / degree;
  int index;
  double a, b, c, part;
  for (int i = 0; i < degree + 1; i++) {
    for (int j = 0; j < i + 1; j++) {
      index = i * (i + 1) / 2 + j;
      a = 1 - i * delta_x;
      c = j * delta_x;
      b = 1 - a - c;
      b = 0.5 * (1 + sin(M_PI / 2 * (2 * b - 1)));
      c = 0.5 * (1 + sin(M_PI / 2 * (2 * c - 1)));
      a = 0.5 * (1 + sin(M_PI / 2 * (2 * a - 1)));
      part = a + b + c;
      points[index][0] = a / part;
      points[index][1] = b / part;
      points[index][2] = c / part;
    }
  }
}

void interp_mat_init_sbb(
    std::vector<double> &mat, const std::vector<std::vector<double>> &points,
    const int degree,
    const int point_count) { // sets up matrix to interpolate with fekete points
  // uses spherical bezier bernstein polynomials
  // for example, for deg 2, evaluates s^2, t^2, u^2, st, su, tu at interpolation points
  int index = 0, place;
  double s, t, u, si, tj, uk, comp, val, spart = 1, tpart=1;
  for (int k = 0; k < point_count; k++) {
    s = points[k][0];
    t = points[k][1];
    u = points[k][2];
    index = 0;
    spart = 1;
    for (int i = 0; i < degree + 1; i++) {
      tpart = 1;
      for (int j = 0; j < degree+1-i; j++) {
        place = point_count * index + k;
        mat[place] = spart * tpart * pow(u, degree-i-j);
        index++;
        tpart *= t;
      }
      spart *= s;
    }
  }
}

std::vector<double> interp_vals_sbb(const double s, const double t, const double u, const int degree) {
  // returns vector of SBB basis values of s, t, u
  int count = (degree + 1) * (degree + 2) / 2;
  std::vector<double> out_vals (count, 0);
  double val, factor, spart, tpart, upart;
  factor = t / u;
  int index = 0;
  spart = 1;
  for (int i = 0; i < degree + 1; i++) { // degree of s
    tpart = 1;
    for (int j = 0; j < degree + 1 - i; j++) {
      upart = pow(u, degree-i-j);
      out_vals[index] = spart*tpart*upart;
      index += 1;
      tpart *= t;
    }
    spart *= s;
  }
  return out_vals;
}

double interp_eval_sbb(
    const std::vector<double> &alphas, const double s, const double t, const double u,
    const int degree) { // evaluate SBB interpolation polynomial with coefficients
                        // alpha and barycentric point (s, t, u)
  double accum = 0;
  int index = 0;
  double val, factor, spart = 1, tpart, upart;
  for (int i = 0; i < degree + 1; i++) { // degree of s
    tpart = 1;
    for (int j = 0; j < degree + 1 - i; j++) {
      upart = pow(u, degree-i-j);
      val = spart * tpart * upart;
      accum += val * alphas[index];
      index += 1;
      tpart *= t;
    }
    spart *= s;
  }
  return accum;
}

double bli_coeff(const int j, const int degree) {
  // bli weight
  if (degree == 0) {
    return 1;
  } else if (j == 0) {
    return 0.5;
  } else if (j == degree) {
    return 0.5 * pow(-1, j);
  } else {
    return pow(-1, j);
  }
}

std::vector<std::vector<double>> interp_vals_bli(const double xi, const double eta, const double min_xi, const double max_xi, const double min_eta, const double max_eta, const int degree) {
  // returns matrix of BLI basis values
  std::vector<double> cheb_points (degree+1, 0), bli_weights (degree+1, 0);
  std::vector<std::vector<double>> func_vals (degree+1, std::vector<double> (degree+1, 0));
  if (degree == 0) {
    cheb_points[0] = 0;
    bli_weights[0] = 1;
  } else {
    for (int i = 0; i < degree+1; i++) {
      cheb_points[i] = cos(M_PI / degree * i);
      bli_weights[i] = bli_coeff(i, degree);
    }
  }
  // translate weights
  std::vector<double> cheb_points_xi (degree+1, 0);
  std::vector<double> cheb_points_eta (degree+1, 0);
  double xi_range = 0.5*(max_xi - min_xi);
  double eta_range = 0.5*(max_eta - min_eta);
  double xi_offset = 0.5*(max_xi + min_xi);
  double eta_offset = 0.5*(max_eta + min_eta);
  for (int i = 0; i < degree+1; i++) {
    cheb_points_xi[i] = xi_range*cheb_points[i]+xi_offset;
    cheb_points_eta[i] = eta_range*cheb_points[i]+eta_offset;
  }
  // compute xi basis vals
  bool found_xi_point = false;
  double denom_xi, val;
  std::vector<double> xi_func_vals (degree+1, 0);
  for (int i = 0; i < degree+1; i++) {
    if (std::abs(xi-cheb_points_xi[i]) < 1e-16) {
      found_xi_point = true;
      xi_func_vals[i] = 1;
    }
  }
  if (not found_xi_point) {
    denom_xi = 0;
    for (int i = 0; i < degree+1; i++) {
      val = bli_weights[i] / (xi - cheb_points_xi[i]);
      xi_func_vals[i] = val;
      denom_xi += val;
    }
    for (int i = 0; i < degree+1; i++) {
      xi_func_vals[i] /= denom_xi;
    }
  }
  // compute eta basis vals
  bool found_eta_point = false;
  double denom_eta;
  std::vector<double> eta_func_vals (degree+1, 0);
  for (int i = 0; i < degree+1; i++) {
    if (std::abs(eta-cheb_points_eta[i]) < 1e-16) {
      found_eta_point = true;
      eta_func_vals[i] = 1;
    }
  }
  if (not found_eta_point) {
    denom_eta = 0;
    for (int i = 0; i < degree+1; i++) {
      val = bli_weights[i] / (eta - cheb_points_eta[i]);
      eta_func_vals[i] = val;
      denom_eta += val;
    }
    for (int i = 0; i < degree+1; i++) {
      eta_func_vals[i] /= denom_eta;
    }
  }
  // compute outer product
  for (int i = 0; i < degree+1; i++) { // xi loop
    for (int j = 0; j < degree+1; j++) { // eta loop
      func_vals[i][j] = xi_func_vals[i] * eta_func_vals[j];
    }
  }
  return func_vals;
}

std::vector<double> bli_interp_points_shift(const double min_x, const double max_x, const int degree) {
  // returns xi points in this range
  std::vector<double> cheb_points (degree+1, 0);
  double x_range = 0.5*(max_x - min_x);
  double x_offset = 0.5*(max_x + min_x);
  if (degree == 0) {
    cheb_points[0] = x_offset;
  } else {
    for (int i = 0; i < degree+1; i++) {
      cheb_points[i] = cos(M_PI / degree * i)*x_range + x_offset;
    }
  }
  return cheb_points;
}

std::tuple<int, int> find_leaf_tri(double tx, double ty, double tz, const std::vector<double> &xcos, const std::vector<double> &ycos, const std::vector<double> &zcos,
                     const std::vector<std::vector<std::vector<int>>> &dynamics_triangles, const std::vector<std::vector<bool>> &dynamics_triangles_is_leaf, const int max_level) {
  bool found_leaf_tri = false, found_curr_level;
  double epsilon = pow(10, -10);
  int curr_level = -1;
  int lb = 0;
  int ub = 20;
  int iv1, iv2, iv3, tri_loc = -1;
  std::vector<double> v1, v2, v3, bary_cords;
  std::vector<double> tri_cent;
  double curr_best_dist = INT_MAX, dist, found_tri_radius;
  int curr_best_tri;

  for (int level = 0; level < max_level; level++) {
    found_curr_level = false;

    for (int j = lb; j < ub; j++) {
      iv1 = dynamics_triangles[level][j][0];
      iv2 = dynamics_triangles[level][j][1];
      iv3 = dynamics_triangles[level][j][2];
      v1 = {xcos[iv1], ycos[iv1], zcos[iv1]};
      v2 = {xcos[iv2], ycos[iv2], zcos[iv2]};
      v3 = {xcos[iv3], ycos[iv3], zcos[iv3]};
      // bary_cords = barycoords(v1, v2, v3, tx, ty, tz);

      if (check_in_tri_thresh(v1, v2, v3, tx, ty, tz, epsilon)) {
        found_curr_level = true;
        curr_level = level;
        tri_loc = j;
        if (dynamics_triangles_is_leaf[level][j]) {
          found_leaf_tri = true;
          tri_loc = j;
          curr_level = level;
          break;
        } else {
          lb = 4 * j;
          ub = 4 * j + 4;
          break;
        }
      }
    }
    if (found_leaf_tri)
      break;
    if (not found_curr_level) {
      for (int j = 0; j < 20 * pow(4, level); j++) {
        iv1 = dynamics_triangles[level][j][0];
        iv2 = dynamics_triangles[level][j][1];
        iv3 = dynamics_triangles[level][j][2];
        v1 = {xcos[iv1], ycos[iv1], zcos[iv1]};
        v2 = {xcos[iv2], ycos[iv2], zcos[iv2]};
        v3 = {xcos[iv3], ycos[iv3], zcos[iv3]};
        // bary_cords = barycoords(v1, v2, v3, tx, ty, tz);
        if (check_in_tri_thresh(v1, v2, v3, tx, ty, tz, epsilon)) {
          found_curr_level = true;
          curr_level = level;
          tri_loc = j;
          if (dynamics_triangles_is_leaf[level][j]) {
            found_leaf_tri = true;
            tri_loc = j;
            curr_level = level;
            break;
          } else {
            lb = 4 * j;
            ub = 4 * j + 4;
            break;
          }
        }
      }
    }
    if (found_leaf_tri)
      break;
  }

  if (found_leaf_tri) {
    return std::make_tuple(curr_level, tri_loc);
  } else {
    for (int i = max_level - 1; i > 0; i--) {
      for (int j = 0; j < 20 * pow(4, i); j++) {

        iv1 = dynamics_triangles[i][j][0];
        iv2 = dynamics_triangles[i][j][1];
        iv3 = dynamics_triangles[i][j][2];
        v1 = {xcos[iv1], ycos[iv1], zcos[iv1]};
        v2 = {xcos[iv2], ycos[iv2], zcos[iv2]};
        v3 = {xcos[iv3], ycos[iv3], zcos[iv3]};
        // bary_cords = barycoords(v1, v2, v3, target_point);
        if (check_in_tri_thresh(v1, v2, v3, tx, ty, tz, epsilon)) {
          curr_level = i;
          tri_loc = j;
          return std::make_tuple(curr_level, tri_loc);
        }
      }
    }
    throw std::runtime_error("Failed to locate point in a leaf triangle");
    return std::make_tuple(-1, -1);
  }
}
