#include <vector>
#include <cmath>

#include "interp_utils.hpp"
#include "general_utils.hpp"
#include "structs.hpp"

void pp_interaction_inverse_laplacian_icos(const RunConfig& run_information, const int index_target, const int index_source, const std::vector<IcosPanel>& icos_panels,
                    const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos,
                    const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral) {
  // compute particle particle interaction
  int target_i, source_j, source_count;
  double tx, ty, tz, sx, sy, sz, gfval;

  source_count = icos_panels[index_source].point_count;

  std::vector<double> sxs (source_count, 0), sys (source_count, 0), szs (source_count, 0), pots (source_count, 0), areas (source_count, 0);

  for (int j = 0; j < source_count; j++) {
    source_j = icos_panels[index_source].points_inside[j];
    sxs[j] = xcos[source_j];
    sys[j] = ycos[source_j];
    szs[j] = zcos[source_j];
    pots[j] = potential[source_j];
    areas[j] = area[source_j];
  }

  for (int i = 0; i < icos_panels[index_target].point_count; i++) {
    target_i = icos_panels[index_target].points_inside[i];
    tx = xcos[target_i];
    ty = ycos[target_i];
    tz = zcos[target_i];
    for (int j = 0; j < source_count; j++) {
      sx = sxs[j], sy = sys[j], sz = szs[j];
      gfval = -1.0/(4.0*M_PI)*log(1-tx*sx-ty*sy-tz*sz);
      integral[target_i] += gfval * areas[j] * pots[j];
    }
  }
}

void pc_interaction_inverse_laplacian_icos(const RunConfig& run_information, const int index_target, const int index_source, const std::vector<IcosPanel>& icos_panels,
                    const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos,
                    const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral) {
  // computes particle cluster interaction
  int point_index, info, dim = run_information.interp_point_count, count_target = icos_panels[index_target].point_count, count_source = icos_panels[index_source].point_count;
  std::vector<double> v1s, v2s, v3s, target_particle, bary_cord, source_particle, interp_matrix(dim * dim, 0), proxy_weights(dim, 0),
                      basis_vals, func_vals(dim * count_target, 0), txs (count_target, 0), tys (count_target, 0), tzs (count_target, 0);
  double pot, us, vs, ws, tx, ty, tz, cx, cy, cz, sx, sy, sz, alpha, scalar;
  std::vector<std::vector<double>> interp_points(dim, std::vector<double>(3, 0));

  fekete_init(interp_points, run_information.interp_degree);
  v1s = icos_panels[index_source].vertex_1;
  v2s = icos_panels[index_source].vertex_2;
  v3s = icos_panels[index_source].vertex_3;

  for (int i = 0; i < count_target; i++) {
    point_index = icos_panels[index_target].points_inside[i];
    txs[i] = xcos[point_index];
    tys[i] = ycos[point_index];
    tzs[i] = zcos[point_index];
  }

  for (int i = 0; i < count_source; i++) { // compute proxy weights
    point_index = icos_panels[index_source].points_inside[i];
    sx = xcos[point_index];
    sy = ycos[point_index];
    sz = zcos[point_index];
    pot = potential[point_index];
    bary_cord = barycoords(v1s, v2s, v3s, sx, sy, sz);
    basis_vals = interp_vals_sbb(bary_cord[0], bary_cord[1], bary_cord[2], run_information.interp_degree);
    for (int j = 0; j < dim; j++) {
      proxy_weights[j] += basis_vals[j] * pot * area[point_index];
    }
  }

  for (int i = 0; i < dim; i++) { // set up interpolation matrix, loop over proxy source points
    us = interp_points[i][0];
    vs = interp_points[i][1];
    ws = interp_points[i][2];
    cx = us * v1s[0] + vs * v2s[0] + ws * v3s[0];
    cy = us * v1s[1] + vs * v2s[1] + ws * v3s[1];
    cz = us * v1s[2] + vs * v2s[2] + ws * v3s[2];
    scalar = run_information.radius / sqrt(cx * cx + cy * cy + cz * cz);
    cx *= scalar;
    cy *= scalar;
    cz *= scalar;
    bary_cord = barycoords(v1s, v2s, v3s, cx, cy, cz);
    interp_points[i]=bary_cord;
    // auto [clat, clon] = xyz_to_latlon(cx, cy, cz);
    for (int j = 0; j < count_target; j++) { // loop over targets
      tx = txs[j];
      ty = tys[j];
      tz = tzs[j];
      func_vals[j * dim + i] = -1.0/(4.0*M_PI)*log(1-tx*cx-ty*cy-tz*cz);
    }
  }

  interp_mat_init_sbb(interp_matrix, interp_points, run_information.interp_degree, dim);

  info = linear_solve(interp_matrix, func_vals, dim, count_target, 3);
  if (info > 0) {
    throw std::runtime_error("Error with linear solve in pc vel computation");
  }

  for (int i = 0; i < count_target; i++) {
    point_index = icos_panels[index_target].points_inside[i];
    for (int j = 0; j < dim; j++) {
      alpha = func_vals[dim*i+j];
      integral[point_index] += alpha * proxy_weights[j];
    }
  }
}

void cp_interaction_inverse_laplacian_icos(const RunConfig& run_information, const int index_target, const int index_source, const std::vector<IcosPanel>& icos_panels,
                    const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos,
                    const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral) {
  // computes cluster particle interaction
  int iv1, iv2, iv3, point_index, dim = run_information.interp_point_count, count_target = icos_panels[index_target].point_count, count_source = icos_panels[index_source].point_count;
  std::vector<double> v1, v2, v3, bary_cord, interptargets(dim, 0), interp_matrix(dim * dim, 0), sxs (count_source, 0), sys (count_source, 0), szs (count_source, 0), areas (count_source, 0), pots (count_source, 0);
  double u, v, w, vor, cx, cy, cz, slon, slat, sx, sy, sz, tx, ty, tz, scalar;
  std::vector<std::vector<double>> interp_points(dim, std::vector<double>(3, 0));

  for (int j = 0; j < count_source; j++) {
    point_index = icos_panels[index_source].points_inside[j];
    sxs[j] = xcos[point_index];
    sys[j] = ycos[point_index];
    szs[j] = zcos[point_index];
    areas[j] = area[point_index];
    pots[j] = potential[point_index];
  }

  fekete_init(interp_points, run_information.interp_degree);

  v1 = icos_panels[index_target].vertex_1;
  v2 = icos_panels[index_target].vertex_2;
  v3 = icos_panels[index_target].vertex_3;
  for (int i = 0; i < dim; i++) {
    u = interp_points[i][0];
    v = interp_points[i][1];
    w = 1.0 - u - v;
    cx = u * v1[0] + v * v2[0] + w * v3[0];
    cy = u * v1[1] + v * v2[1] + w * v3[1];
    cz = u * v1[2] + v * v2[2] + w * v3[2];
    scalar = run_information.radius / sqrt(cx * cx + cy * cy + cz * cz);
    cx *= scalar;
    cy *= scalar;
    cz *= scalar;
    bary_cord = barycoords(v1, v2, v3, cx, cy, cz);
    interp_points[i] = bary_cord;
    for (int j = 0; j < count_source; j++) {
      point_index = icos_panels[index_source].points_inside[j];
      sx = sxs[j], sy = sys[j], sz = szs[j];
      interptargets[i] += -1.0/(4.0*M_PI)*log(1-cx*sx-cy*sy-cz*sz)*pots[j]*areas[j];
    }
  }

  interp_mat_init_sbb(interp_matrix, interp_points, run_information.interp_degree, dim);

  int info = linear_solve(interp_matrix, interptargets, dim, 1, 3);
  if (info > 0) {
    throw std::runtime_error("Error with linear solve in cp vel computation");
  }

  for (int i = 0; i < count_target; i++) {
    point_index = icos_panels[index_target].points_inside[i];
    tx = xcos[point_index];
    ty = ycos[point_index];
    tz = zcos[point_index];
    bary_cord = barycoords(v1, v2, v3, tx, ty, tz);
    integral[point_index] += interp_eval_sbb(interptargets, bary_cord[0], bary_cord[1], bary_cord[2], run_information.interp_degree);
  }
}

void cc_interaction_inverse_laplacian_icos(const RunConfig& run_information, const int index_target, const int index_source, const std::vector<IcosPanel>& icos_panels,
                    const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos,
                    const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral) {
  // computes cluster cluster interaction
  int point_index, info, dim = run_information.interp_point_count, count_target = icos_panels[index_target].point_count, count_source = icos_panels[index_source].point_count;
  std::vector<double> v1, v2, v3, v1s, v2s, v3s, func_vals(dim * dim, 0), potential_val (dim, 0), bary_cord, target_particle, source_particle,
      proxy_weights(dim, 0), basis_vals, interptargets(dim, 0), source_interp_matrix(dim * dim, 0), target_interp_matrix(dim * dim, 0),
      sxs (count_source, 0), sys (count_source, 0), szs (count_source, 0), areas (count_source, 0), pots (count_source, 0), txs (count_target, 0), tys (count_target, 0), tzs (count_target, 0);
  double u, v, w, us, vs, ws, vor, func_val, cx, cy, cz, scalar, ssh, sx, sy, sz, tx, ty, tz, pot, cxt, cyt, czt, cxs, cys, czs;
  std::vector<std::vector<double>> proxy_source_points(dim, std::vector<double>(3, 0)), proxy_target_points(dim, std::vector<double>(3, 0)),
      source_interp_points(dim, std::vector<double>(3, 0)), target_interp_points(dim, std::vector<double>(3, 0));

  fekete_init(source_interp_points, run_information.interp_degree);
  fekete_init(target_interp_points, run_information.interp_degree);
  v1s = icos_panels[index_source].vertex_1;
  v2s = icos_panels[index_source].vertex_2;
  v3s = icos_panels[index_source].vertex_3;
  v1 = icos_panels[index_target].vertex_1;
  v2 = icos_panels[index_target].vertex_2;
  v3 = icos_panels[index_target].vertex_3;

  for (int i = 0; i < count_target; i++) {
    point_index = icos_panels[index_target].points_inside[i];
    txs[i] = xcos[point_index];
    tys[i] = ycos[point_index];
    tzs[i] = zcos[point_index];
  }

  for (int j = 0; j < count_source; j++) {
    point_index = icos_panels[index_source].points_inside[j];
    sxs[j] = xcos[point_index];
    sys[j] = ycos[point_index];
    szs[j] = zcos[point_index];
    areas[j] = area[point_index];
    pots[j] = potential[point_index];
  }

  for (int i = 0; i < count_source; i++) { // compute proxy weights
    point_index = icos_panels[index_source].points_inside[i];
    sx = sxs[i];
    sy = sys[i];
    sz = szs[i];
    pot = pots[point_index];
    bary_cord = barycoords(v1s, v2s, v3s, sx, sy, sz);
    basis_vals = interp_vals_sbb(bary_cord[0], bary_cord[1], bary_cord[2], run_information.interp_degree);
    for (int j = 0; j < dim; j++) {
      proxy_weights[j] += basis_vals[j] * pot * area[point_index];
    }
  }

  for (int i = 0; i < dim; i++) { // set up source interpolation matrix
    us = source_interp_points[i][0];
    vs = source_interp_points[i][1];
    ws = 1.0 - us - vs;
    cx = us * v1s[0] + vs * v2s[0] + ws * v3s[0];
    cy = us * v1s[1] + vs * v2s[1] + ws * v3s[1];
    cz = us * v1s[2] + vs * v2s[2] + ws * v3s[2];
    scalar = run_information.radius / sqrt(cx * cx + cy * cy + cz * cz);
    cx *= scalar;
    cy *= scalar;
    cz *= scalar;
    bary_cord = barycoords(v1s, v2s, v3s, cx, cy, cz);
    source_interp_points[i]=bary_cord;
    proxy_source_points[i][0] = cx;
    proxy_source_points[i][1] = cy;
    proxy_source_points[i][2] = cz;
  }

  interp_mat_init_sbb(source_interp_matrix, source_interp_points, run_information.interp_degree, dim);

  for (int i = 0; i < dim; i++) {
    u = target_interp_points[i][0];
    v = target_interp_points[i][1];
    w = 1.0 - u - v;
    cx = u * v1[0] + v * v2[0] + w * v3[0];
    cy = u * v1[1] + v * v2[1] + w * v3[1];
    cz = u * v1[2] + v * v2[2] + w * v3[2];
    scalar = run_information.radius / sqrt(cx * cx + cy * cy + cz * cz);
    cx *= scalar;
    cy *= scalar;
    cz *= scalar;
    bary_cord = barycoords(v1, v2, v3, cx, cy, cz);
    target_interp_points[i] = bary_cord;
    proxy_target_points[i][0] = cx;
    proxy_target_points[i][1] = cy;
    proxy_target_points[i][2] = cz;
  }

  interp_mat_init_sbb(target_interp_matrix, target_interp_points, run_information.interp_degree, dim);

  for (int i = 0; i < dim; i++) { // loop over proxy target particles
    target_particle = proxy_target_points[i];
    cxt = target_particle[0];
    cyt = target_particle[1];
    czt = target_particle[2];
    for (int j = 0; j < dim; j++) { // loop over proxy source particles
      source_particle = proxy_source_points[j];
      cxs = source_particle[0];
      cys = source_particle[1];
      czs = source_particle[2];
      func_vals[dim*i+j] = -1.0/(4.0*M_PI)*log(1-cxs*cxt-cys*cyt-czs*czt);
    }
  }

  info = linear_solve(source_interp_matrix, func_vals, dim, dim, 3);
  if (info > 0) {
    throw std::runtime_error("Error with linear solve in cc source computation");
  }

  for (int i = 0; i < dim; i++) { // do PC interaction with proxy target points
    for (int j = 0; j < dim; j++) {
      potential_val[i] += func_vals[dim * i + j] * proxy_weights[j];
    }
  }

  info = linear_solve(target_interp_matrix, potential_val, dim, 1, 3);
  if (info > 0) {
    throw std::runtime_error("Error with linear solve in cc target computation");
  }

  for (int i = 0; i < count_target; i++) {
    // point_index = icos_panels[index_target].points_inside[i];
    tx = txs[i];
    ty = tys[i];
    tz = tzs[i];
    bary_cord = barycoords(v1, v2, v3, tx, ty, tz);
    integral[point_index] += interp_eval_sbb(potential_val, bary_cord[0], bary_cord[1], bary_cord[2], run_information.interp_degree);
  }
}

void fast_sum_inverse_laplacian_icos(const RunConfig& run_information, const std::vector<InteractPair>& interactions, const std::vector<IcosPanel>& icos_panels,
                                        const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos,
                                        const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral) {
  // computes fast sum with icos panels
  for (int i = 0; i < interactions.size(); i++) {
    if (interactions[i].interact_type == 0) {
      // PP
      pp_interaction_inverse_laplacian_icos(run_information, interactions[i].index_target, interactions[i].index_source, icos_panels, xcos, ycos, zcos, area, potential, integral);
    } else if (interactions[i].interact_type == 1) {
      // PC
      pc_interaction_inverse_laplacian_icos(run_information, interactions[i].index_target, interactions[i].index_source, icos_panels, xcos, ycos, zcos, area, potential, integral);
    } else if (interactions[i].interact_type == 2) {
      // CP
      cp_interaction_inverse_laplacian_icos(run_information, interactions[i].index_target, interactions[i].index_source, icos_panels, xcos, ycos, zcos, area, potential, integral);
    } else {
      // CC
      cc_interaction_inverse_laplacian_icos(run_information, interactions[i].index_target, interactions[i].index_source, icos_panels, xcos, ycos, zcos, area, potential, integral);
    }
  }
}

void pp_interaction_inverse_laplacian_cube(const RunConfig& run_information, const int index_target, const int index_source, const std::vector<CubePanel>& cube_panels,
                    const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos,
                    const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral) {
  // compute particle particle interaction
  int target_i, source_j, source_count;
  double tx, ty, tz, sx, sy, sz, gfval;

  source_count = cube_panels[index_source].point_count;

  std::vector<double> sxs (source_count, 0), sys (source_count, 0), szs (source_count, 0), pots (source_count, 0), areas (source_count, 0);

  for (int j = 0; j < source_count; j++) {
    source_j = cube_panels[index_source].points_inside[j];
    sxs[j] = xcos[source_j];
    sys[j] = ycos[source_j];
    szs[j] = zcos[source_j];
    pots[j] = potential[source_j];
    areas[j] = area[source_j];
  }

  for (int i = 0; i < cube_panels[index_target].point_count; i++) {
    target_i = cube_panels[index_target].points_inside[i];
    tx = xcos[target_i];
    ty = ycos[target_i];
    tz = zcos[target_i];
    for (int j = 0; j < source_count; j++) {
      if (target_i != cube_panels[index_source].points_inside[j]) {
        sx = sxs[j], sy = sys[j], sz = szs[j];
        gfval = -1.0/(4.0*M_PI)*log(1-tx*sx-ty*sy-tz*sz);
        integral[target_i] += gfval * areas[j] * pots[j];
      }
    }
  }
}

void pc_interaction_inverse_laplacian_cube(const RunConfig& run_information, const int index_target, const int index_source, const std::vector<CubePanel>& cube_panels,
                    const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos,
                    const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral) {
  // pc interaction
  int point_index, info, dim = run_information.interp_point_count, count_target = cube_panels[index_target].point_count, count_source = cube_panels[index_source].point_count, degree=run_information.interp_degree;
  std::vector<std::vector<double>> proxy_weights (degree+1, std::vector<double> (degree+1, 0)), basis_vals, func_vals (degree+1, std::vector<double> (degree+1, 0));
  std::vector<double> txs (count_target, 0), tys (count_target, 0), tzs (count_target, 0), xieta, cheb_xi, cheb_eta, xyz;
  double pot, tx, ty, tz, cx, cy, cz, sx, sy, sz, xi, eta;

  for (int i = 0; i < count_target; i++) {
    point_index = cube_panels[index_target].points_inside[i];
    txs[i] = xcos[point_index];
    tys[i] = ycos[point_index];
    tzs[i] = zcos[point_index];
  }

  cheb_xi = bli_interp_points_shift(cube_panels[index_source].min_xi, cube_panels[index_source].max_xi, degree);
  cheb_eta = bli_interp_points_shift(cube_panels[index_source].min_eta, cube_panels[index_source].max_eta, degree);

  // compute proxy weights
  for (int i = 0; i < count_source; i++) {
    // loop over source particles
    point_index = cube_panels[index_source].points_inside[i];
    sx = xcos[point_index], sy = ycos[point_index], sz = zcos[point_index];
    xieta = xieta_from_xyz(sx, sy, sz, cube_panels[index_source].face);
    basis_vals = interp_vals_bli(xieta[0], xieta[1], cube_panels[index_source].min_xi, cube_panels[index_source].max_xi, cube_panels[index_source].min_eta, cube_panels[index_source].max_eta, degree);
    for (int j = 0; j < degree+1; j++) { // loop over xi
      for (int k = 0; k < degree+1; k++) { // loop over eta
        proxy_weights[j][k] += basis_vals[j][k] * area[point_index] * potential[point_index];
      }
    }
  }

  // interpolate to target points
  for (int i = 0; i < count_target; i++) {
    // loop over target particles
    tx = txs[i], ty = tys[i], tz = tzs[i];
    point_index = cube_panels[index_target].points_inside[i];
    for (int j = 0; j < degree+1; j++) { // xi loop
      xi = cheb_xi[j];
      for (int k = 0; k < degree+1; k++) { // eta loop
        eta = cheb_eta[k];
        xyz = xyz_from_xieta(xi, eta, cube_panels[index_source].face);
        integral[point_index] += -1.0/(4.0*M_PI)*log(1-tx*xyz[0]-ty*xyz[1]-tz*xyz[2]) * proxy_weights[j][k];
      }
    }
  }
}

void cp_interaction_inverse_laplacian_cube(const RunConfig& run_information, const int index_target, const int index_source, const std::vector<CubePanel>& cube_panels,
                    const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos,
                    const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral) {
  // cp interaction
  int point_index, count_target = cube_panels[index_target].point_count, count_source = cube_panels[index_source].point_count, degree=run_information.interp_degree;
  std::vector<double> cheb_xi, cheb_eta, sxs (count_source, 0), sys (count_source, 0), szs (count_source, 0), areas (count_source, 0), pots (count_source, 0), xyz, xieta;
  std::vector<std::vector<double>> func_points (degree+1, std::vector<double> (degree+1, 0)), basis_vals;
  double sx, sy, sz, tx, ty, tz;

  for (int j = 0; j < count_source; j++) {
    point_index = cube_panels[index_source].points_inside[j];
    sxs[j] = xcos[point_index];
    sys[j] = ycos[point_index];
    szs[j] = zcos[point_index];
    areas[j] = area[point_index];
    pots[j] = potential[point_index];
  }

  // compute func vals with proxy target points
  cheb_xi = bli_interp_points_shift(cube_panels[index_target].min_xi, cube_panels[index_target].max_xi, degree);
  cheb_eta = bli_interp_points_shift(cube_panels[index_target].min_eta, cube_panels[index_target].max_eta, degree);
  for (int i = 0; i < degree+1; i++) { // xi loop
    for (int j = 0; j < degree+1; j++) { // eta loop
      xyz = xyz_from_xieta(cheb_xi[i], cheb_eta[j], cube_panels[index_target].face);
      for (int k = 0; k < count_source; k++) {
        // loop over source points
        sx = sxs[k], sy = sys[k], sz = szs[k];
        func_points[i][j] += -1.0/(4.0*M_PI)*log(1-xyz[0]*sx-xyz[1]*sy-xyz[2]*sz) * areas[k]*pots[k];
      }
    }
  }

  // interpolate from proxy target points to target points
  for (int i = 0; i < count_target; i++) {
    point_index = cube_panels[index_target].points_inside[i];
    tx = xcos[point_index];
    ty = ycos[point_index];
    tz = zcos[point_index];
    xieta = xieta_from_xyz(tx, ty, tz, cube_panels[index_target].face);
    basis_vals = interp_vals_bli(xieta[0], xieta[1], cube_panels[index_target].min_xi, cube_panels[index_target].max_xi, cube_panels[index_target].min_eta, cube_panels[index_target].max_eta, degree);
    for (int j = 0; j < degree+1; j++) {
      for (int k = 0; k < degree+1; k++) {
        integral[point_index] += basis_vals[j][k] * func_points[j][k];
      }
    }
  }
}

void cc_interaction_inverse_laplacian_cube(const RunConfig& run_information, const int index_target, const int index_source, const std::vector<CubePanel>& cube_panels,
                    const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos,
                    const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral) {
  // cc interaction
  int point_index, count_target = cube_panels[index_target].point_count, count_source = cube_panels[index_source].point_count, degree=run_information.interp_degree;
  std::vector<std::vector<double>> proxy_weights (degree+1, std::vector<double> (degree+1, 0)), basis_vals, func_vals (degree+1, std::vector<double> (degree+1, 0)), func_points (degree+1, std::vector<double> (degree+1, 0));
  std::vector<double> txs (count_target, 0), tys (count_target, 0), tzs (count_target, 0), sxs (count_source, 0), sys (count_source, 0), szs (count_source, 0), areas (count_source, 0), pots (count_source, 0), xieta, cheb_xi_t, cheb_eta_t, cheb_xi_s, cheb_eta_s, xyz_t, xyz_s;
  double pot, tx, ty, tz, cxs, cys, czs, cxt, cyt, czt, sx, sy, sz, xi_t, eta_t, xi_s, eta_s;

  for (int i = 0; i < count_target; i++) {
    point_index = cube_panels[index_target].points_inside[i];
    txs[i] = xcos[point_index];
    tys[i] = ycos[point_index];
    tzs[i] = zcos[point_index];
  }

  for (int j = 0; j < count_source; j++) {
    point_index = cube_panels[index_source].points_inside[j];
    sxs[j] = xcos[point_index];
    sys[j] = ycos[point_index];
    szs[j] = zcos[point_index];
    areas[j] = area[point_index];
    pots[j] = potential[point_index];
  }

  cheb_xi_s = bli_interp_points_shift(cube_panels[index_source].min_xi, cube_panels[index_source].max_xi, degree);
  cheb_eta_s = bli_interp_points_shift(cube_panels[index_source].min_eta, cube_panels[index_source].max_eta, degree);
  cheb_xi_t = bli_interp_points_shift(cube_panels[index_target].min_xi, cube_panels[index_target].max_xi, degree);
  cheb_eta_t = bli_interp_points_shift(cube_panels[index_target].min_eta, cube_panels[index_target].max_eta, degree);

  // compute proxy weights
  for (int i = 0; i < count_source; i++) {
    // loop over source particles
    point_index = cube_panels[index_source].points_inside[i];
    sx = xcos[point_index], sy = ycos[point_index], sz = zcos[point_index];
    xieta = xieta_from_xyz(sx, sy, sz, cube_panels[index_source].face);
    basis_vals = interp_vals_bli(xieta[0], xieta[1], cube_panels[index_source].min_xi, cube_panels[index_source].max_xi, cube_panels[index_source].min_eta, cube_panels[index_source].max_eta, degree);
    for (int j = 0; j < degree+1; j++) { // loop over xi
      for (int k = 0; k < degree+1; k++) { // loop over eta
        proxy_weights[j][k] += basis_vals[j][k] * area[point_index] * potential[point_index];
      }
    }
  }

  // computes func vals at proxy target points from proxy source points
  for (int i = 0; i < degree+1; i++) { // target xi
    xi_t = cheb_xi_t[i];
    for (int j = 0; j < degree+1; j++) { // target eta
      eta_t = cheb_eta_t[j];
      xyz_t = xyz_from_xieta(xi_t, eta_t, cube_panels[index_target].face);
      cxt = xyz_t[0], cyt = xyz_t[1], czt = xyz_t[2];
      for (int k = 0; k < degree+1; k++) { // source xi
        xi_s = cheb_xi_s[k];
        for (int l = 0; l < degree+1; l++) { // source eta
          eta_s = cheb_eta_s[l];
          xyz_s = xyz_from_xieta(xi_s, eta_s, cube_panels[index_source].face);
          cxs = xyz_s[0], cys = xyz_s[1], czs = xyz_s[2];
          func_points[i][j] += -1.0/(4.0*M_PI)*log(1-cxs*cxt-cys*cyt-czs*czt)*proxy_weights[k][l];
        }
      }
    }
  }

  for (int i = 0; i < count_target; i++) {
    point_index = cube_panels[index_target].points_inside[i];
    tx = xcos[point_index];
    ty = ycos[point_index];
    tz = zcos[point_index];
    xieta = xieta_from_xyz(tx, ty, tz, cube_panels[index_target].face);
    basis_vals = interp_vals_bli(xieta[0], xieta[1], cube_panels[index_target].min_xi, cube_panels[index_target].max_xi, cube_panels[index_target].min_eta, cube_panels[index_target].max_eta, degree);
    for (int j = 0; j < degree+1; j++) {
      for (int k = 0; k < degree+1; k++) {
        integral[point_index] += basis_vals[j][k] * func_points[j][k];
      }
    }
  }
}

void fast_sum_inverse_laplacian_cube(const RunConfig& run_information, const std::vector<InteractPair>& interactions, const std::vector<CubePanel>& cube_panels,
                                        const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos,
                                        const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral) {
  // computes fast sum with cubed sphere panels
  for (int i = 0; i < interactions.size(); i++) {
    if (interactions[i].interact_type == 0) {
      // PP
      pp_interaction_inverse_laplacian_cube(run_information, interactions[i].index_target, interactions[i].index_source, cube_panels, xcos, ycos, zcos, area, potential, integral);
    } else if (interactions[i].interact_type == 1) {
      // PC
      pc_interaction_inverse_laplacian_cube(run_information, interactions[i].index_target, interactions[i].index_source, cube_panels, xcos, ycos, zcos, area, potential, integral);
    } else if (interactions[i].interact_type == 2) {
      // CP
      cp_interaction_inverse_laplacian_cube(run_information, interactions[i].index_target, interactions[i].index_source, cube_panels, xcos, ycos, zcos, area, potential, integral);
    } else {
      // CC
      cc_interaction_inverse_laplacian_cube(run_information, interactions[i].index_target, interactions[i].index_source, cube_panels, xcos, ycos, zcos, area, potential, integral);
    }
  }
}
