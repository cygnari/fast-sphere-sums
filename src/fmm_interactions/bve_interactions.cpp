#include <vector>
#include <cmath>
#include <iostream>

#include "../interp_utils.hpp"
#include "../general_utils.hpp"
#include "../structs.hpp"

void fmm_pp_interaction_bve(const RunConfig& run_information, const CubePanel& cube_panel_source, const CubePanel& cube_panel_target,
                    const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos,
                    const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral1, 
                    std::vector<double>& integral2, std::vector<double>& integral3) {
  // compute particle particle interaction
  int target_i, source_j, source_count;
  double tx, ty, tz, sx, sy, sz, gfval;

  source_count = cube_panel_source.point_count;

  std::vector<double> sxs (source_count, 0), sys (source_count, 0), szs (source_count, 0), pots (source_count, 0), areas (source_count, 0);

  for (int j = 0; j < source_count; j++) {
    source_j = cube_panel_source.points_inside[j];
    sxs[j] = xcos[source_j];
    sys[j] = ycos[source_j];
    szs[j] = zcos[source_j];
    pots[j] = potential[source_j];
    areas[j] = area[source_j];
  }

  for (int i = 0; i < cube_panel_target.point_count; i++) {
    target_i = cube_panel_target.points_inside[i];
    tx = xcos[target_i];
    ty = ycos[target_i];
    tz = zcos[target_i];
    for (int j = 0; j < source_count; j++) {
      if (target_i != cube_panel_source.points_inside[j]) {
        sx = sxs[j], sy = sys[j], sz = szs[j];
        gfval = -1.0/(4.0*M_PI*(1.0-tx*sx-ty*sy-tz*sz))*pots[j]*areas[j];
        integral1[target_i] += (ty*sz-tz*sy)*gfval;
        integral2[target_i] += (tz*sx-tx*sz)*gfval;
        integral3[target_i] += (tx*sy-ty*sx)*gfval;
      }
    }
  }
}

void fmm_pc_interaction_bve(const RunConfig& run_information, const CubePanel& cube_panel_source, const CubePanel& cube_panel_target,
                    const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos,
                    const std::vector<double>& area, const std::vector<std::vector<double>>& proxy_source_weights, 
                    std::vector<double>& integral1, std::vector<double>& integral2, std::vector<double>& integral3) {
  // pc interaction
  int point_index, count_target = cube_panel_target.point_count, count_source = cube_panel_source.point_count, degree=run_information.interp_degree;
  std::vector<double> cheb_xi, cheb_eta, xyz;
  double tx, ty, tz, xi, eta, gfval;

  for (int i = 0; i < count_target; i++) {
    point_index = cube_panel_target.points_inside[i];
    tx = xcos[point_index];
    ty = ycos[point_index];
    tz = zcos[point_index];
    cheb_xi = bli_interp_points_shift(cube_panel_source.min_xi, cube_panel_source.max_xi, degree);
    cheb_eta = bli_interp_points_shift(cube_panel_source.min_eta, cube_panel_source.max_eta, degree);
    for (int j = 0; j < degree+1; j++) { // xi loop
      xi = cheb_xi[j];
      for (int k = 0; k < degree+1; k++) { // eta loop
        eta = cheb_eta[k];
        xyz = xyz_from_xieta(xi, eta, cube_panel_source.face);
        gfval = -1.0/(4.0*M_PI*(1.0-tx*xyz[0]-ty*xyz[1]-tz*xyz[2]))*proxy_source_weights[j][k];
        integral1[point_index] += (ty*xyz[2]-tz*xyz[1])*gfval;
        integral2[point_index] += (tz*xyz[0]-tx*xyz[2])*gfval;
        integral3[point_index] += (tx*xyz[1]-ty*xyz[0])*gfval;
      }
    }
  }
}

void fmm_cp_interaction_bve(const RunConfig& run_information, const CubePanel& cube_panel_source, const CubePanel& cube_panel_target,
                    const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos,
                    const std::vector<double>& area, const std::vector<double>& potential, 
                    std::vector<std::vector<double>>& proxy_target_potenti1, std::vector<std::vector<double>>& proxy_target_potenti2, 
                    std::vector<std::vector<double>>& proxy_target_potenti3) {
  // cp interaction
  int point_index, count_target = cube_panel_target.point_count, count_source = cube_panel_source.point_count, degree=run_information.interp_degree;
  std::vector<double> cheb_xi, cheb_eta, sxs (count_source, 0), sys (count_source, 0), szs (count_source, 0), areas (count_source, 0), pots (count_source, 0), xyz;
  double sx, sy, sz, gfval;

  for (int j = 0; j < count_source; j++) {
    point_index = cube_panel_source.points_inside[j];
    sxs[j] = xcos[point_index];
    sys[j] = ycos[point_index];
    szs[j] = zcos[point_index];
    areas[j] = area[point_index];
    pots[j] = potential[point_index];
  }

  // compute func vals with proxy target points
  cheb_xi = bli_interp_points_shift(cube_panel_target.min_xi, cube_panel_target.max_xi, degree);
  cheb_eta = bli_interp_points_shift(cube_panel_target.min_eta, cube_panel_target.max_eta, degree);
  for (int i = 0; i < degree+1; i++) { // xi loop
    for (int j = 0; j < degree+1; j++) { // eta loop
      xyz = xyz_from_xieta(cheb_xi[i], cheb_eta[j], cube_panel_target.face);
      for (int k = 0; k < count_source; k++) {
        // loop over source points
        sx = sxs[k], sy = sys[k], sz = szs[k];
        gfval = -1.0/(4.0*M_PI*(1.0-xyz[0]*sx-xyz[1]*sy-xyz[2]*sz))*pots[k]*areas[k];
        proxy_target_potenti1[i][j] += (xyz[1]*sz-xyz[2]*sy)*gfval;
        proxy_target_potenti2[i][j] += (xyz[2]*sx-xyz[0]*sz)*gfval;
        proxy_target_potenti3[i][j] += (xyz[0]*sy-xyz[1]*sx)*gfval;
      }
    }
  }
}

void fmm_cc_interaction_bve(const RunConfig& run_information, const CubePanel& cube_panel_source, const CubePanel& cube_panel_target,
                    const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos,
                    const std::vector<double>& area, const std::vector<std::vector<double>>& proxy_source_weights, 
                    std::vector<std::vector<double>>& proxy_target_potenti1, std::vector<std::vector<double>>& proxy_target_potenti2, 
                    std::vector<std::vector<double>>& proxy_target_potenti3) {
  // cc interaction
  int point_index, count_target = cube_panel_target.point_count, count_source = cube_panel_source.point_count, degree=run_information.interp_degree;
  std::vector<std::vector<double>> proxy_weights (degree+1, std::vector<double> (degree+1, 0)), basis_vals, func_vals (degree+1, std::vector<double> (degree+1, 0)), func_points (degree+1, std::vector<double> (degree+1, 0));
  std::vector<double> txs (count_target, 0), tys (count_target, 0), tzs (count_target, 0), sxs (count_source, 0), sys (count_source, 0), szs (count_source, 0), areas (count_source, 0), pots (count_source, 0), xieta, cheb_xi_t, cheb_eta_t, cheb_xi_s, cheb_eta_s, xyz_t, xyz_s;
  double pot, tx, ty, tz, cxs, cys, czs, cxt, cyt, czt, sx, sy, sz, xi_t, eta_t, xi_s, eta_s, gfval;

  cheb_xi_s = bli_interp_points_shift(cube_panel_source.min_xi, cube_panel_source.max_xi, degree);
  cheb_eta_s = bli_interp_points_shift(cube_panel_source.min_eta, cube_panel_source.max_eta, degree);
  cheb_xi_t = bli_interp_points_shift(cube_panel_target.min_xi, cube_panel_target.max_xi, degree);
  cheb_eta_t = bli_interp_points_shift(cube_panel_target.min_eta, cube_panel_target.max_eta, degree);

  // computes func vals at proxy target points from proxy source points
  for (int i = 0; i < degree+1; i++) { // target xi
    xi_t = cheb_xi_t[i];
    for (int j = 0; j < degree+1; j++) { // target eta
      eta_t = cheb_eta_t[j];
      xyz_t = xyz_from_xieta(xi_t, eta_t, cube_panel_target.face);
      cxt = xyz_t[0], cyt = xyz_t[1], czt = xyz_t[2];
      for (int k = 0; k < degree+1; k++) { // source xi
        xi_s = cheb_xi_s[k];
        for (int l = 0; l < degree+1; l++) { // source eta
          eta_s = cheb_eta_s[l];
          xyz_s = xyz_from_xieta(xi_s, eta_s, cube_panel_source.face);
          cxs = xyz_s[0], cys = xyz_s[1], czs = xyz_s[2];
          gfval = -1.0/(4.0*M_PI*(1.0-cxs*cxt-cys*cyt-czs*czt))*proxy_source_weights[k][l];
          // proxy_target_potenti[i][j] += -1.0/(4.0*M_PI)*log(1-cxs*cxt-cys*cyt-czs*czt)*proxy_source_weights[k][l];
          proxy_target_potenti1[i][j] += (cyt*czs-czt*cys)*gfval;
          proxy_target_potenti2[i][j] += (czt*cxs-cxt*czs)*gfval;
          proxy_target_potenti3[i][j] += (cxt*cys-cyt*cxs)*gfval;
        }
      }
    }
  }
}
