#include <vector>
#include <cmath>
#include <iostream>

#include "interp_utils.hpp"
#include "general_utils.hpp"
#include "structs.hpp"
#include "./fast_sum_interactions/inverse_laplacian_interactions.hpp"
#include "./fast_sum_interactions/inverse_biharmonic_interactions.hpp"
#include "./fast_sum_interactions/laplacian_interactions.hpp"
#include "./fast_sum_interactions/sal_interactions.hpp"
#include "./fast_sum_interactions/sal_interactions_lat_deriv.hpp"
#include "./fast_sum_interactions/sal_interactions_lon_deriv.hpp"

void fast_sum_inverse_laplacian(const RunConfig& run_information, const std::vector<InteractPair>& interactions, const std::vector<CubePanel>& cube_panels,
                                        const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos,
                                        const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral) {
  // computes fast sum with cubed sphere panels
  int i_t, i_s;
  for (int i = 0; i < interactions.size(); i++) {
    if (i % run_information.mpi_P == run_information.mpi_ID) {
      i_t = interactions[i].index_target;
      i_s = interactions[i].index_source;
      if (interactions[i].interact_type == 0) { // PP
        pp_interaction_inverse_laplacian(run_information, cube_panels[i_s], cube_panels[i_t], xcos, ycos, zcos, area, potential, integral);
      } else if (interactions[i].interact_type == 1) { // PC
        pc_interaction_inverse_laplacian(run_information, cube_panels[i_s], cube_panels[i_t], xcos, ycos, zcos, area, potential, integral);
      } else if (interactions[i].interact_type == 2) { // CP
        cp_interaction_inverse_laplacian(run_information, cube_panels[i_s], cube_panels[i_t], xcos, ycos, zcos, area, potential, integral);
      } else { // CC
        cc_interaction_inverse_laplacian(run_information, cube_panels[i_s], cube_panels[i_t], xcos, ycos, zcos, area, potential, integral);
      }
    }
  }
}

void fast_sum_inverse_biharmonic(const RunConfig& run_information, const std::vector<InteractPair>& interactions, const std::vector<CubePanel>& cube_panels_source, const std::vector<CubePanel>& cube_panels_target,
                                        const std::vector<double>& xcos_t, const std::vector<double>& ycos_t, const std::vector<double>& zcos_t,
                                        const std::vector<double>& xcos_s, const std::vector<double>& ycos_s, const std::vector<double>& zcos_s,
                                        const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral) {
  // computes fast sum with cubed sphere panels
  int i_t, i_s;
  for (int i = 0; i < interactions.size(); i++) {
    if (i % run_information.mpi_P == run_information.mpi_ID) {
      i_t = interactions[i].index_target;
      i_s = interactions[i].index_source;
      if (interactions[i].interact_type == 0) { // PP
        pp_interaction_inverse_biharmonic(run_information, cube_panels_source[i_s], cube_panels_target[i_t], xcos_t, ycos_t, zcos_t, xcos_s, ycos_s, zcos_s, area, potential, integral);
      } else if (interactions[i].interact_type == 1) { // PC
        // pp_interaction_inverse_biharmonic(run_information, cube_panels_source[i_s], cube_panels_target[i_t], xcos_t, ycos_t, zcos_t, xcos_s, ycos_s, zcos_s, area, potential, integral);
        pc_interaction_inverse_biharmonic(run_information, cube_panels_source[i_s], cube_panels_target[i_t], xcos_t, ycos_t, zcos_t, xcos_s, ycos_s, zcos_s, area, potential, integral);
      } else if (interactions[i].interact_type == 2) { // CP
        // pp_interaction_inverse_biharmonic(run_information, cube_panels_source[i_s], cube_panels_target[i_t], xcos_t, ycos_t, zcos_t, xcos_s, ycos_s, zcos_s, area, potential, integral);
        cp_interaction_inverse_biharmonic(run_information, cube_panels_source[i_s], cube_panels_target[i_t], xcos_t, ycos_t, zcos_t, xcos_s, ycos_s, zcos_s, area, potential, integral);
      } else { // CC
        // pp_interaction_inverse_biharmonic(run_information, cube_panels_source[i_s], cube_panels_target[i_t], xcos_t, ycos_t, zcos_t, xcos_s, ycos_s, zcos_s, area, potential, integral);
        // pc_interaction_inverse_biharmonic(run_information, cube_panels_source[i_s], cube_panels_target[i_t], xcos_t, ycos_t, zcos_t, xcos_s, ycos_s, zcos_s, area, potential, integral);
        cc_interaction_inverse_biharmonic(run_information, cube_panels_source[i_s], cube_panels_target[i_t], xcos_t, ycos_t, zcos_t, xcos_s, ycos_s, zcos_s, area, potential, integral);
      }
    }
  }
}

void fast_sum_laplacian(const RunConfig& run_information, const std::vector<InteractPair>& interactions, const std::vector<CubePanel>& cube_panels,
                                        const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos,
                                        const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral) {
  // computes fast sum with cubed sphere panels
  int i_t, i_s;
  for (int i = 0; i < interactions.size(); i++) {
    if (i % run_information.mpi_P == run_information.mpi_ID) {
      i_t = interactions[i].index_target;
      i_s = interactions[i].index_source;
      if (interactions[i].interact_type == 0) { // PP
        pp_interaction_laplacian(run_information, cube_panels[i_s], cube_panels[i_t], xcos, ycos, zcos, area, potential, integral);
      } else if (interactions[i].interact_type == 1) { // PC
        pc_interaction_laplacian(run_information, cube_panels[i_s], cube_panels[i_t], xcos, ycos, zcos, area, potential, integral);
      } else if (interactions[i].interact_type == 2) { // CP
        // cp_interaction_laplacian(run_information, cube_panels[i_s], cube_panels[i_t], xcos, ycos, zcos, area, potential, integral);
        pp_interaction_laplacian(run_information, cube_panels[i_s], cube_panels[i_t], xcos, ycos, zcos, area, potential, integral);
      } else { // CC
        // cc_interaction_laplacian(run_information, cube_panels[i_s], cube_panels[i_t], xcos, ycos, zcos, area, potential, integral);
        pc_interaction_laplacian(run_information, cube_panels[i_s], cube_panels[i_t], xcos, ycos, zcos, area, potential, integral);
      }
    }
  }
}

void fast_sum_sal(const RunConfig& run_information, const std::vector<InteractPair>& interactions, const std::vector<CubePanel>& cube_panels,
                                        const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos,
                                        const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral) {
  // computes fast sum with cubed sphere panels
  int i_t, i_s;
  for (int i = 0; i < interactions.size(); i++) {
    if (i % run_information.mpi_P == run_information.mpi_ID) {
      i_t = interactions[i].index_target;
      i_s = interactions[i].index_source;
      if (interactions[i].interact_type == 0) { // PP
        pp_interaction_sal(run_information, cube_panels[i_s], cube_panels[i_t], xcos, ycos, zcos, area, potential, integral);
      } else if (interactions[i].interact_type == 1) { // PC
        pc_interaction_sal(run_information, cube_panels[i_s], cube_panels[i_t], xcos, ycos, zcos, area, potential, integral);
      } else if (interactions[i].interact_type == 2) { // CP
        cp_interaction_sal(run_information, cube_panels[i_s], cube_panels[i_t], xcos, ycos, zcos, area, potential, integral);
      } else { // CC
        cc_interaction_sal(run_information, cube_panels[i_s], cube_panels[i_t], xcos, ycos, zcos, area, potential, integral);
      }
    }
  }
}

void fast_sum_sal_lat_deriv(const RunConfig& run_information, const std::vector<InteractPair>& interactions, const std::vector<CubePanel>& cube_panels,
                                        const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos,
                                        const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral) {
  // computes fast sum with cubed sphere panels
  int i_t, i_s;
  for (int i = 0; i < interactions.size(); i++) {
    if (i % run_information.mpi_P == run_information.mpi_ID) {
      i_t = interactions[i].index_target;
      i_s = interactions[i].index_source;
      if (interactions[i].interact_type == 0) { // PP
        pp_interaction_sal_lat_deriv(run_information, cube_panels[i_s], cube_panels[i_t], xcos, ycos, zcos, area, potential, integral);
      } else if (interactions[i].interact_type == 1) { // PC
        pc_interaction_sal_lat_deriv(run_information, cube_panels[i_s], cube_panels[i_t], xcos, ycos, zcos, area, potential, integral);
      } else if (interactions[i].interact_type == 2) { // CP
        cp_interaction_sal_lat_deriv(run_information, cube_panels[i_s], cube_panels[i_t], xcos, ycos, zcos, area, potential, integral);
        // pp_interaction_sal_lat_deriv(run_information, cube_panels[i_s], cube_panels[i_t], xcos, ycos, zcos, area, potential, integral);
      } else { // CC
        cc_interaction_sal_lat_deriv(run_information, cube_panels[i_s], cube_panels[i_t], xcos, ycos, zcos, area, potential, integral);
        // pc_interaction_sal_lat_deriv(run_information, cube_panels[i_s], cube_panels[i_t], xcos, ycos, zcos, area, potential, integral);
      }
    }
  }
}

void fast_sum_sal_lon_deriv(const RunConfig& run_information, const std::vector<InteractPair>& interactions, const std::vector<CubePanel>& cube_panels,
                                        const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos,
                                        const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral) {
  // computes fast sum with cubed sphere panels
  int i_t, i_s;
  for (int i = 0; i < interactions.size(); i++) {
    if (i % run_information.mpi_P == run_information.mpi_ID) {
      i_t = interactions[i].index_target;
      i_s = interactions[i].index_source;
      if (interactions[i].interact_type == 0) { // PP
        pp_interaction_sal_lon_deriv(run_information, cube_panels[i_s], cube_panels[i_t], xcos, ycos, zcos, area, potential, integral);
      } else if (interactions[i].interact_type == 1) { // PC
        pc_interaction_sal_lon_deriv(run_information, cube_panels[i_s], cube_panels[i_t], xcos, ycos, zcos, area, potential, integral);
      } else if (interactions[i].interact_type == 2) { // CP
        cp_interaction_sal_lon_deriv(run_information, cube_panels[i_s], cube_panels[i_t], xcos, ycos, zcos, area, potential, integral);
        // pp_interaction_sal_lon_deriv(run_information, cube_panels[i_s], cube_panels[i_t], xcos, ycos, zcos, area, potential, integral);

      } else { // CC
        cc_interaction_sal_lon_deriv(run_information, cube_panels[i_s], cube_panels[i_t], xcos, ycos, zcos, area, potential, integral);
        // pc_interaction_sal_lon_deriv(run_information, cube_panels[i_s], cube_panels[i_t], xcos, ycos, zcos, area, potential, integral);

      }
    }
  }
}
