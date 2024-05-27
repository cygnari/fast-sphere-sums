#include <vector>
#include <cmath>
#include <iostream>

#include "interp_utils.hpp"
#include "general_utils.hpp"
#include "structs.hpp"
#include "./fast_sum_interactions/inverse_laplacian_interactions.hpp"

void fast_sum_inverse_laplacian_cube(const RunConfig& run_information, const std::vector<InteractPair>& interactions, const std::vector<CubePanel>& cube_panels,
                                        const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos,
                                        const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral) {
  // computes fast sum with cubed sphere panels
  for (int i = 0; i < interactions.size(); i++) {
    if (i % run_information.mpi_P == run_information.mpi_ID) {
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
}
