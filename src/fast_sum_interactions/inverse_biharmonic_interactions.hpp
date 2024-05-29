#ifndef H_FAST_SUM_MIN_CURVE_INTERP_H
#define H_FAST_SUM_MIN_CURVE_INTERP_H

#include <vector>

#include "../structs.hpp"

void pp_interaction_inverse_biharmonic(const RunConfig& run_information, const CubePanel& cube_panel_source, const CubePanel& cube_panel_target,
                    const std::vector<double>& xcos_t, const std::vector<double>& ycos_t, const std::vector<double>& zcos_t,
                    const std::vector<double>& xcos_s, const std::vector<double>& ycos_s, const std::vector<double>& zcos_s,
                    const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral);

void pc_interaction_inverse_biharmonic(const RunConfig& run_information, const CubePanel& cube_panel_source, const CubePanel& cube_panel_target,
                    const std::vector<double>& xcos_t, const std::vector<double>& ycos_t, const std::vector<double>& zcos_t,
                    const std::vector<double>& xcos_s, const std::vector<double>& ycos_s, const std::vector<double>& zcos_s,
                    const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral);

void cp_interaction_inverse_biharmonic(const RunConfig& run_information, const CubePanel& cube_panel_source, const CubePanel& cube_panel_target,
                    const std::vector<double>& xcos_t, const std::vector<double>& ycos_t, const std::vector<double>& zcos_t,
                    const std::vector<double>& xcos_s, const std::vector<double>& ycos_s, const std::vector<double>& zcos_s,
                    const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral);

void cc_interaction_inverse_biharmonic(const RunConfig& run_information, const CubePanel& cube_panel_source, const CubePanel& cube_panel_target,
                    const std::vector<double>& xcos_t, const std::vector<double>& ycos_t, const std::vector<double>& zcos_t,
                    const std::vector<double>& xcos_s, const std::vector<double>& ycos_s, const std::vector<double>& zcos_s,
                    const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral);
#endif
