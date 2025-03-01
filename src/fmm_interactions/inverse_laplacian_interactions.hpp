#ifndef H_FAST_SUM_INVERSE_LAPLACIAN_H
#define H_FAST_SUM_INVERSE_LAPLACIAN_H

#include <vector>

#include "../structs.hpp"

void fmm_pp_interaction_inverse_laplacian(const RunConfig& run_information, const CubePanel& cube_panel_source, const CubePanel& cube_panel_target,
                    const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos,
                    const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral);

void fmm_pc_interaction_inverse_laplacian(const RunConfig& run_information, const CubePanel& cube_panel_source, const CubePanel& cube_panel_target,
                    const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos,
                    const std::vector<double>& area, const std::vector<std::vector<double>>& proxy_source_weights, 
                    std::vector<double>& integral);

void fmm_cp_interaction_inverse_laplacian(const RunConfig& run_information, const CubePanel& cube_panel_source, const CubePanel& cube_panel_target,
                    const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos,
                    const std::vector<double>& area, const std::vector<double>& potential, 
                    std::vector<std::vector<double>>& proxy_target_potenti);

void fmm_cc_interaction_inverse_laplacian(const RunConfig& run_information, const CubePanel& cube_panel_source, const CubePanel& cube_panel_target,
                    const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos,
                    const std::vector<double>& area, const std::vector<std::vector<double>>& proxy_source_weights, 
                    std::vector<std::vector<double>>& proxy_target_potenti);
#endif
