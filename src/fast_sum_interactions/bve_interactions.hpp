#ifndef H_FAST_SUM_BVE_H
#define H_FAST_SUM_BVE_H

#include <vector>

#include "../structs.hpp"

void pp_interaction_bve(const RunConfig& run_information, const CubePanel& cube_panel_source, const CubePanel& cube_panel_target,
                    const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos,
                    const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral_1, 
                    std::vector<double>& integral_2, std::vector<double>& integral_3);

void pc_interaction_bve(const RunConfig& run_information, const CubePanel& cube_panel_source, const CubePanel& cube_panel_target,
                    const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos,
                    const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral_1, 
                    std::vector<double>& integral_2, std::vector<double>& integral_3);

void cp_interaction_bve(const RunConfig& run_information, const CubePanel& cube_panel_source, const CubePanel& cube_panel_target,
                    const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos,
                    const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral_1, 
                    std::vector<double>& integral_2, std::vector<double>& integral_3);

void cc_interaction_bve(const RunConfig& run_information, const CubePanel& cube_panel_source, const CubePanel& cube_panel_target,
                    const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos,
                    const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral_1, 
                    std::vector<double>& integral_2, std::vector<double>& integral_3);
#endif