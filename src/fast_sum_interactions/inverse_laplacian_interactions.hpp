#ifndef H_FAST_SUM_INVERSE_LAPLACIAN_H
#define H_FAST_SUM_INVERSE_LAPLACIAN_H

#include <vector>

#include "../structs.hpp"

void pp_interaction_inverse_laplacian_cube(const RunConfig& run_information, const int index_target, const int index_source, const std::vector<CubePanel>& cube_panels,
                    const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos,
                    const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral);

void pc_interaction_inverse_laplacian_cube(const RunConfig& run_information, const int index_target, const int index_source, const std::vector<CubePanel>& cube_panels,
                    const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos,
                    const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral);

void cp_interaction_inverse_laplacian_cube(const RunConfig& run_information, const int index_target, const int index_source, const std::vector<CubePanel>& cube_panels,
                    const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos,
                    const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral);

void cc_interaction_inverse_laplacian_cube(const RunConfig& run_information, const int index_target, const int index_source, const std::vector<CubePanel>& cube_panels,
                    const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos,
                    const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral);                    
#endif
