#ifndef H_FAST_SUM_FUNCS_H
#define H_FAST_SUM_FUNCS_H

#include <vector>

#include "structs.hpp"

void fast_sum_inverse_laplacian(const RunConfig& run_information, const std::vector<InteractPair>& interactions, const std::vector<CubePanel>& cube_panels,
                                        const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos,
                                        const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral);


void fast_sum_inverse_biharmonic(const RunConfig& run_information, const std::vector<InteractPair>& interactions, const std::vector<CubePanel>& cube_panels_target, const std::vector<CubePanel>& cube_panels_source,
                                        const std::vector<double>& xcos_t, const std::vector<double>& ycos_t, const std::vector<double>& zcos_t,
                                        const std::vector<double>& xcos_s, const std::vector<double>& ycos_s, const std::vector<double>& zcos_s,
                                        const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral);
#endif
