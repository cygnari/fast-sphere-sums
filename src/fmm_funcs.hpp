#ifndef H_FMM_FUNCS_H
#define H_FMM_FUNCS_H

#include <vector>

void fmm_inverse_laplacian(const RunConfig& run_information, const std::vector<InteractPair>& interactions, const std::vector<CubePanel>& source_tree,
						   const std::vector<CubePanel>& target_tree, const std::vector<double>& xcos, const std::vector<double>& ycos, 
						   const std::vector<double>& zcos, const std::vector<double>& area, const std::vector<double>& potential, 
						   const std::vector<int>& point_source_leaf, std::vector<double>& integral);

void fmm_bve(const RunConfig& run_information, const std::vector<InteractPair>& interactions, const std::vector<CubePanel>& source_tree,
						   const std::vector<CubePanel>& target_tree, const std::vector<double>& xcos, const std::vector<double>& ycos, 
						   const std::vector<double>& zcos, const std::vector<double>& area, const std::vector<double>& potential, 
						   const std::vector<int>& point_source_leaf, std::vector<double>& integral1, std::vector<double>& integral2, 
						   std::vector<double>& integral3);

#endif