#ifndef H_FMM_FUNCS_H
#define H_FMM_FUNCS_H

#include <vector>

void assign_points_to_source_leafs(const RunConfig& run_information, const std::vector<CubePanel>& source_tree, const std::vector<double>& xcos, 
				 const std::vector<double>& ycos, const std::vector<double>& zcos, std::vector<int>& point_source_leaf);

void fmm_inverse_laplacian(const RunConfig& run_information, const std::vector<InteractPair>& interactions, const std::vector<CubePanel>& source_tree,
						   const std::vector<CubePanel>& target_tree, const std::vector<double>& xcos, const std::vector<double>& ycos, 
						   const std::vector<double>& zcos, const std::vector<double>& area, const std::vector<double>& potential, 
						   const std::vector<int>& point_source_leaf, std::vector<double>& integral);

#endif