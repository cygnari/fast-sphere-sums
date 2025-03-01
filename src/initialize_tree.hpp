#ifndef H_TREE_FUNCS_H
#define H_TREE_FUNCS_H

#include <vector>
#include "structs.hpp"

void initialize_cube_tree(const RunConfig& run_information, std::vector<CubePanel>& cube_panels, const std::vector<double>& x,
        const std::vector<double>& y, const std::vector<double>& z, std::vector<int>& point_source_leaf);

void dual_tree_traversal(const RunConfig& run_information, std::vector<InteractPair>& interactions, const std::vector<CubePanel>& cube_panels);

void dual_tree_traversal(const RunConfig& run_information, std::vector<InteractPair>& interactions, const std::vector<CubePanel>& cube_panels_target, const std::vector<CubePanel>& cube_panels_source);

#endif
