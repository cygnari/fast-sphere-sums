#ifndef H_TREE_FUNCS_H
#define H_TREE_FUNCS_H

#include <vector>
#include "structs.hpp"

void initialize_icosahedron_tree(const RunConfig& run_information, std::vector<IcosPanel>& icos_panels, const std::vector<double>& x,
        const std::vector<double>& y, const std::vector<double>& z);

void initialize_cube_tree(const RunConfig& run_information, std::vector<CubePanel>& cube_panels, const std::vector<double>& x,
        const std::vector<double>& y, const std::vector<double>& z);

void dual_tree_traversal_icos(const RunConfig& run_information, std::vector<InteractPair>& interactions, const std::vector<IcosPanel>& icos_panels);

void dual_tree_traversal_cube(const RunConfig& run_information, std::vector<InteractPair>& interactions, const std::vector<CubePanel>& cube_panels);

#endif
