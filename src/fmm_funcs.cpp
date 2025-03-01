#include <vector>
#include <cmath>
#include <iostream>

#include "interp_utils.hpp"
#include "general_utils.hpp"
#include "structs.hpp"

#include "./fmm_interactions/inverse_laplacian_interactions.hpp"

void assign_points_to_source_leafs(const RunConfig& run_information, const std::vector<CubePanel>& source_tree, const std::vector<double>& xcos, 
				 const std::vector<double>& ycos, const std::vector<double>& zcos, std::vector<int>& point_source_leaf) {
	int lb, ub, new_lb, new_ub;
	bool found_leaf, found;

	for (int i = 0; i < run_information.point_count; i++) { // loop through points
		lb = 0, ub = 6;
		found_leaf = false;
		while (!found_leaf) {
			for (int j = lb; j <= ub; j++) { // find the panel at a given level which contains the point
				found = false;
				if (source_tree[j].point_count > 0) {
					for (int k = 0; k < source_tree[j].point_count; k++) {
						if (i == source_tree[j].points_inside[k]) {
							found = true;
							if (source_tree[j].is_leaf) {
								point_source_leaf[i] = j;
								found_leaf = true;
								break;
							} else {
								new_lb = source_tree[j].child_panel_1;
								new_ub = source_tree[j].child_panel_4;
								break;
							}	
						}
					}
				}
				if (found) break;
			}
			if (found_leaf) break;
			lb = new_lb;
			ub = new_ub;
		}
	}
}

void upward_pass(const RunConfig& run_information, const std::vector<CubePanel>& source_tree, const std::vector<double>& xcos, 
				 const std::vector<double>& ycos, const std::vector<double>& zcos, const std::vector<double>& area, 
				 const std::vector<double>& potential, const std::vector<int>& point_source_leaf, 
				 std::vector<std::vector<std::vector<double>>>& proxy_source_weights) {
	// performs the upward pass to interpolate from points to weights in the leaf panels, then from leaf panels up the source tree
	// std::vector<int> point_source_leaf (run_information.point_count, -1);
	std::vector<double> xieta, cheb_xi, cheb_eta;
	std::vector<std::vector<double>> basis_vals;
	double x, y, z, min_xc, max_xc, min_ec, max_ec;
	int leaf, deg, source_size, parent, i_s;

	deg = run_information.interp_degree;
	source_size = source_tree.size();
	
	for (int i = 0; i < run_information.point_count; i++) {
		// compute proxy source weights
		x = xcos[i];
		y = ycos[i];
		z = zcos[i];
		leaf = point_source_leaf[i];
		xieta = xieta_from_xyz(x, y, z, source_tree[leaf].face);
		basis_vals = interp_vals_bli(xieta[0], xieta[1], source_tree[leaf].min_xi, source_tree[leaf].max_xi, source_tree[leaf].min_eta, source_tree[leaf].max_eta, deg);
		for (int j = 0; j < deg+1; j++) {
			for (int k = 0; k < deg+1; k++) {
				proxy_source_weights[leaf][j][k] += basis_vals[j][k]*potential[i]*area[i];
			}
		}
	}

	for (int i = source_size-1; i>-1;i--) {
		// interpolate from child to parent
		parent = source_tree[i].parent_panel;
		if ((parent > -1) and (source_tree[i].point_count > 0)) {
			min_xc = source_tree[i].min_xi;
			max_xc = source_tree[i].max_xi;
			min_ec = source_tree[i].min_eta;
			max_ec = source_tree[i].max_eta;
			cheb_xi = bli_interp_points_shift(min_xc, max_xc, deg);
			cheb_eta = bli_interp_points_shift(min_ec, max_ec, deg);
			min_xc = source_tree[parent].min_xi;
			max_xc = source_tree[parent].max_xi;
			min_ec = source_tree[parent].min_eta;
			max_ec = source_tree[parent].max_eta;
			for (int j = 0; j < deg+1; j++) {
				for (int k = 0; k < deg+1; k++) {
					basis_vals = interp_vals_bli(cheb_xi[j], cheb_eta[k], min_xc, max_xc, min_ec, max_ec, deg);
					for (int l = 0; l < deg+1; l++) {
						for (int m = 0; m < deg+1; m++) {
							proxy_source_weights[parent][l][m] += basis_vals[l][m]*proxy_source_weights[i][j][k];
						}
					}
				}
			}
		}
	}
}

void downward_pass(const RunConfig& run_information, const std::vector<CubePanel>& target_tree, const std::vector<double>& xcos, const std::vector<double>& ycos, 
				   const std::vector<double>& zcos, std::vector<std::vector<std::vector<double>>>& proxy_target_potenti, std::vector<double>& integral) {
	// downward pass to interpolate from parent to child panels
	int target_size, i_t, deg;
	double min_xc, max_xc, min_ec, max_ec, x, y, z, min_xp, max_xp, min_ep, max_ep;
	std::vector<double> xieta, cheb_xc, cheb_ec;
	std::vector<std::vector<double>> basis_vals;

	target_size = target_tree.size();
	deg = run_information.interp_degree;

	for (int i = 0; i < target_size; i++) {
		min_xp = target_tree[i].min_xi;
		max_xp = target_tree[i].max_xi;
		min_ep = target_tree[i].min_eta;
		max_ep = target_tree[i].max_eta;
		if (!target_tree[i].is_leaf) { // interpolate from parent panel to child panel
			for (int j = target_tree[i].child_panel_1; j <= target_tree[i].child_panel_4; j++) {
				min_xc = target_tree[j].min_xi;
				max_xc = target_tree[j].max_xi;
				min_ec = target_tree[j].min_eta;
				max_ec = target_tree[j].max_eta;
				cheb_xc = bli_interp_points_shift(min_xc, max_xc, deg);
				cheb_ec = bli_interp_points_shift(min_ec, max_ec, deg);
				for (int k = 0; k < deg+1; k++) {
					for (int l = 0; l < deg+1; l++) {
						basis_vals = interp_vals_bli(cheb_xc[k], cheb_ec[l], min_xp, max_xp, min_ep, max_ep, deg);
						for (int m = 0; m<deg+1; m++) {
							for (int n = 0; n<deg+1; n++) {
								proxy_target_potenti[j][k][l] += proxy_target_potenti[i][m][n]*basis_vals[m][n];
							}
						}
					}
				}
			}
		} else { // interpolate from leaf panel to points inside
			for (int j = 0; j < target_tree[i].point_count; j++) {
				i_t = target_tree[i].points_inside[j];
				x = xcos[i_t];
				y = ycos[i_t];
				z = zcos[i_t];
				xieta = xieta_from_xyz(x, y, z, target_tree[i].face);
				basis_vals = interp_vals_bli(xieta[0], xieta[1], min_xp, max_xp, min_ep, max_ep, deg);
				for (int k = 0; k < deg+1; k++) {
					for (int l = 0; l < deg+1; l++) {
						integral[i_t] += basis_vals[k][l]*proxy_target_potenti[i][k][l];
					}
				}
			}	
		}
	}
}

void fmm_inverse_laplacian(const RunConfig& run_information, const std::vector<InteractPair>& interactions, const std::vector<CubePanel>& source_tree,
						   const std::vector<CubePanel>& target_tree, const std::vector<double>& xcos, const std::vector<double>& ycos, 
						   const std::vector<double>& zcos, const std::vector<double>& area, const std::vector<double>& potential, 
						   const std::vector<int>& point_source_leaf, std::vector<double>& integral) {
	// use an upward pass and downward pass to compute proxy source weights/proxy target potentials
	// four types of interactions
	int source_tree_count = source_tree.size();
	int target_tree_count = target_tree.size();
	int id = run_information.interp_degree;
	int i_t, i_s;
	std::vector<std::vector<std::vector<double>>> proxy_source_weights (source_tree_count, 
																		std::vector<std::vector<double>> (id+1, 
																										  std::vector<double> (id+1, 0)));
	std::vector<std::vector<std::vector<double>>> proxy_target_potenti (target_tree_count, 
																		std::vector<std::vector<double>> (id+1, 
																										  std::vector<double> (id+1, 0)));

	upward_pass(run_information, source_tree, xcos, ycos, zcos, area, potential, point_source_leaf, proxy_source_weights);

	for (int i = 0; i < interactions.size(); i++) {
		if (i % run_information.mpi_P == run_information.mpi_ID) {
			i_t = interactions[i].index_target;
			i_s = interactions[i].index_source;
			// std::cout << i << " " << i_t << " " << i_s << std::endl;
			if (interactions[i].interact_type == 0) { // PP
				fmm_pp_interaction_inverse_laplacian(run_information, source_tree[i_s], target_tree[i_t], xcos, ycos, zcos, area, potential, integral);
			} else if (interactions[i].interact_type == 1) { // PC
				// fmm_pp_interaction_inverse_laplacian(run_information, source_tree[i_s], target_tree[i_t], xcos, ycos, zcos, area, potential, integral);
				fmm_pc_interaction_inverse_laplacian(run_information, source_tree[i_s], target_tree[i_t], xcos, ycos, zcos, area, proxy_source_weights[i_s], integral);
			} else if (interactions[i].interact_type == 2) { // CP
				// fmm_pp_interaction_inverse_laplacian(run_information, source_tree[i_s], target_tree[i_t], xcos, ycos, zcos, area, potential, integral);
				fmm_cp_interaction_inverse_laplacian(run_information, source_tree[i_s], target_tree[i_t], xcos, ycos, zcos, area, potential, proxy_target_potenti[i_t]);
			} else { // CC
				// fmm_pp_interaction_inverse_laplacian(run_information, source_tree[i_s], target_tree[i_t], xcos, ycos, zcos, area, potential, integral);
				fmm_cc_interaction_inverse_laplacian(run_information, source_tree[i_s], target_tree[i_t], xcos, ycos, zcos, area, proxy_source_weights[i_s], proxy_target_potenti[i_t]);
				// fmm_pc_interaction_inverse_laplacian(run_information, source_tree[i_s], target_tree[i_t], xcos, ycos, zcos, area, proxy_source_weights[i_s], integral);
			}
		}
	}

	downward_pass(run_information, target_tree, xcos, ycos, zcos, proxy_target_potenti, integral);
}