#include <iostream>
#include <mpi.h>
#include <vector>
#include <cmath>
#include <chrono>

#include "fast-sphere-sums-config.h"
#include "direct_sum_funcs.hpp"
#include "fmm_funcs.hpp"
#include "general_utils.hpp"
#include "initial_conditions.hpp"
#include "initialize_tree.hpp"
#include "interp_utils.hpp"
#include "io_utils.hpp"
#include "mpi_utils.hpp"
#include "structs.hpp"

int check_point_exist(const std::vector<std::vector<int>> &parent_points,
                      const int point_count, const int iv1, const int iv2) {
  for (int i = 0; i < point_count; i++) {
	if ((parent_points[i][0] == iv1) and (parent_points[i][1] == iv2)) return i;
  }
  return -1;
}

void dynamics_points_initialize(RunConfig& run_information, std::vector<std::vector<double>>& coordinates, std::vector<double>& areas,
								std::vector<std::vector<std::vector<int>>>& dynamics_triangles, std::vector<std::vector<bool>>& dynamics_triangles_is_leaf) {
	int levels = run_information.levels;
	dynamics_triangles.clear();
	dynamics_triangles_is_leaf.clear();
	coordinates.clear();
	areas.clear();
	int point_count = 10 * pow(4, levels-1) + 2;
	run_information.point_count = point_count;
	int tri_count = 2 * (point_count - 2);
	std::vector<double> xco (point_count, 0);
	std::vector<double> yco (point_count, 0);
	std::vector<double> zco (point_count, 0);
	// std::vector<double> areas (point_count, 0);
	double phi = (1 + sqrt(5)) / 2;
	double point_norm = sqrt(1 + phi*phi);
	double ipn = 1.0/point_norm;
	double phipn = phi/point_norm;

	// initial twelve points
	xco[0] = 0, yco[0] = ipn, zco[0] = phipn;
	xco[1] = 0, yco[1] = -ipn, zco[1] = phipn;
	xco[2] = 0, yco[2] = ipn, zco[2] = -phipn;
	xco[3] = 0, yco[3] = -ipn, zco[3] = -phipn;
	xco[4] = ipn, yco[4] = phipn, zco[4] = 0;
	xco[5] = ipn, yco[5] = -phipn, zco[5] = 0;
	xco[6] = -ipn, yco[6] = phipn, zco[6] = 0;
	xco[7] = -ipn, yco[7] = -phipn, zco[7] = 0;
	xco[8] = phipn, yco[8] = 0, zco[8] = ipn;
	xco[9] = phipn, yco[9] = 0, zco[9] = -ipn;
	xco[10] = -phipn, yco[10] = 0, zco[10] = ipn;
	xco[11] = -phipn, yco[11] = 0, zco[11] = -ipn;

	std::vector<std::vector<int>> dynamics_points_parents;

	dynamics_triangles.resize(levels);
	dynamics_triangles[0] = std::vector<std::vector<int>>(20, std::vector<int>(4, 0));
	dynamics_triangles_is_leaf.resize(levels);
	dynamics_points_parents.resize(point_count, std::vector<int>(2, -1));
	areas.resize(point_count);
	coordinates.resize(3);

	// 0, 1, 2 are indices of the three vertices
	// 20 starting faces
	// make sure the points are in CCW order
	dynamics_triangles[0][0].insert(dynamics_triangles[0][0].begin(), {0, 1, 8, 0});
	dynamics_triangles[0][1].insert(dynamics_triangles[0][1].begin(), {0, 10, 1, 0});
	dynamics_triangles[0][2].insert(dynamics_triangles[0][2].begin(), {0, 4, 6, 0});
	dynamics_triangles[0][3].insert(dynamics_triangles[0][3].begin(), {0, 8, 4, 0});
	dynamics_triangles[0][4].insert(dynamics_triangles[0][4].begin(), {0, 6, 10, 0});
	dynamics_triangles[0][5].insert(dynamics_triangles[0][5].begin(), {1, 7, 5, 0});
	dynamics_triangles[0][6].insert(dynamics_triangles[0][6].begin(), {1, 5, 8, 0});
	dynamics_triangles[0][7].insert(dynamics_triangles[0][7].begin(), {1, 10, 7, 0});
	dynamics_triangles[0][8].insert(dynamics_triangles[0][8].begin(), {2, 9, 3, 0});
	dynamics_triangles[0][9].insert(dynamics_triangles[0][9].begin(), {2, 3, 11, 0});
	dynamics_triangles[0][10].insert(dynamics_triangles[0][10].begin(), {2, 6, 4, 0});
	dynamics_triangles[0][11].insert(dynamics_triangles[0][11].begin(), {2, 4, 9, 0});
	dynamics_triangles[0][12].insert(dynamics_triangles[0][12].begin(), {2, 11, 6, 0});
	dynamics_triangles[0][13].insert(dynamics_triangles[0][13].begin(), {3, 5, 7, 0});
	dynamics_triangles[0][14].insert(dynamics_triangles[0][14].begin(), {3, 9, 5, 0});
	dynamics_triangles[0][15].insert(dynamics_triangles[0][15].begin(), {3, 7, 11, 0});
	dynamics_triangles[0][16].insert(dynamics_triangles[0][16].begin(), {4, 8, 9, 0});
	dynamics_triangles[0][17].insert(dynamics_triangles[0][17].begin(), {5, 9, 8, 0});
	dynamics_triangles[0][18].insert(dynamics_triangles[0][18].begin(), {6, 11, 10, 0});
	dynamics_triangles[0][19].insert(dynamics_triangles[0][19].begin(), {7, 10, 11, 0});

	int curr_point_count = 12;

	int iv1, iv2, iv3, iv12, iv23, iv31;
	std::vector<double> v1 (3, 0), v2 (3, 0), v3 (3, 0), v12 (3, 0), v23 (3, 0), v31 (3, 0);

	for (int i = 0; i < levels-1; i++) {
		// iteratively refine
		dynamics_triangles[i + 1] = std::vector<std::vector<int>>(20 * pow(4, i + 1), std::vector<int>(4, 0));
		for (int j = 0; j < 20 * pow(4, i); j++) {
			iv1 = dynamics_triangles[i][j][0];
			iv2 = dynamics_triangles[i][j][1];
			iv3 = dynamics_triangles[i][j][2];
			v1[0] = xco[iv1], v1[1] = yco[iv1], v1[2] = zco[iv1];
			v2[0] = xco[iv2], v2[1] = yco[iv2], v2[2] = zco[iv2];
			v3[0] = xco[iv3], v3[1] = yco[iv3], v3[2] = zco[iv3];
			v12[0] = 0.5*(xco[iv1]+xco[iv2]), v12[1] = 0.5*(yco[iv1]+yco[iv2]), v12[2] = 0.5*(zco[iv1]+zco[iv2]);
			v23[0] = 0.5*(xco[iv2]+xco[iv3]), v23[1] = 0.5*(yco[iv2]+yco[iv3]), v23[2] = 0.5*(zco[iv2]+zco[iv3]);
			v31[0] = 0.5*(xco[iv3]+xco[iv1]), v31[1] = 0.5*(yco[iv3]+yco[iv1]), v31[2] = 0.5*(zco[iv3]+zco[iv1]);
			project_to_sphere(v12, 1);
			project_to_sphere(v23, 1);
			project_to_sphere(v31, 1);
			iv12 = check_point_exist(dynamics_points_parents, curr_point_count, std::min(iv1, iv2), std::max(iv1, iv2));
			iv23 = check_point_exist(dynamics_points_parents, curr_point_count, std::min(iv2, iv3), std::max(iv2, iv3));
			iv31 = check_point_exist(dynamics_points_parents, curr_point_count, std::min(iv3, iv1), std::max(iv3, iv1));
			if (iv12 == -1) {
				iv12 = curr_point_count;
				curr_point_count += 1;
				dynamics_points_parents[iv12] = {std::min(iv1, iv2), std::max(iv1, iv2)};
				xco[iv12] = v12[0], yco[iv12] = v12[1], zco[iv12] = v12[2];
			}
			if (iv23 == -1) {
				iv23 = curr_point_count;
				curr_point_count += 1;
				dynamics_points_parents[iv23] = {std::min(iv2, iv3), std::max(iv2, iv3)};
				xco[iv23] = v23[0], yco[iv23] = v23[1], zco[iv23] = v23[2];
			}
			if (iv31 == -1) {
				iv31 = curr_point_count;
				curr_point_count += 1;
				dynamics_points_parents[iv31] = {std::min(iv3, iv1), std::max(iv3, iv1)};
				xco[iv31] = v31[0], yco[iv31] = v31[1], zco[iv31] = v31[2];
			}
			dynamics_triangles[i + 1][4 * j].insert(dynamics_triangles[i + 1][4 * j].begin(), {iv1, iv12, iv31, i + 1});
			dynamics_triangles[i + 1][4 * j + 1].insert(dynamics_triangles[i + 1][4 * j + 1].begin(), {iv2, iv23, iv12, i + 1});
			dynamics_triangles[i + 1][4 * j + 2].insert(dynamics_triangles[i + 1][4 * j + 2].begin(), {iv3, iv31, iv23, i + 1});
			dynamics_triangles[i + 1][4 * j + 3].insert(dynamics_triangles[i + 1][4 * j + 3].begin(), {iv12, iv23, iv31, i + 1});
		}
	}

	double tri_area;
	for (int i = 0; i < tri_count; i++) {
		// compute area associated to each point
		iv1 = dynamics_triangles[levels-1][i][0];
		iv2 = dynamics_triangles[levels-1][i][1];
		iv3 = dynamics_triangles[levels-1][i][2];
		v1[0] = xco[iv1], v1[1] = yco[iv1], v1[2] = zco[iv1];
		v2[0] = xco[iv2], v2[1] = yco[iv2], v2[2] = zco[iv2];
		v3[0] = xco[iv3], v3[1] = yco[iv3], v3[2] = zco[iv3];
		tri_area = sphere_tri_area(v1, v2, v3, 1.0);
		areas[iv1] += tri_area / 3.0;
		areas[iv2] += tri_area / 3.0;
		areas[iv3] += tri_area / 3.0;
	}

	for (int i = 0; i < levels; i++) {
		if (i == levels - 1) {
			dynamics_triangles_is_leaf[i].resize(20 * pow(4, i), true);
		} else {
			dynamics_triangles_is_leaf[i].resize(20 * pow(4, i), false);
		}
	}
	coordinates[0] = xco;
	coordinates[1] = yco;
	coordinates[2] = zco;
}

void remesh_points(const RunConfig &run_information, const std::vector<std::vector<double>> &target_cos, const std::vector<std::vector<double>> &new_co, std::vector<double>& target_vors, 
    const std::vector<double>& new_vors, std::vector<double>& target_tra, const std::vector<double>& new_tra, const std::vector<std::vector<std::vector<int>>> &dynamics_triangles, 
    const std::vector<std::vector<bool>> &dynamics_triangles_is_leaf, const int point_count, const double omega) {
	int curr_level, tri_loc, super_tri_loc;
	std::vector<double> v1, v2, v3, v4, v5, v6;
	double vor1, vor2, vor3, vor4, vor5, vor6, vormax, vormin, vor, tra, vormed, vors, tra1, tra2, tra3;
	std::vector<double> interp_matrix(36, 0);
	std::vector<int> ipiv(6, 0);
	std::vector<double> b_vec(12, 0);
	std::vector<std::vector<double>> points(6, std::vector<double>(3, 0));
	std::vector<double> ivs (6, 0);
	std::vector<double> bary_cords, curr_alphas;

	std::fill(target_vors.begin(), target_vors.end(), 0);
	std::fill(target_tra.begin(), target_tra.end(), 0);

	for (int i = run_information.point_lb; i < run_information.point_ub; i++) {
		std::tie(curr_level, tri_loc) = find_leaf_tri(target_cos[0][i], target_cos[1][i], target_cos[2][i], new_co[0], new_co[1], new_co[2],
                 	dynamics_triangles, dynamics_triangles_is_leaf, run_information.levels);
		super_tri_loc = floor(tri_loc / 4.0);

		ivs[0] = dynamics_triangles[curr_level - 1][super_tri_loc][0];
		ivs[1] = dynamics_triangles[curr_level - 1][super_tri_loc][1];
		ivs[2] = dynamics_triangles[curr_level - 1][super_tri_loc][2];
		ivs[3] = dynamics_triangles[curr_level][4 * super_tri_loc + 3][0];
		ivs[4] = dynamics_triangles[curr_level][4 * super_tri_loc + 3][1];
		ivs[5] = dynamics_triangles[curr_level][4 * super_tri_loc + 3][2];
		v1 = {new_co[0][ivs[0]], new_co[1][ivs[0]], new_co[2][ivs[0]]};
		v2 = {new_co[0][ivs[1]], new_co[1][ivs[1]], new_co[2][ivs[1]]};
		v3 = {new_co[0][ivs[2]], new_co[1][ivs[2]], new_co[2][ivs[2]]};
		v4 = {new_co[0][ivs[3]], new_co[1][ivs[3]], new_co[2][ivs[3]]};
		v5 = {new_co[0][ivs[4]], new_co[1][ivs[4]], new_co[2][ivs[4]]};
		v6 = {new_co[0][ivs[5]], new_co[1][ivs[5]], new_co[2][ivs[5]]};

		points[0][0] = 1;
		points[1][1] = 1;
		points[2][2] = 1;
		points[3] = normalized_barycoords(v1, v2, v3, v4);
		points[4] = normalized_barycoords(v1, v2, v3, v5);
		points[5] = normalized_barycoords(v1, v2, v3, v6);

		for (int j = 0; j < 6; j++) {
			b_vec[j] = new_vors[ivs[j]];
			b_vec[j+6] = new_tra[ivs[j]];
		}

		interp_mat_init_sbb(interp_matrix, points, 2, 6);

		int info = linear_solve(interp_matrix, b_vec, 6, 2, 1);
		if (info > 0) throw std::runtime_error("Biquadratic interpolation linear solve failed at line 136");
		bary_cords = barycoords(v1, v2, v3, target_cos[0][i], target_cos[1][i], target_cos[2][i]);
		curr_alphas = {b_vec[0], b_vec[1], b_vec[2], b_vec[3], b_vec[4], b_vec[5]};
		// target_vors[i]
		vor = interp_eval_sbb(curr_alphas, bary_cords[0], bary_cords[1], bary_cords[2], 2);
		curr_alphas = {b_vec[6], b_vec[7], b_vec[8], b_vec[9], b_vec[10], b_vec[11]};
		// target_tra[i] 
		tra = interp_eval_sbb(curr_alphas, bary_cords[0], bary_cords[1], bary_cords[2], 2);
		// check monotonicity 
		// vor1 = new_vors[ivs[0]]+2*omega*new_co[2][ivs[0]];
		// vor2 = new_vors[ivs[1]]+2*omega*new_co[2][ivs[1]];
		// vor3 = new_vors[ivs[2]]+2*omega*new_co[2][ivs[2]];
		// vor4 = new_vors[ivs[3]]+2*omega*new_co[2][ivs[3]];
		// vor5 = new_vors[ivs[4]]+2*omega*new_co[2][ivs[4]];
		// vor6 = new_vors[ivs[5]]+2*omega*new_co[2][ivs[5]];
		// vormax = std::max(vor1, std::max(vor2, std::max(vor3, std::max(vor4, std::max(vor5, vor6)))));
		// vormin = std::min(vor1, std::min(vor2, std::min(vor3, std::min(vor4, std::min(vor5, vor6)))));
		// vormed = 0.5*(vormax+vormin);
		// vors = 0.5*(vormax-vormin);
		// vormax = vormed+1.5*vors; // a bit of leeway
		// vormin = vormed-1.5*vors;
		// if ((vor > vormax) or (vor < vormin)) {
		// 	// monotonicity violation, bilinear interp
		// 	ivs[0] = dynamics_triangles[curr_level][tri_loc][0];
		// 	ivs[1] = dynamics_triangles[curr_level][tri_loc][1];
		// 	ivs[2] = dynamics_triangles[curr_level][tri_loc][2];
		// 	v1 = {new_co[0][ivs[0]], new_co[1][ivs[0]], new_co[2][ivs[0]]};
		// 	v2 = {new_co[0][ivs[1]], new_co[1][ivs[1]], new_co[2][ivs[1]]};
		// 	v3 = {new_co[0][ivs[2]], new_co[1][ivs[2]], new_co[2][ivs[2]]};
		// 	vor1 = new_vors[ivs[0]]+2*omega*new_co[2][ivs[0]];
		// 	vor2 = new_vors[ivs[1]]+2*omega*new_co[2][ivs[1]];
		// 	vor3 = new_vors[ivs[2]]+2*omega*new_co[2][ivs[2]];
		// 	tra1 = new_tra[ivs[0]];
		// 	tra2 = new_tra[ivs[1]];
		// 	tra3 = new_tra[ivs[2]];
		// 	bary_cords = barycoords(v1, v2, v3, target_cos[0][i], target_cos[1][i], target_cos[2][i]);
		// 	vor = bary_cords[0]*vor1 + bary_cords[1]*vor2 + bary_cords[2]*vor3 - 2*omega*target_cos[2][i];
		// 	tra = bary_cords[0]*tra1 + bary_cords[1]*tra2 + bary_cords[2]*tra3;
		// }
		target_vors[i] = vor;
		target_tra[i] = tra;
	}
}

int main(int argc, char **argv) {
	MPI_Init(&argc, &argv);
	int P, ID;
	MPI_Status status;
	MPI_Comm_size(MPI_COMM_WORLD, &P);
	MPI_Comm_rank(MPI_COMM_WORLD, &ID);
	MPI_Win win_velxt, win_velyt, win_velzt, win_vor, win_tra;

	std::chrono::steady_clock::time_point begin, end;

	RunConfig run_information;
	const std::string namelist_file = std::string(NAMELIST_DIR) + std::string("namelist.txt");
	read_run_config(namelist_file, run_information);
	run_information.mpi_P = P;
	run_information.mpi_ID = ID;

	double omega = 7.292115855e-5; // 2 pi per sidereal day in 1/seconds

	std::vector<std::vector<double>> coordinates;
	std::vector<double> area;
	std::vector<std::vector<std::vector<int>>> dynamics_triangles;
	std::vector<std::vector<bool>> dynamics_triangles_is_leaf;

	dynamics_points_initialize(run_information, coordinates, area, dynamics_triangles, dynamics_triangles_is_leaf);

	std::vector<double> vorticity (run_information.point_count, 0);
	std::vector<double> velx (run_information.point_count, 0);
	std::vector<double> vely (run_information.point_count, 0);
	std::vector<double> velz (run_information.point_count, 0);
	std::vector<double> passive_tracer (run_information.point_count, 0);

	std::vector<std::vector<double>> inter_state (3, std::vector<double> (run_information.point_count, 0));
	std::vector<double> inter_vor (run_information.point_count, 0);
	std::vector<double> inter_tra (run_information.point_count, 0);
	std::vector<double> velxt (run_information.point_count, 0), velyt (run_information.point_count, 0), velzt (run_information.point_count, 0);

	initialize_condition(run_information, coordinates[0], coordinates[1], coordinates[2], vorticity);
	for (int i = 0; i < vorticity.size(); i++) {
		vorticity[i] /= 86400; // convert vorticity [1/day] to vorticity [1/second]
		passive_tracer[i] = coordinates[2][i];
	}

	MPI_Win_create(&velxt[0], run_information.point_count * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_velxt);
	MPI_Win_create(&velyt[0], run_information.point_count * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_velyt);
	MPI_Win_create(&velzt[0], run_information.point_count * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_velzt);
	MPI_Win_create(&vorticity[0], run_information.point_count * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_vor);
	MPI_Win_create(&passive_tracer[0], run_information.point_count * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_tra);

	std::vector<CubePanel> cube_panels;
	std::vector<int> point_source_leaf (run_information.point_count, -1);
	std::vector<InteractPair> interactions;

	if (run_information.use_fast) {
		initialize_cube_tree(run_information, cube_panels, coordinates[0], coordinates[1], coordinates[2], point_source_leaf);
		dual_tree_traversal(run_information, interactions, cube_panels);
	} else {
		bounds_determine_2d(run_information, P, ID);
	}
	bounds_determine_1d(run_information, P, ID);

	std::string output_folder = create_config(run_information) + "_bve_evo";
	std::string outpath = run_information.out_path + "/" + output_folder + "/";
	if (ID == 0) {
		std::string filename = NAMELIST_DIR + std::string("initialize.py ") + run_information.out_path + "/" + output_folder;
		std::string command = "python ";
		command += filename;
		system(command.c_str());
		write_state(vorticity, outpath, "vor_init.csv");
		write_state(coordinates[0], outpath, "x_init.csv");
		write_state(coordinates[1], outpath, "y_init.csv");
		write_state(coordinates[2], outpath, "z_init.csv");
		write_state(area, outpath, "areas_init.csv");
		write_state(passive_tracer, outpath, "tracer_init.csv");
	}

	int time;
	double ix, iy, iz, norm;
	
	for (int i = 0; i < run_information.time_steps; i++) {
		inter_state = coordinates;
		time = (i+1) * run_information.delta_t;
		if (ID == 0) std::cout << time << std::endl;
		if (run_information.use_fast) { // RK4 step 1
			fmm_bve(run_information, interactions, cube_panels, cube_panels, coordinates[0], coordinates[1], coordinates[2], area, vorticity, point_source_leaf, velxt, velyt, velzt);
		} else {
			direct_sum_bve(run_information, coordinates[0], coordinates[1], coordinates[2], area, vorticity, velxt, velyt, velzt);
		}
		sync_updates<double>(velxt, P, ID, &win_velxt, MPI_DOUBLE);
		sync_updates<double>(velyt, P, ID, &win_velyt, MPI_DOUBLE);
		sync_updates<double>(velzt, P, ID, &win_velzt, MPI_DOUBLE);
		for (int j = 0; j < run_information.point_count; j++) {
			ix = coordinates[0][j] + velxt[j]*run_information.delta_t/2.0;
			iy = coordinates[1][j] + velyt[j]*run_information.delta_t/2.0;
			iz = coordinates[2][j] + velzt[j]*run_information.delta_t/2.0;
			norm = sqrt(ix*ix+iy*iy+iz*iz);
			inter_state[0][j] = ix/norm;
			inter_state[1][j] = iy/norm;
			inter_state[2][j] = iz/norm;
			inter_vor[j] = vorticity[j]+2*omega*(coordinates[2][j]-iz/norm);
			velx[j] = velxt[j];
			vely[j] = velyt[j];
			velz[j] = velzt[j];
		}
		if (run_information.use_fast) { // RK4 step 2
			fmm_bve(run_information, interactions, cube_panels, cube_panels, inter_state[0], inter_state[1], inter_state[2], area, inter_vor, point_source_leaf, velxt, velyt, velzt);
		} else {
			direct_sum_bve(run_information, inter_state[0], inter_state[1], inter_state[2], area, inter_vor, velxt, velyt, velzt);
		}
		sync_updates<double>(velxt, P, ID, &win_velxt, MPI_DOUBLE);
		sync_updates<double>(velyt, P, ID, &win_velyt, MPI_DOUBLE);
		sync_updates<double>(velzt, P, ID, &win_velzt, MPI_DOUBLE);
		for (int j = 0; j < run_information.point_count; j++) {
			ix = coordinates[0][j] + velxt[j]*run_information.delta_t/2.0;
			iy = coordinates[1][j] + velyt[j]*run_information.delta_t/2.0;
			iz = coordinates[2][j] + velzt[j]*run_information.delta_t/2.0;
			norm = sqrt(ix*ix+iy*iy+iz*iz);
			inter_state[0][j] = ix/norm;
			inter_state[1][j] = iy/norm;
			inter_state[2][j] = iz/norm;
			inter_vor[j] = vorticity[j]+2*omega*(coordinates[2][j]-iz/norm);
			velx[j] += 2*velxt[j];
			vely[j] += 2*velyt[j];
			velz[j] += 2*velzt[j];
		}
		if (run_information.use_fast) { // RK4 step 3
			fmm_bve(run_information, interactions, cube_panels, cube_panels, inter_state[0], inter_state[1], inter_state[2], area, inter_vor, point_source_leaf, velxt, velyt, velzt);
		} else {
			direct_sum_bve(run_information, inter_state[0], inter_state[1], inter_state[2], area, inter_vor, velxt, velyt, velzt);
		}
		sync_updates<double>(velxt, P, ID, &win_velxt, MPI_DOUBLE);
		sync_updates<double>(velyt, P, ID, &win_velyt, MPI_DOUBLE);
		sync_updates<double>(velzt, P, ID, &win_velzt, MPI_DOUBLE);
		for (int j = 0; j < run_information.point_count; j++) {
			ix = coordinates[0][j] + velxt[j]*run_information.delta_t;
			iy = coordinates[1][j] + velyt[j]*run_information.delta_t;
			iz = coordinates[2][j] + velzt[j]*run_information.delta_t;
			norm = sqrt(ix*ix+iy*iy+iz*iz);
			inter_state[0][j] = ix/norm;
			inter_state[1][j] = iy/norm;
			inter_state[2][j] = iz/norm;
			inter_vor[j] = vorticity[j]+2*omega*(coordinates[2][j]-iz/norm);
			velx[j] += 2*velxt[j];
			vely[j] += 2*velyt[j];
			velz[j] += 2*velzt[j];
		}
		if (run_information.use_fast) { // RK4 step 4
			fmm_bve(run_information, interactions, cube_panels, cube_panels, inter_state[0], inter_state[1], inter_state[2], area, inter_vor, point_source_leaf, velxt, velyt, velzt);
		} else {
			direct_sum_bve(run_information, inter_state[0], inter_state[1], inter_state[2], area, inter_vor, velxt, velyt, velzt);
		}
		sync_updates<double>(velxt, P, ID, &win_velxt, MPI_DOUBLE);
		sync_updates<double>(velyt, P, ID, &win_velyt, MPI_DOUBLE);
		sync_updates<double>(velzt, P, ID, &win_velzt, MPI_DOUBLE);
		for (int j = 0; j < run_information.point_count; j++) { // RK4 update
			velx[j] += velxt[j];
			vely[j] += velyt[j];
			velz[j] += velzt[j];
			velx[j] *= run_information.delta_t / 6.0;
			vely[j] *= run_information.delta_t / 6.0;
			velz[j] *= run_information.delta_t / 6.0;
			ix = coordinates[0][j] + velx[j];
			iy = coordinates[1][j] + vely[j];
			iz = coordinates[2][j] + velz[j];
			norm = sqrt(ix*ix + iy*iy + iz*iz);
			inter_vor[j] = vorticity[j] + 2*omega*(coordinates[2][j]-iz/norm);
			inter_state[0][j] = ix/norm;
			inter_state[1][j] = iy/norm;
			inter_state[2][j] = iz/norm;
			inter_tra[j] = passive_tracer[j];
		}
		if (run_information.use_remesh) {
			remesh_points(run_information, coordinates, inter_state, vorticity, inter_vor, passive_tracer, inter_tra, dynamics_triangles, dynamics_triangles_is_leaf, run_information.point_count, omega);
			sync_updates<double>(vorticity, P, ID, &win_vor, MPI_DOUBLE);
			sync_updates<double>(passive_tracer, P, ID, &win_tra, MPI_DOUBLE);
		} else {
			for (int j = 0; j < run_information.point_count; j++) {
				coordinates[0][j] = inter_state[0][j];
				coordinates[1][j] = inter_state[1][j];
				coordinates[2][j] = inter_state[2][j];
				vorticity[j] = inter_vor[j];
			}
			if (run_information.use_fast) {
				cube_panels.clear();
				point_source_leaf.clear();
				interactions.clear();
				initialize_cube_tree(run_information, cube_panels, coordinates[0], coordinates[1], coordinates[2], point_source_leaf);
				dual_tree_traversal(run_information, interactions, cube_panels);
			}	
		}
		if (ID == 0) {
			write_state(vorticity, outpath, "vor_" + std::to_string(time) + ".csv");
			write_state(coordinates[0], outpath, "x_" + std::to_string(time) + ".csv");
			write_state(coordinates[1], outpath, "y_" + std::to_string(time) + ".csv");
			write_state(coordinates[2], outpath, "z_" + std::to_string(time) + ".csv");
			write_state(passive_tracer, outpath, "tracer_" + std::to_string(time) + ".csv");
		}
	}

	MPI_Finalize();
	return 0;
}