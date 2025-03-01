#include <iostream>
#include <mpi.h>
#include <vector>
#include <cmath>
#include <chrono>

#include "fast-sphere-sums-config.h"
#include "direct_sum_funcs.hpp"
#include "fast_sum_funcs.hpp"
#include "general_utils.hpp"
#include "initial_conditions.hpp"
#include "initialize_tree.hpp"
#include "io_utils.hpp"
#include "mpi_utils.hpp"
#include "structs.hpp"

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  int P, ID;
  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD, &P);
  MPI_Comm_rank(MPI_COMM_WORLD, &ID);
  MPI_Win win_integrated;

  MPI_Datatype dt_interaction;
  MPI_Type_contiguous(3, MPI_INT, &dt_interaction);
  MPI_Type_commit(&dt_interaction);

  std::chrono::steady_clock::time_point begin, end;

  RunConfig run_information;
  const std::string namelist_file = std::string(NAMELIST_DIR) + std::string("namelist.txt");
  read_run_config(namelist_file, run_information);
  run_information.mpi_P = P;
  run_information.mpi_ID = ID;

  std::vector<double> xcos_t (run_information.point_count, 0);
  std::vector<double> ycos_t (run_information.point_count, 0);
  std::vector<double> zcos_t (run_information.point_count, 0);
  std::vector<double> xcos_s (run_information.point_count, 0);
  std::vector<double> ycos_s (run_information.point_count, 0);
  std::vector<double> zcos_s (run_information.point_count, 0);
  std::vector<double> area (run_information.point_count, 0);
  std::vector<double> potential (run_information.point_count, 0);
  std::vector<double> integrated (run_information.point_count, 0);

  MPI_Win_create(&integrated[0], run_information.point_count * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_integrated);

  std::string data_pre = DATA_DIR + std::to_string(run_information.point_count) + "_" + run_information.grid + "_";

  read_data_field(run_information.point_count, xcos_s, data_pre + "x.csv");
  read_data_field(run_information.point_count, ycos_s, data_pre + "y.csv");
  read_data_field(run_information.point_count, zcos_s, data_pre + "z.csv");
  read_data_field(run_information.point_count, area, data_pre + "areas.csv");

  xcos_t = xcos_s;
  ycos_t = ycos_s;
  zcos_t = zcos_s;

  if (run_information.rotate) {
    rotate_points(xcos_t, ycos_t, zcos_t, run_information.alph, run_information.beta, run_information.gamm);
    rotate_points(xcos_s, ycos_s, zcos_s, -run_information.alph, -run_information.beta, -run_information.gamm);
  }

  initialize_condition(run_information, xcos_s, ycos_s, zcos_s, potential);
  if (run_information.balance_condition) {
    balance_conditions(potential, area);
  }

  if (run_information.use_fast) {
    // use fast summation
    std::vector<CubePanel> cube_panels_target, cube_panels_source;
    std::vector<int> point_source_leaf (run_information.point_count, -1);
    initialize_cube_tree(run_information, cube_panels_target, xcos_t, ycos_t, zcos_t, point_source_leaf);
    initialize_cube_tree(run_information, cube_panels_source, xcos_s, ycos_s, zcos_s, point_source_leaf);

    std::vector<InteractPair> interactions;
    begin = std::chrono::steady_clock::now();
    dual_tree_traversal(run_information, interactions, cube_panels_target, cube_panels_source);
    end = std::chrono::steady_clock::now();
    if (ID == 0) {
      std::cout << "tree traversal time: " << std::chrono::duration<double>(end - begin).count() << " seconds" << std::endl;
      // std::cout << "interactions: " << interactions.size() << std::endl;
    }
    begin = std::chrono::steady_clock::now();
    fast_sum_inverse_biharmonic(run_information, interactions, cube_panels_source, cube_panels_target, xcos_t, ycos_t, zcos_t, xcos_s, ycos_s, zcos_s, area, potential, integrated);
    sync_updates<double>(integrated, P, ID, &win_integrated, MPI_DOUBLE);
    end = std::chrono::steady_clock::now();
  } else {
    // direct summation
    bounds_determine_2d(run_information, P, ID);
    begin = std::chrono::steady_clock::now();
    direct_sum_inverse_biharmonic(run_information, xcos_t, ycos_t, zcos_t, xcos_s, ycos_s, zcos_s, area, potential, integrated);
    sync_updates<double>(integrated, P, ID, &win_integrated, MPI_DOUBLE);
    end = std::chrono::steady_clock::now();
  }

  if (ID == 0) {
    std::cout << "integration time: " << std::chrono::duration<double>(end - begin).count() << " seconds" << std::endl;

    if (run_information.write_output) {
      std::string output_folder = create_config(run_information) + "_inv_biharm";
      std::string filename = NAMELIST_DIR + std::string("initialize.py ") + run_information.out_path + "/" + output_folder;
      std::string command = "python ";
      command += filename;
      system(command.c_str());
      std::string outpath = run_information.out_path + "/" + output_folder + "/";
      write_state(integrated, outpath, "output.csv");
      write_state(potential, outpath, "potential.csv");
      write_state(xcos_t, outpath, "x_t.csv");
      write_state(ycos_t, outpath, "y_t.csv");
      write_state(zcos_t, outpath, "z_t.csv");
      write_state(xcos_s, outpath, "x_s.csv");
      write_state(ycos_s, outpath, "y_s.csv");
      write_state(zcos_s, outpath, "z_s.csv");
      write_state(area, outpath, "areas.csv");
    }
  }

  MPI_Finalize();
  return 0;
}
