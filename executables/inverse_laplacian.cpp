#include <iostream>
#include <mpi.h>
#include <vector>
#include <cmath>

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

  std::chrono::steady_clock::time_point begin, end;

  RunConfig run_information;
  const std::string namelist_file = std::string(NAMELIST_DIR) + std::string("namelist.txt");
  read_run_config(namelist_file, run_information);

  std::vector<double> xcos (run_information.point_count, 0);
  std::vector<double> ycos (run_information.point_count, 0);
  std::vector<double> zcos (run_information.point_count, 0);
  std::vector<double> area (run_information.point_count, 0);
  std::vector<double> potential (run_information.point_count, 0);
  std::vector<double> integrated (run_information.point_count, 0);

  std::string data_pre = DATA_DIR + std::to_string(run_information.point_count) + "_" + run_information.grid + "_";

  read_data_field(run_information.point_count, xcos, data_pre + "x.csv");
  read_data_field(run_information.point_count, ycos, data_pre + "y.csv");
  read_data_field(run_information.point_count, zcos, data_pre + "z.csv");
  read_data_field(run_information.point_count, area, data_pre + "areas.csv");

  if (run_information.rotate) {
    rotate_points(xcos, ycos, zcos, run_information.alph, run_information.beta, run_information.gamm);
  }

  initialize_condition(run_information, xcos, ycos, zcos, potential);
  if (run_information.balance_condition) {
    balance_conditions(potential, area);
  }

  if (run_information.use_fast) {
    // use fast summation
    if (run_information.use_icos) {
      // icosahedral tree code
      std::vector<IcosPanel> icos_panels;
      initialize_icosahedron_tree(run_information, icos_panels, xcos, ycos, zcos);

      std::vector<InteractPair> interactions;
      dual_tree_traversal_icos(run_information, interactions, icos_panels);

      begin = std::chrono::steady_clock::now();
      fast_sum_inverse_laplacian_icos(run_information, interactions, icos_panels, xcos, ycos, zcos, area, potential, integrated);
      end = std::chrono::steady_clock::now();
    } else {
      // cubed sphere tree code
      std::vector<CubePanel> cube_panels;
      initialize_cube_tree(run_information, cube_panels, xcos, ycos, zcos);

      std::vector<InteractPair> interactions;
      dual_tree_traversal_cube(run_information, interactions, cube_panels);

      begin = std::chrono::steady_clock::now();
      fast_sum_inverse_laplacian_cube(run_information, interactions, cube_panels, xcos, ycos, zcos, area, potential, integrated);
      end = std::chrono::steady_clock::now();
    }
  } else {
    // direct summation
    begin = std::chrono::steady_clock::now();
    direct_sum_invert_laplacian(xcos, ycos, zcos, area, potential, integrated);
    end = std::chrono::steady_clock::now();
  }



  std::cout << "fast sum time: " << std::chrono::duration<double>(end - begin).count()
              << " seconds" << std::endl;

  std::string output_folder = create_config(run_information);

  std::string filename = NAMELIST_DIR + std::string("initialize.py ") + run_information.out_path + "/" + output_folder;
  std::string command = "python ";
  command += filename;
  system(command.c_str());
  std::string outpath = run_information.out_path + "/" + output_folder + "/";
  write_state(integrated, outpath, "output.csv");
  write_state(potential, outpath, "potential.csv");
  write_state(xcos, outpath, "x.csv");
  write_state(ycos, outpath, "y.csv");
  write_state(zcos, outpath, "z.csv");
  write_state(area, outpath, "areas.csv");
  // std::string outpath = run_information.out_path + "/" + output_folder + "/output.csv";
  // write_state(integrated, outpath);
  // std::string potpath = run_information.out_path + "/" + output_folder + "/potential.csv";
  // write_state(potential, potpath);
  // std::cout << interactions.size() << std::endl;


  MPI_Finalize();
  return 0;
}
