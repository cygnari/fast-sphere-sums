#include <iostream>
#include <mpi.h>
#include <vector>
#include <cmath>
#include <chrono>

#include "fast-sphere-sums-config.h"
#include "direct_sum_funcs.hpp"
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

  // if (ID == 0) {
  //   begin = std::chrono::steady_clock::now();
  // }
  begin = std::chrono::steady_clock::now();

  direct_sum_invert_laplacian(xcos, ycos, zcos, area, potential, integrated);

  end = std::chrono::steady_clock::now();
  std::cout << "direct sum time: " << std::chrono::duration<double>(end - begin).count()
              << " seconds" << std::endl;

  std::string output_folder = create_config(run_information);

  std::string filename = NAMELIST_DIR + std::string("initialize.py ") + run_information.out_path + "/" + output_folder;
  std::string command = "python ";
  command += filename;
  system(command.c_str());
  std::string outpath = run_information.out_path + "/" + output_folder + "/output.csv";
  write_state(integrated, outpath);
  std::string potpath = run_information.out_path + "/" + output_folder + "/potential.csv";
  write_state(potential, potpath);

  MPI_Finalize();
  return 0;
}
