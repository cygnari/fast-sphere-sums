#include <iostream>
#include <mpi.h>
#include <vector>
#include <cmath>

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

  // double total_area;
  // for (int i = 0; i < run_information.point_count; i++) {
  //   total_area += area[i];
  // }
  // if (abs(4*M_PI - total_area) > 1e-16) {
  //   // throw std::runtime_error("Incorrect surface area");
  //   std::cout << 4*M_PI - total_area << std::endl;
  // }

  initialize_condition(run_information, xcos, ycos, zcos, potential);
  // balance_conditions(potential, area);

  direct_sum_invert_laplacian(xcos, ycos, zcos, area, potential, integrated);

  std::string output_folder = create_config(run_information);

  // std::cout << output_folder << std::endl;

  std::string filename = NAMELIST_DIR + std::string("initialize.py ") + run_information.out_path + "/" + output_folder;
  std::string command = "python ";
  command += filename;
  system(command.c_str());
  std::string outpath = run_information.out_path + "/" + output_folder + "/output.csv";
  write_state(integrated, outpath);
  std::string potpath = run_information.out_path + "/" + output_folder + "/potential.csv";
  write_state(potential, potpath);

  // std::cout << xcos[0] << "," << ycos[0] << "," << zcos[0] << std::endl;

  // double l2_error = 0, li_error = 0;
  // double diff;
  // for (int i = 0; i < run_information.point_count; i++) {
  //   // computes l2, linf error
  //   diff = integrated[i] -2.0*potential[i];
  //   diff = integrated[i];
  //   l2_error += diff*diff*area[i];
  //   li_error = std::max(li_error, abs(diff));
  // }
  // double l2_norm = 0, li_norm = 0;
  // for (int i = 0; i < run_information.point_count; i++) {
  //   l2_norm += area[i]*pow(2.0*potential[i], 2);
  //   // l2_norm += area[i]*pow(potential[i], 2);
  //   li_norm = std::max(li_norm, abs(potential[i]));
  // }
  //
  // std::cout << "Point count: " << run_information.point_count << std::endl;
  // std::cout << "abs l2 error: " << l2_error << std::endl;
  // std::cout << "abs li error: " << li_error << std::endl;
  // std::cout << "rel l2 error: " << l2_error / l2_norm << std::endl;
  // std::cout << "rel li error: " << li_error / li_norm << std::endl;

  MPI_Finalize();
  return 0;
}
