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
  MPI_Win win_sal;

  MPI_Datatype dt_interaction;
  MPI_Type_contiguous(3, MPI_INT, &dt_interaction);
  MPI_Type_commit(&dt_interaction);

  std::chrono::steady_clock::time_point begin, end;

  RunConfig run_information;
  const std::string namelist_file = std::string(NAMELIST_DIR) + std::string("namelist.txt");
  read_run_config(namelist_file, run_information);
  run_information.mpi_P = P;
  run_information.mpi_ID = ID;
  // run_information.initial_condition = "ssh";

  std::vector<double> xcos (run_information.point_count, 0);
  std::vector<double> ycos (run_information.point_count, 0);
  std::vector<double> zcos (run_information.point_count, 0);
  std::vector<double> area (run_information.point_count, 0);
  std::vector<double> sshs (run_information.point_count, 0);
  std::vector<double> sals (run_information.point_count, 0);

  MPI_Win_create(&sals[0], run_information.point_count * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_sal);

  std::string data_pre = DATA_DIR + std::to_string(run_information.point_count) + "_" + run_information.grid + "_";

  read_data_field(run_information.point_count, xcos, data_pre + "x.csv");
  read_data_field(run_information.point_count, ycos, data_pre + "y.csv");
  read_data_field(run_information.point_count, zcos, data_pre + "z.csv");
  read_data_field(run_information.point_count, area, data_pre + "areas.csv");

  if (run_information.initial_condition != "ssh") {
    initialize_condition(run_information, xcos, ycos, zcos, sshs);
    if (run_information.balance_condition) {
      balance_conditions(sshs, area);
    }
  } else {
    read_data_field(run_information.point_count, sshs, data_pre + "ssh" + std::to_string(run_information.time) + ".csv");
  }

  if (run_information.use_fast) {
    // use fast summation
    std::vector<CubePanel> cube_panels;
    initialize_cube_tree(run_information, cube_panels, xcos, ycos, zcos);

    std::vector<InteractPair> interactions;
    begin = std::chrono::steady_clock::now();
    dual_tree_traversal(run_information, interactions, cube_panels);
    end = std::chrono::steady_clock::now();
    if (ID == 0) {
      std::cout << "tree traversal time: " << std::chrono::duration<double>(end - begin).count() << " seconds" << std::endl;
    }
    begin = std::chrono::steady_clock::now();
    fast_sum_sal(run_information, interactions, cube_panels, xcos, ycos, zcos, area, sshs, sals);
    sync_updates<double>(sals, P, ID, &win_sal, MPI_DOUBLE);
    end = std::chrono::steady_clock::now();
  } else {
    // direct summation
    bounds_determine_2d(run_information, P, ID);
    begin = std::chrono::steady_clock::now();
    direct_sum_sal(run_information, xcos, ycos, zcos, area, sshs, sals);
    sync_updates<double>(sals, P, ID, &win_sal, MPI_DOUBLE);
    end = std::chrono::steady_clock::now();
  }

  if (ID == 0) {
    std::cout << "integration time: " << std::chrono::duration<double>(end - begin).count() << " seconds" << std::endl;

    std::string output_folder = create_config(run_information) + "_sal";
    std::string filename = NAMELIST_DIR + std::string("initialize.py ") + run_information.out_path + "/" + output_folder;
    std::string command = "python ";
    command += filename;
    system(command.c_str());
    std::string outpath = run_information.out_path + "/" + output_folder + "/";
    write_state(sals, outpath, "sals.csv");
    write_state(sshs, outpath, "sshs.csv");
    write_state(xcos, outpath, "x.csv");
    write_state(ycos, outpath, "y.csv");
    write_state(zcos, outpath, "z.csv");
    write_state(area, outpath, "areas.csv");
  }

  MPI_Finalize();
  return 0;
}
