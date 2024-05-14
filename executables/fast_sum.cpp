#include <iostream>
#include <mpi.h>
#include <vector>
#include <cmath>

#include "fast-sphere-sums-config.h"
#include "general_utils.hpp"
#include "initialize_tree.hpp"
#include "mpi_utils.hpp"
#include "structs.hpp"

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  int P, ID;
  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD, &P);
  MPI_Comm_rank(MPI_COMM_WORLD, &ID);
  // std::
  std::vector<double> xyz;
  std::vector<double> xieta;
  for (int i = 1; i < 7; i++) {
    xyz = xyz_from_xieta(0, 0, i);
    std::cout << "face " << i << " xyz " << xyz[0] << "," << xyz[1] << "," << xyz[2] << std::endl;
    xieta = xieta_from_xyz(xyz[0], xyz[1], xyz[2]);
    // xieta = xieta_from_xyz(xyz[0], xyz[1], xyz[2], i);
    std::cout << xieta[0] << "," << xieta[1] << std::endl;
  }


  MPI_Finalize();
  return 0;
}
