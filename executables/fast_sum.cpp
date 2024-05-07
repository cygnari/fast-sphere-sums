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

  std::vector<double> test_vec {-0.5, -0.4, 1};
  std::vector<double> on_sphere = cube_to_sphere(test_vec);
  std::cout << on_sphere[0] << "," << on_sphere[1] << "," << on_sphere[2] << std::endl;
  std::vector<double> on_cube = sphere_to_cube(on_sphere);
  std::cout << on_cube[0] << "," << on_cube[1] << "," << on_cube[2] << std::endl;


  MPI_Finalize();
  return 0;
}
