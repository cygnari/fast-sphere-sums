#include <cassert>
#include <mpi.h>
#include <vector>
#include <iostream>
#include "structs.hpp"
#include <cmath>
#define assertm(exp, msg) assert(((void)msg, exp))

bool test_is_same(const int x, MPI_Comm mpi_communicator = MPI_COMM_WORLD) { // test if all processes have the same value for a variable
  int p[2];
  p[0] = -x;
  p[1] = x;
  MPI_Allreduce(MPI_IN_PLACE, p, 2, MPI_INT, MPI_MIN, mpi_communicator);
  return (p[0] == -p[1]);
}

void bounds_determine_1d(RunConfig& run_information, const int P, const int ID) {
  // 1d parallel layout
  std::vector<int> particles(P, int(run_information.point_count / P));
  std::vector<int> lb(P, 0);
  std::vector<int> ub(P, 0);
  int total = P * int(run_information.point_count / P);
  int gap = run_information.point_count - total;
  for (int i = 1; i < gap + 1; i++) {
    particles[i] += 1;
  }
  total = 0;
  for (int i = 0; i < P; i++) {
    total += particles[i];
  }

  assertm(total == run_information.point_count, "Particle count not correct");

  ub[0] = particles[0];
  for (int i = 1; i < P; i++) {
    lb[i] = ub[i - 1];
    ub[i] = lb[i] + particles[i];
  }
  run_information.point_lb = lb[ID];
  run_information.point_ub = ub[ID];
  run_information.point_own = particles[ID];
}

void bounds_determine_2d(RunConfig& run_information, const int P, const int ID) {
  // 2d parallel layout
  int pow_2 = std::log2(P);
  int outer_count = pow(2, pow_2 / 2);
  int inner_count = P / outer_count;
  // assertm(P == outer_count*inner_count, "Processors not split up correctly");

  std::vector<int> particles_outer(outer_count, run_information.point_count / outer_count);
  std::vector<int> lb_outer(outer_count, 0);
  std::vector<int> ub_outer(outer_count, 0);
  int total_outer = outer_count * (run_information.point_count / outer_count);
  int gap_outer = run_information.point_count - total_outer;
  for (int i = 1; i < gap_outer+1; i++) {
    particles_outer[i] += 1;
  }
  total_outer = 0;
  for (int i = 0; i < outer_count; i++) {
    total_outer += particles_outer[i];
  }
  assertm(total_outer == run_information.point_count, "Particle count not correct");
  ub_outer[0] = particles_outer[0];
  for (int i = 1; i < outer_count; i++) {
    lb_outer[i] = ub_outer[i-1];
    ub_outer[i] = lb_outer[i]+particles_outer[i];
  }

  std::vector<int> particles_inner(inner_count, run_information.point_count / inner_count);
  std::vector<int> lb_inner(inner_count, 0);
  std::vector<int> ub_inner(inner_count, 0);
  int total_inner = inner_count * (run_information.point_count / inner_count);
  int gap_inner = run_information.point_count - total_inner;
  for (int i = 1; i < gap_inner+1; i++) {
    particles_inner[i] += 1;
  }
  total_inner = 0;
  for (int i = 0; i < inner_count; i++) {
    total_inner += particles_inner[i];
  }
  assertm(total_inner == run_information.point_count, "Particle count not correct");
  ub_inner[0] = particles_inner[0];
  for (int i = 1; i < inner_count; i++) {
    lb_inner[i] = ub_inner[i-1];
    ub_inner[i] = lb_inner[i]+particles_inner[i];
  }

  int index = 0;
  for (int i = 0; i < outer_count; i++) {
    for (int j = 0; j < inner_count; j++) {
      if (index == ID) {
        run_information.two_d_one_lb = lb_outer[i];
        run_information.two_d_one_ub = ub_outer[i];
        run_information.two_d_two_lb = lb_inner[j];
        run_information.two_d_two_ub = ub_inner[j];
      }
      index++;
    }
  }
}
