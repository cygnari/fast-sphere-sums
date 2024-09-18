#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>

#include "fast-sphere-sums-config.h"
#include "general_utils.hpp"
#include "initialize_tree.hpp"
#include "io_utils.hpp"
#include "structs.hpp"

int check_point_exist(const std::vector<std::vector<int>> &parent_points,
                      const int point_count, const int iv1, const int iv2) {
  for (int i = 0; i < point_count; i++) {
    if ((parent_points[i][0] == iv1) and (parent_points[i][1] == iv2))
      return i;
  }
  return -1;
}

int main(int argc, char ** argv) {
  // generates an icosahedral grid with levels
  int levels = 2;
  int point_count = 6 * pow(4, levels);
  // int points_per_face = pow(4, levels);
  int points_per_side = pow(2, levels);
  // int tri_count = 2 * (point_count - 2);
  // std::cout << point_count << std::endl;
  std::vector<double> xco (point_count, 0);
  std::vector<double> yco (point_count, 0);
  std::vector<double> zco (point_count, 0);
  std::vector<double> areas (point_count, 0);
  std::vector<double> theta_phi (points_per_side, 0);

  double dist_between = M_PI / (2.0*points_per_side);
  double small_offset = dist_between/2.0;
  for (int i = 0; i < points_per_side; i++) {
    theta_phi[i] = -M_PI/4.0+small_offset+dist_between*i;
  }

  int index = 0;
  double theta, phi, area1, area2;
  std::vector<double> xyz (3, 0), xyz1 (3, 0), xyz2 (3, 0), xyz3 (3, 0), xyz4 (3, 0);
  for (int i = 0; i < points_per_side; i++) {
    theta = theta_phi[i];
    for (int j = 0; j < points_per_side; j++) {
      phi = theta_phi[j];
      for (int k = 1; k <= 6; k++) {
        xyz = xyz_from_xieta(theta, phi, k);
        xco[index] = xyz[0];
        yco[index] = xyz[1];
        zco[index] = xyz[2];
        xyz1 = xyz_from_xieta(theta-small_offset, phi-small_offset, k);
        xyz2 = xyz_from_xieta(theta+small_offset, phi-small_offset, k);
        xyz3 = xyz_from_xieta(theta+small_offset, phi+small_offset, k);
        xyz4 = xyz_from_xieta(theta-small_offset, phi+small_offset, k);
        area1 = sphere_tri_area(xyz1, xyz2, xyz3, 1.0);
        area2 = sphere_tri_area(xyz1, xyz3, xyz4, 1.0);
        areas[index] = area1+area2;
        index++;
      }
    }
  }

  std::string prefix = DATA_DIR + std::to_string(point_count) + "_cube_grid_";
  write_state(xco, prefix, "x.csv");
  write_state(yco, prefix, "y.csv");
  write_state(zco, prefix, "z.csv");
  write_state(areas, prefix, "areas.csv");

  return 0;
}
