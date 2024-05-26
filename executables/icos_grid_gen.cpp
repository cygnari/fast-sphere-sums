#include <iostream>
#include <vector>
#include <cmath>

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
  int levels = 8;
  int point_count = 10 * pow(4, levels-1) + 2;
  int tri_count = 2 * (point_count - 2);
  // std::cout << point_count << std::endl;
  std::vector<double> xco (point_count, 0);
  std::vector<double> yco (point_count, 0);
  std::vector<double> zco (point_count, 0);
  std::vector<double> areas (point_count, 0);
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

  std::vector<std::vector<std::vector<int>>> dynamics_triangles;
  std::vector<std::vector<bool>> dynamics_triangles_is_leaf;
  std::vector<std::vector<int>> dynamics_points_parents;

  dynamics_triangles.resize(levels);
  dynamics_triangles[0] = std::vector<std::vector<int>>(20, std::vector<int>(4, 0));
  dynamics_triangles_is_leaf.resize(levels);
  dynamics_points_parents.resize(point_count, std::vector<int>(2, -1));

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

  std::string prefix = DATA_DIR + std::to_string(point_count) + "_icos_grid_";
  write_state(xco, prefix, "x.csv");
  write_state(yco, prefix, "y.csv");
  write_state(zco, prefix, "z.csv");
  write_state(areas, prefix, "areas.csv");

  return 0;
}
