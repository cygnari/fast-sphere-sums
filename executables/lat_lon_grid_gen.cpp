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
  // generates an lat lon grid with grid spacing deg_sep
  double deg_sep_between_points = 0.5;
  int lat_points = static_cast<int>(180/deg_sep_between_points);
  int lon_points = static_cast<int>(360/deg_sep_between_points);
  std::vector<double> lat_vals (lat_points, 0);
  std::vector<double> lon_vals (lon_points, 0);
  std::vector<double> lat_area_part (lat_points, 0);
  std::vector<double> lon_area_part (lon_points, 0);
  int point_count = lat_points*lon_points;
  std::cout << "lat size: " << lat_points << ", lon size: " << lon_points << ", total size: " << point_count << std::endl;
  std::vector<double> xco (point_count, 0);
  std::vector<double> yco (point_count, 0);
  std::vector<double> zco (point_count, 0);
  std::vector<double> areas (point_count, 0);

  double min_lat = -90.0+deg_sep_between_points/2.0;
  double min_lon = deg_sep_between_points/2.0;
  for (int i = 0; i < lat_points; i++) {
    lat_vals[i] = min_lat + i*deg_sep_between_points;
    lat_area_part[i] = std::abs(sin(M_PI/180.0*(lat_vals[i]-deg_sep_between_points/2.0))-sin(M_PI/180.0*(lat_vals[i]+deg_sep_between_points/2.0)));
  }
  for (int i = 0; i < lon_points; i++) {
    lon_vals[i] = min_lon + i*deg_sep_between_points;
    lon_area_part[i] = M_PI*deg_sep_between_points/180.0;
  }

  int index = 0;
  double lon, lat, lon_area, lat_area, total_area = 0;
  std::vector<double> xyz (3, 0);
  for (int i = 0; i < lat_points; i++) {
    lat = lat_vals[i]*M_PI/180.0;
    lat_area = lat_area_part[i];
    for (int j = 0; j < lon_points; j++) {
      // if (index == 128525) {
      //   std::cout << "lat: " << lat_vals[i] << ", lon: " << lon_vals[j] << std::endl;
      //   std::cout << i << "," << j << std::endl;
      // }
      // if (index == 130325) {
      //   std::cout << "lat: " << lat_vals[i] << ", lon: " << lon_vals[j] << std::endl;
      //   std::cout << i << "," << j << std::endl;
      // }
      lon = lon_vals[j]*M_PI/180.0;
      lon_area = lon_area_part[j];
      xyz = latlon_to_xyz(lat, lon, 1.0);
      xco[index]=xyz[0];
      yco[index]=xyz[1];
      zco[index]=xyz[2];
      areas[index]=lat_area*lon_area;
      total_area += lat_area*lon_area;
      index++;
    }
  }

  std::cout << "total_area: " << total_area << ", 4pi: " << 4*M_PI << std::endl;

  std::string prefix = DATA_DIR + std::to_string(point_count) + "_ll_grid_";
  write_state(xco, prefix, "x.csv");
  write_state(yco, prefix, "y.csv");
  write_state(zco, prefix, "z.csv");
  write_state(areas, prefix, "areas.csv");

  return 0;
}
