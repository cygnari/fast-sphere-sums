#ifndef H_GENERAL_UTIL_H
#define H_GENERAL_UTIL_H

#include <vector>
#include "structs.hpp"

int linear_solve(const std::vector<double> &a_matrix, std::vector<double> &b_vec, int size, int nrhs, const int solver);

std::vector<double> latlon_to_xyz(const double lat, const double lon, const double radius);

std::vector<double> xyz_to_latlon(const double x, const double y, const double z);

std::vector<double> xyz_to_latlon(const std::vector<double>& point);

std::vector<double> project_to_sphere(double x, double y, double z, const double radius);

double gcdist(const double x1, const double y1, const double z1, const double x2, const double y2, const double z2, const double radius);

double gcdist(const double lat1, const double lon1, const double lat2, const double lon2, const double radius);

double gcdist(const std::vector<double> p1, const std::vector<double> p2, const double radius);

std::vector<double> barycoords(const std::vector<double> &p1,
                               const std::vector<double> &p2,
                               const std::vector<double> &p3,
                               const double x, const double y, const double z);

int face_from_xyz(const double x, const double y, const double z);

std::vector<double> xyz_from_xieta_1(const double xi, const double eta);

std::vector<double> xyz_from_xieta_2(const double xi, const double eta);

std::vector<double> xyz_from_xieta_3(const double xi, const double eta);

std::vector<double> xyz_from_xieta_4(const double xi, const double eta);

std::vector<double> xyz_from_xieta_5(const double xi, const double eta);

std::vector<double> xyz_from_xieta_6(const double xi, const double eta);

std::vector<double> xyz_from_xieta(const double xi, const double eta, const int face);

std::vector<double> xieta_from_xyz_1(const double x, const double y, const double z);

std::vector<double> xieta_from_xyz_2(const double x, const double y, const double z);

std::vector<double> xieta_from_xyz_3(const double x, const double y, const double z);

std::vector<double> xieta_from_xyz_4(const double x, const double y, const double z);

std::vector<double> xieta_from_xyz_5(const double x, const double y, const double z);

std::vector<double> xieta_from_xyz_6(const double x, const double y, const double z);

std::vector<double> xieta_from_xyz(const double x, const double y, const double z);

std::vector<double> xieta_from_xyz(const double x, const double y, const double z, const int face);

// std::vector<double> cube_to_sphere(const double x, const double y, const double z);
//
// std::vector<double> cube_to_sphere(const std::vector<double>& cube_point);
//
// std::vector<double> sphere_to_cube(const double x, const double y, const double z);
//
// std::vector<double> sphere_to_cube(const std::vector<double>& sphere_point);
//
// std::vector<double> cube_non_one (const double x, const double y, const double z);
//
// std::vector<double> cube_non_one(const std::vector<double>& cube_point);
//
// bool check_in_cube_square(const double x, const double y, const double z, const CubePanel& cube_panel);
//
// bool check_in_cube_square(const std::vector<double> point, const CubePanel& cube_panel);
//
// bool check_in_cube_square(const double x, const double y, const double z, const double min_side_1, const double max_side_1, const double min_side_2, const double max_side_2);

#endif
