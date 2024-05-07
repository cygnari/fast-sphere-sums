#ifndef H_GENERAL_UTIL_H
#define H_GENERAL_UTIL_H

#include <vector>

int linear_solve(const std::vector<double> &a_matrix, std::vector<double> &b_vec, int size, int nrhs, const int solver);

std::tuple<double, double, double> latlon_to_xyz(const double lat, const double lon, const double radius);

std::tuple<double, double> xyz_to_latlon(const double x, const double y, const double z);

std::tuple<double, double> xyz_to_latlon(const std::vector<double>& point);

std::vector<double> project_to_sphere(double x, double y, double z, const double radius);

double gcdist(const double x1, const double y1, const double z1, const double x2, const double y2, const double z2, const double radius);

double gcdist(const double lat1, const double lon1, const double lat2, const double lon2, const double radius);

std::vector<double> barycoords(const std::vector<double> &p1,
                               const std::vector<double> &p2,
                               const std::vector<double> &p3,
                               const double x, const double y, const double z);

std::vector<double> cube_to_sphere(const double x, const double y, const double z);

std::vector<double> cube_to_sphere(const std::vector<double>& cube_point);

std::vector<double> sphere_to_cube(const double x, const double y, const double z);

std::vector<double> sphere_to_cube(const std::vector<double>& sphere_point);

std::tuple<double, double> cube_non_one (const double x, const double y, const double z);

std::tuple<double, double> cube_non_one(const std::vector<double>& cube_point);

#endif
