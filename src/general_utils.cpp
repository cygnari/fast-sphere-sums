#include <vector>
#include <string>
#include <cmath>
#include <tuple>
#include <cassert>
#include <iostream>

extern "C" { // lapack
extern int dgesv_(int *, int *, double *, int *, int *, double *, int *, int *);
extern int dgelsd_(int *, int *, int *, double *, int *, double *, int *, double *, double *, int *, double *, int *, int *, int *);
extern int dgelss_(int *, int *, int *, double *, int *, double *, int *, double *, double *, int *, double *, int *, int *);
extern int ilaenv_(int *, char *, char *, int *, int *, int *, int *);
}

int linear_solve(const std::vector<double> &a_matrix, std::vector<double> &b_vec, int size, int nrhs, const int solver) {
  // solve Ax=b via dgesv or dgelsd or other options
  std::vector<double> interp_mat_copy = a_matrix; // copy matrix, some methods destroy the matrix
  int info;
  if (solver == 1) { // dgesv
    std::vector<int> ipiv (size, 0);
    dgesv_(&size, &nrhs, &*interp_mat_copy.begin(), &size, &*ipiv.begin(), &*b_vec.begin(), &size, &info);
  } else if (solver == 2) { // dgelss
    std::vector<double> sing_vals(size, 0);
    double rcond = -1;
    int lwork = std::max(3 * size + std::max(2*size, nrhs), 1);
    lwork *= 2;
    std::vector<double> work(lwork, 0);
    int rank;
    dgelss_(&size, &size, &nrhs, &*interp_mat_copy.begin(), &size, &*b_vec.begin(),
        &size, &*sing_vals.begin(), &rcond, &rank, &*work.begin(), &lwork, &info);
  } else if (solver == 3) { // dgelsd
    std::vector<double> sing_vals(size, 0);
    int rank;
    double rcond = -1;
    std::string name = "DGELSD";
    std::string empty = " ";
    int nine = 9;
    int zero = 0;
    int smallsize = ilaenv_(&nine, &name[0], &empty[0], &zero, &zero, &zero, &zero);
    int nlvl = std::max(1, int(log2(size/(smallsize+1)))+1);
    int lwork = 12*size+2*size*smallsize+8*size*nlvl+size*nrhs+pow(smallsize+1,2);
    lwork *= 2;
    int liwork = std::max(1, 3*size*nlvl+11*size);
    std::vector<double> work(lwork, 0);
    std::vector<int> iwork(liwork, 0);
    dgelsd_(&size, &size, &nrhs, &*interp_mat_copy.begin(), &size, &*b_vec.begin(), &size, &*sing_vals.begin(), &rcond, &rank, &*work.begin(), &lwork, &*iwork.begin(), &info);
  }
  return info;
}

std::vector<double> mat_vec_mult_3_3_3(const std::vector<std::vector<double>> &Amat,
                                        const double x, const double y, const double z) {
  std::vector<double> vec (3, 0);
  for (int i = 0; i < 3; i++) {
    vec[i] = Amat[i][0] * x + Amat[i][1] * y + Amat[i][2] * z;
  }
  return vec;
}

void rotate_points(std::vector<double>& xcos, std::vector<double>& ycos, std::vector<double>& zcos,
                    const double alph, const double beta, const double gamm) {
  // rotates points
  std::vector<std::vector<double>> rot_mat(3, std::vector<double>(3, 0));
  rot_mat[0][0] = cos(beta) * cos(gamm);
  rot_mat[0][1] = sin(alph) * sin(beta) * cos(gamm) - cos(alph) * sin(gamm);
  rot_mat[0][2] = cos(alph) * sin(beta) * cos(gamm) + sin(alph) * sin(gamm);
  rot_mat[1][0] = cos(beta) * sin(gamm);
  rot_mat[1][1] = sin(alph) * sin(beta) * sin(gamm) + cos(alph) * cos(gamm);
  rot_mat[1][2] = cos(alph) * sin(beta) * sin(gamm) - sin(alph) * cos(gamm);
  rot_mat[2][0] = -sin(beta);
  rot_mat[2][1] = sin(alph) * cos(beta);
  rot_mat[2][2] = cos(alph) * cos(beta);
  std::vector<double> rotated;
  for (int i = 0; i < xcos.size(); i++) {
    rotated = mat_vec_mult_3_3_3(rot_mat, xcos[i], ycos[i], zcos[i]);
    xcos[i] = rotated[0];
    ycos[i] = rotated[1];
    zcos[i] = rotated[2];
  }
}

void project_to_sphere(std::vector<double> &p1, const double radius) {
  // projects a point to the surface of a sphere of radius, modifies p1
  double norm = sqrt(p1[0]*p1[0] + p1[1]*p1[1] + p1[2]*p1[2]);
  for (int i = 0; i < p1.size(); i++) p1[i] *= radius / norm;
}

std::vector<double> project_to_sphere_2(const std::vector<double> p1, const double radius) {
  // projects a point to the surface of a sphere of radius,
  // returns the new coordinates
  double norm = sqrt(p1[0]*p1[0] + p1[1]*p1[1] + p1[2]*p1[2]);
  std::vector<double> newcords(p1.size());
  for (int i = 0; i < p1.size(); i++) newcords[i] = radius * p1[i] / norm;
  return newcords;
}

std::vector<double> latlon_to_xyz(const double lat, const double lon, const double radius) {
  // converts latitude, longitude to x y z coords
  std::vector<double> xyz (3, 0);
  double x, y, z, dist;
  x = radius * cos(lat) * cos(lon);
  y = radius * cos(lat) * sin(lon);
  z = radius * sin(lat);
  dist = sqrt(x * x + y * y + z * z);
  xyz[0] = x / dist;
  xyz[1] = y / dist;
  xyz[2] = z / dist;
  return xyz;
}

std::vector<double> xyz_to_latlon(const double x, const double y, const double z) {
  // turns cartesian coordinates to spherical coordinates
  double colat, lon;
  std::vector<double> latlon {0, 0};
  colat = atan2(sqrt(x * x + y * y), z); // colatitude
  lon = atan2(y, x);                   // longitude
  latlon[0] = M_PI / 2.0 - colat;
  latlon[1] = lon;
  return latlon;
}

std::vector<double> xyz_to_latlon(const std::vector<double>& point) {
  return xyz_to_latlon(point[0], point[1], point[2]);
}

std::vector<double> project_to_sphere(double x, double y, double z, const double radius) {
  // projects (x, y, z) to sphere of radius
  double point_dist = sqrt(x * x + y * y + z * z);
  x /= point_dist;
  y /= point_dist;
  z /= point_dist;
  return std::vector<double> {x, y, z};
}

double gcdist(const double x1, const double y1, const double z1, const double x2, const double y2, const double z2, const double radius) {
  double s = (x1 * x2 + y1 * y2 + z1 * z2) / (sqrt(x1 * x1 + y1 * y1 + z1 * z1) * sqrt(x2 * x2 + y2 * y2 + z2 * z2));
  double theta = acos(std::min(std::max(s, -1.0), 1.0));
  return theta * radius;
}

double gcdist(const double lat1, const double lon1, const double lat2, const double lon2, const double radius) {
  return radius * acos(std::min(1.0, std::max(-1.0, sin(lat1) * sin(lat2) +
                                cos(lat1) * cos(lat2) * cos(lon2 - lon1))));
}

double gcdist(const std::vector<double> p1, const std::vector<double> p2, const double radius) {
  return gcdist(p1[0], p1[1], p1[2], p2[0], p2[1], p2[2], radius);
}

std::vector<double> barycoords(const std::vector<double> &p1, const std::vector<double> &p2,
                               const std::vector<double> &p3, const double x, const double y, const double z) {
  // finds triangle barycentric coordinates of point p
  assert(p1.size() == 3);
  assert(p2.size() == 3);
  assert(p3.size() == 3);
  std::vector<double> coords {x, y, z};
  std::vector<double> mat{p1[0], p1[1], p1[2], p2[0], p2[1],
                          p2[2], p3[0], p3[1], p3[2]};
  int info = linear_solve(mat, coords, 3, 1, 1);
  if (info != 0) {
    throw std::runtime_error("Error in barycentric coordinate computation, line 94");
  }
  return coords;
}

double sphere_tri_area(const std::vector<double> &p1, const std::vector<double> &p2,
                       const std::vector<double> &p3, const double radius) {
  // finds the area of a spherical triangle
  std::vector<double> p1n, p2n, p3n;
  double a, b, c, s, z, area;
  p1n = p1; // don't modify p1, modify p1n
  p2n = p2;
  p3n = p3;
  p2n[0] -= p3[0], p2n[1] -= p3[1], p2n[2] -= p3[2];
  p3n[0] -= p1[0], p3n[1] -= p1[1], p3n[2] -= p1[2];
  p1n[0] -= p2[0], p1n[1] -= p2[1], p1n[2] -= p2[2];
  a = acos(1 - 0.5 * (p2n[0]*p2n[0]+p2n[1]*p2n[1]+p2n[2]*p2n[2]));
  b = acos(1 - 0.5 * (p3n[0]*p3n[0]+p3n[1]*p3n[1]+p3n[2]*p3n[2]));
  c = acos(1 - 0.5 * (p1n[0]*p1n[0]+p1n[1]*p1n[1]+p1n[2]*p1n[2]));
  s = (a + b + c) / 2;
  z = tan(s / 2) * tan((s - a) / 2) * tan((s - b) / 2) * tan((s - c) / 2);
  area = 4 * radius * radius * atan(sqrt(z));
  return area;
}

int face_from_xyz(const double x, const double y, const double z) {
  double ax = abs(x);
  double ay = abs(y);
  double az = abs(z);
  if ((ax >= ay) and (ax >= az)) {
    if (x >= 0) {
      return 1;
    } else {
      return 3;
    }
  } else if ((ay >= ax) and (ay >= az)) {
    if (y >= 0) {
      return 2;
    } else {
      return 4;
    }
  } else {
    if (z >= 0) {
      return 5;
    } else {
      return 6;
    }
  }
}

std::vector<double> xyz_from_xieta_1(const double xi, const double eta) {
  double X = tan(xi);
  double Y = tan(eta);
  std::vector<double> xyz (3, 0);
  xyz[0] = 1/sqrt(1+X*X+Y*Y);
  xyz[1] = X*xyz[0];
  xyz[2] = Y*xyz[0];
  return xyz;
}

std::vector<double> xyz_from_xieta_2(const double xi, const double eta) {
  double X = tan(xi);
  double Y = tan(eta);
  std::vector<double> xyz (3, 0);
  xyz[1] = 1/sqrt(1+X*X+Y*Y);
  xyz[0] = -X*xyz[1];
  xyz[2] = Y*xyz[1];
  return xyz;
}

std::vector<double> xyz_from_xieta_3(const double xi, const double eta) {
  double X = tan(xi);
  double Y = tan(eta);
  std::vector<double> xyz (3, 0);
  xyz[0] = -1/sqrt(1+X*X+Y*Y);
  xyz[1] = X*xyz[0];
  xyz[2] = -Y*xyz[0];
  return xyz;
}

std::vector<double> xyz_from_xieta_4(const double xi, const double eta) {
  double X = tan(xi);
  double Y = tan(eta);
  std::vector<double> xyz (3, 0);
  xyz[1] = -1/sqrt(1+X*X+Y*Y);
  xyz[0] = -X*xyz[1];
  xyz[2] = -Y*xyz[1];
  return xyz;
}

std::vector<double> xyz_from_xieta_5(const double xi, const double eta) {
  double X = tan(xi);
  double Y = tan(eta);
  std::vector<double> xyz (3, 0);
  xyz[2] = 1/sqrt(1+X*X+Y*Y);
  xyz[0] = -Y*xyz[2];
  xyz[1] = X*xyz[2];
  return xyz;
}

std::vector<double> xyz_from_xieta_6(const double xi, const double eta) {
  double X = tan(xi);
  double Y = tan(eta);
  std::vector<double> xyz (3, 0);
  xyz[2] = -1/sqrt(1+X*X+Y*Y);
  xyz[0] = -Y*xyz[2];
  xyz[1] = -X*xyz[2];
  return xyz;
}

std::vector<double> xyz_from_xieta(const double xi, const double eta, const int face) {
  std::vector<double> xyz;
  if (face == 1) {
    xyz = xyz_from_xieta_1(xi, eta);
  } else if (face == 2) {
    xyz = xyz_from_xieta_2(xi, eta);
  } else if (face == 3) {
    xyz =  xyz_from_xieta_3(xi, eta);
  } else if (face == 4) {
    xyz = xyz_from_xieta_4(xi, eta);
  } else if (face == 5) {
    xyz = xyz_from_xieta_5(xi, eta);
  } else if (face == 6) {
    xyz = xyz_from_xieta_6(xi, eta);
  } else {
    throw std::runtime_error("face is not 1 to 6, line 198");
  }
  return xyz;
}

std::vector<double> xieta_from_xyz_1(const double x, const double y, const double z) {
  std::vector<double> xieta (2, 0);
  xieta[0] = atan(y/x);
  xieta[1] = atan(z/x);
  return xieta;
}

std::vector<double> xieta_from_xyz_2(const double x, const double y, const double z) {
  std::vector<double> xieta (2, 0);
  xieta[0] = atan(-x/y);
  xieta[1] = atan(z/y);
  return xieta;
}

std::vector<double> xieta_from_xyz_3(const double x, const double y, const double z) {
  std::vector<double> xieta (2, 0);
  xieta[0] = atan(y/x);
  xieta[1] = atan(-z/x);
  return xieta;
}

std::vector<double> xieta_from_xyz_4(const double x, const double y, const double z) {
  std::vector<double> xieta (2, 0);
  xieta[0] = atan(-x/y);
  xieta[1] = atan(-z/y);
  return xieta;
}

std::vector<double> xieta_from_xyz_5(const double x, const double y, const double z) {
  std::vector<double> xieta (2, 0);
  xieta[0] = atan(y/z);
  xieta[1] = atan(-x/z);
  return xieta;
}

std::vector<double> xieta_from_xyz_6(const double x, const double y, const double z) {
  std::vector<double> xieta (2, 0);
  xieta[0] = atan(-y/z);
  xieta[1] = atan(-x/z);
  return xieta;
}

std::vector<double> xieta_from_xyz(const double x, const double y, const double z) {
  double ax = abs(x);
  double ay = abs(y);
  double az = abs(z);

  if ((ax >= ay) and (ax >= az)) {
    if (x >= 0) {
      return xieta_from_xyz_1(x, y, z);
    } else {
      return xieta_from_xyz_3(x, y, z);
    }
  } else if ((ay >= ax) and (ay >= az)) {
    if (y >= 0) {
      return xieta_from_xyz_2(x, y, z);
    } else {
      return xieta_from_xyz_4(x, y, z);
    }
  } else {
    if (z >= 0) {
      return xieta_from_xyz_5(x, y, z);
    } else {
      return xieta_from_xyz_6(x, y, z);
    }
  }
}

std::vector<double> xieta_from_xyz(const double x, const double y, const double z, const int face) {
  if (face == 1) {
    return xieta_from_xyz_1(x, y, z);
  } else if (face == 2) {
    return xieta_from_xyz_2(x, y, z);
  } else if (face == 3) {
    return xieta_from_xyz_3(x, y, z);
  } else if (face == 4) {
    return xieta_from_xyz_4(x, y, z);
  } else if (face == 5) {
    return xieta_from_xyz_5(x, y, z);
  } else if (face == 6) {
    return xieta_from_xyz_6(x, y, z);
  } else {
    throw std::runtime_error("Face not between 1 and 6, line 342");
  }
}

// std::vector<double> cube_to_sphere(const double x, const double y, const double z) {
//   // maps a point x, y, z on the unit cube to the surface of the sphere
//   std::vector<double> point {x, y, z};
//   double x2 = x*x;
//   double y2 = y*y;
//   double z2 = z*z;
//   point[0] *= sqrt(1-y2/2-z2/2+y2*z2/3);
//   point[1] *= sqrt(1-x2/2-z2/2+x2*z2/3);
//   point[2] *= sqrt(1-x2/2-y2/2+x2*y2/3);
//   return point;
// }
//
// std::vector<double> cube_to_sphere(const std::vector<double>& cube_point) {
//   double x = cube_point[0];
//   double y = cube_point[1];
//   double z = cube_point[2];
//   return cube_to_sphere(x, y, z);
// }
//
// std::vector<double> sphere_to_cube(const double x, const double y, const double z) {
//   // inverts the previous cube to sphere map
//   std::vector<double> cube_point {0, 0, 0};
//   double ax = abs(x);
//   double ay = abs(y);
//   double az = abs(z);
//   double a2, b2, inner, innersqrt;
//   if ((ax >= ay) and (ax >= az)) {
//     // x coordiante is the largest so one of the x faces
//     a2 = 2*y*y;
//     b2 = 2*z*z;
//     inner = -a2+b2-3;
//     innersqrt = -sqrt(inner*inner - 12*a2);
//     // set y coordinate
//     if (y == 0) {
//       cube_point[1] = 0;
//     } else {
//       cube_point[1] = std::copysign(std::min(1.0, sqrt(innersqrt + a2 - b2 + 3.0) / sqrt(2)), y);
//     }
//     // set z coordinate
//     if (z == 0) {
//       cube_point[2] = 0;
//     } else {
//       cube_point[2] = std::copysign(std::min(1.0, sqrt(innersqrt - a2 + b2 + 3.0) / sqrt(2)), z);
//     }
//     // set x coordinate
//     cube_point[0] = std::copysign(1, x);
//   } else if ((ay >= ax) and (ay >= az)) {
//     // y coordinate is the largest so one of the y faces
//     a2 = 2*x*x;
//     b2 = 2*z*z;
//     inner = -a2+b2-3;
//     innersqrt = -sqrt(inner*inner - 12*a2);
//     // set x coordinate
//     if (x == 0) {
//       cube_point[0] = 0;
//     } else {
//       cube_point[0] = std::copysign(std::min(1.0, sqrt(innersqrt + a2 - b2 + 3.0) / sqrt(2)), x);
//     }
//     // set z coordinate
//     if (z == 0) {
//       cube_point[2] = 0;
//     } else {
//       cube_point[2] = std::copysign(std::min(1.0, sqrt(innersqrt - a2 + b2 + 3.0) / sqrt(2)), z);
//     }
//     // set y coordinate
//     cube_point[1] = std::copysign(1, y);
//   } else {
//     // z coordiante is the largest so one of the z faces
//     a2 = 2*x*x;
//     b2 = 2*y*y;
//     inner = -a2+b2-3;
//     innersqrt = -sqrt(inner*inner - 12*a2);
//     // set x coordinate
//     if (x == 0) {
//       cube_point[0] = 0;
//     } else {
//       cube_point[0] = std::copysign(std::min(1.0, sqrt(innersqrt + a2 - b2 + 3.0) / sqrt(2)), x);
//     }
//     // set y coordinate
//     if (y == 0) {
//       cube_point[1] = 0;
//     } else {
//       cube_point[1] = std::copysign(std::min(1.0, sqrt(innersqrt - a2 + b2 + 3.0) / sqrt(2)), y);
//     }
//     // set z coordinate
//     cube_point[2] = std::copysign(1, z);
//   }
//   return cube_point;
// }
//
// std::vector<double> sphere_to_cube(const std::vector<double>& sphere_point) {
//   return sphere_to_cube(sphere_point[0], sphere_point[1], sphere_point[2]);
// }
//
// std::vector<double> cube_non_one (const double x, const double y, const double z) {
//   std::vector<double> out {0, 0};
//   if (x == 1) {
//     out[0] = y;
//     out[1] = z;
//   } else if (x == -1) {
//     out[0] = y;
//     out[1] = z;
//   } else if (y == 1) {
//     out[0] = x;
//     out[1] = z;
//   } else if (y == -1) {
//     out[0] = x;
//     out[1] = z;
//   } else if (z == 1) {
//     out[0] = x;
//     out[1] = y;
//   } else {
//     out[0] = x;
//     out[1] = y;
//   }
//   return out;
// }
//
// std::vector<double> cube_non_one(const std::vector<double>& cube_point) {
//   return cube_non_one(cube_point[0], cube_point[1], cube_point[2]);
// }
//
// bool check_in_cube_square(const double x, const double y, const double z, const CubePanel& cube_panel) {
//   // check if point is in cube panel
//   if (x == 1) {
//     if ((cube_panel.vertex_1[0] == 1) and (cube_panel.vertex_2[0] == 1) and (cube_panel.vertex_3[0] == 1) and (cube_panel.vertex_4[0] == 1)) {
//       if ((y > cube_panel.min_side_1) and (y < cube_panel.max_side_1) and (z > cube_panel.min_side_2) and (z < cube_panel.min_side_2)) {
//         return true;
//       }
//     }
//   } else if (x == -1) {
//     if ((cube_panel.vertex_1[0] == -1) and (cube_panel.vertex_2[0] == -1) and (cube_panel.vertex_3[0] == -1) and (cube_panel.vertex_4[0] == -1)) {
//       if ((y > cube_panel.min_side_1) and (y < cube_panel.max_side_1) and (z > cube_panel.min_side_2) and (z < cube_panel.min_side_2)) {
//         return true;
//       }
//     }
//   } else if (y == 1) {
//     if ((cube_panel.vertex_1[1] == 1) and (cube_panel.vertex_2[1] == 1) and (cube_panel.vertex_3[1] == 1) and (cube_panel.vertex_4[1] == 1)) {
//       if ((x > cube_panel.min_side_1) and (x < cube_panel.max_side_1) and (z > cube_panel.min_side_2) and (z < cube_panel.min_side_2)) {
//         return true;
//       }
//     }
//   } else if (y == -1) {
//     if ((cube_panel.vertex_1[1] == -1) and (cube_panel.vertex_2[1] == -1) and (cube_panel.vertex_3[1] == -1) and (cube_panel.vertex_4[1] == -1)) {
//       if ((x > cube_panel.min_side_1) and (x < cube_panel.max_side_1) and (z > cube_panel.min_side_2) and (z < cube_panel.min_side_2)) {
//         return true;
//       }
//     }
//   } else if (z == 1) {
//     if ((cube_panel.vertex_1[2] == 1) and (cube_panel.vertex_2[2] == 1) and (cube_panel.vertex_3[2] == 1) and (cube_panel.vertex_4[2] == 1)) {
//       if ((x > cube_panel.min_side_1) and (x < cube_panel.max_side_1) and (y > cube_panel.min_side_2) and (y < cube_panel.min_side_2)) {
//         return true;
//       }
//     }
//   } else if (z == -1) {
//     if ((cube_panel.vertex_1[2] == -1) and (cube_panel.vertex_2[2] == -1) and (cube_panel.vertex_3[2] == -1) and (cube_panel.vertex_4[2] == -1)) {
//       if ((x > cube_panel.min_side_1) and (x < cube_panel.max_side_1) and (y > cube_panel.min_side_2) and (y < cube_panel.min_side_2)) {
//         return true;
//       }
//     }
//   } else {
//     // throw std::runtime_error("Point is not on cube, check_in_cube_square, general_utils, line 309 ")
//     return false;
//   }
//   return false;
// }
//
// bool check_in_cube_square(const std::vector<double> point, const CubePanel& cube_panel) {
//   // check if point is in cube panel
//   return check_in_cube_square(point[0], point[1], point[2], cube_panel);
// }
//
// bool check_in_cube_square(const double x, const double y, const double z, const double min_side_1, const double max_side_1, const double min_side_2, const double max_side_2) {
//   std::vector<double> non_ones = cube_non_one(x, y, z);
//   if ((non_ones[0] > min_side_1) and (non_ones[0] < max_side_1)) {
//     if ((non_ones[1] > min_side_2) and (non_ones[1] < max_side_2)) {
//       return true;
//     }
//   }
//   return false;
// }
