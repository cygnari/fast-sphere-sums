#include <vector>
#include <cmath>

void direct_sum_invert_laplacian(const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos, const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral) {
  // perform direct summation to convolve to invert the laplacian
  double tx, ty, tz, sx, sy, sz;
  for (int i = 0; i < xcos.size(); i++) { // loop over targets
    tx = xcos[i], ty = ycos[i], tz = zcos[i];
    for (int j = 0; j < xcos.size(); j++) { // loop over sources
      if (i != j) {
      // skip singularity
      sx = xcos[j], sy = ycos[j], sz = zcos[j];
      integral[i] += -1.0/(4.0*M_PI) * log(1-tx*sx-ty*sy-tz*sz)*potential[j]*area[j];
      }
    }
  }
}
