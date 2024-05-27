#ifndef H_DIRECT_SUM_FUNCS_H
#define H_DIRECT_SUM_FUNCS_H

#include <vector>
#include "structs.hpp"

void direct_sum_invert_laplacian(RunConfig& run_information, const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos, const std::vector<double>& area, const std::vector<double>& potential, std::vector<double>& integral);

#endif
