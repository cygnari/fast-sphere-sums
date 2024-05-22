#ifndef H_IC_H
#define H_IC_H

#include <vector>
#include "structs.hpp"

void initialize_condition(const RunConfig& run_information, const std::vector<double>& xcos, const std::vector<double>& ycos, const std::vector<double>& zcos, std::vector<double>& potential);

void balance_conditions(std::vector<double>& potential, const std::vector<double>& area);

#endif
