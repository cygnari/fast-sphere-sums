#ifndef H_IO_UTILS_H
#define H_IO_UTILS_H

#include "structs.hpp"
#include <string>

void read_run_config(const std::string file_name, RunConfig &run_information);

std::string create_config(const RunConfig &run_information, const std::string extra = "");

// void read_xyz(const RunConfig &run_information, std::vector<double>& xvals, std::vector<double>& yvals, std::vector<double>& zvals,
//                 std::vector<double>& areas);
//
// void read_latlon(const RunConfig &run_information, std::vector<double>& latvals, std::vector<double>& lonvals, std::vector<double>& areas);

void read_data_field(const RunConfig &run_information, std::vector<double>& data, std::string file_name);

void load_llns(int sph_harm_comps, std::vector<double>& llns);

void write_state(const std::vector<double> &data, const std::string path);

#endif
