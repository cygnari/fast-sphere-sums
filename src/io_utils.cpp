#include "structs.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>

#include "fast-sphere-sums-config.h"

void read_run_config(const std::string file_name, RunConfig &run_information) {
  // reads run information of file_name
  std::ifstream config_file(file_name);
  if (config_file.fail()) {
      std::cout << "namelist file at " << file_name << std::endl;
      throw std::runtime_error("namelist not found, try running from inside /bin");
  }
  std::string line, word1, word2;

  while (true) {
    getline(config_file, line);
    std::stringstream str1(line);
    getline(str1, word1, '=');
    getline(str1, word2);
    if (word1 == "use_fast") {
      if (stoi(word2) == 1) {
        run_information.use_fast = true;
      }
    } else if (word1 == "out_path") {
      run_information.out_path = word2;
    } else if (word1 == "write_output") {
      if (stoi(word2) == 1) {
        run_information.write_output = true;
      }
    } else if (word1 == "write_precision") {
      run_information.write_precision = stoi(word2);
    } else if (word1 == "radius") {
      run_information.radius = stod(word2);
    } else if (word1 == "point_count") {
      run_information.point_count = stoi(word2);
    } else if (word1 == "fast_sum_cluster_thresh") {
      run_information.fast_sum_cluster_thresh = stoi(word2);
    } else if (word1 == "theta") {
      run_information.fast_sum_theta = stod(word2);
    } else if (word1 == "rotate") {
      if (stoi(word2) == 0) {
        run_information.rotate = false;
      }
    } else if (word1 == "initial_condition") {
      run_information.initial_condition = word2;
    } else if (word1 == "icp1") {
      run_information.init_cond_param1 = stoi(word2);
    } else if (word1 == "icp2") {
      run_information.init_cond_param2 = stoi(word2);
    } else if (word1 == "time") {
      run_information.time = stoi(word2);
    } else if (word1 == "degree") {
      run_information.interp_degree = stoi(word2);
    } else if (word1 == "alph") {
      run_information.alph = stod(word2);
    } else if (word1 == "beta") {
      run_information.beta = stod(word2);
    } else if (word1 == "gamm") {
      run_information.gamm = stod(word2);
    } else if (word1 == "sph_harm_comps") {
      run_information.sph_harm_comps = stoi(word2);
    } else if (word1 == "use_ces") {
      if (stoi(word2) == 1) {
        run_information.use_cesaro = true;
      }
    } else if (word1 == "tree") {
      if (word2 == "icos") {
        run_information.use_icos = true;
      }
    } else if (word1 == "grid") {
      run_information.grid = word2;
    } else {
      if (run_information.use_icos) {
        run_information.interp_point_count = (run_information.interp_degree + 1) * (run_information.interp_degree + 2) / 2;
      } else {
        run_information.interp_point_count = pow(run_information.interp_degree + 1, 2);
      }
      if (run_information.fast_sum_cluster_thresh == -1) {
        run_information.fast_sum_cluster_thresh = run_information.interp_point_count + 2;
      }

      return;
    }
  }
}

std::string create_config(const RunConfig &run_information, const std::string extra) {
  std::stringstream ss1, ss2, ss3;
  int precision;
  std::string output_filename = std::to_string(run_information.time) + "_" +
      std::to_string(run_information.point_count) + "_" + std::to_string(run_information.sph_harm_comps) + "_";
  if (run_information.use_fast) {
    output_filename +=
        "fast_" + std::to_string(run_information.fast_sum_cluster_thresh) + "_" +
        std::to_string(run_information.fast_sum_theta).substr(0, 3) + "_" + std::to_string(run_information.interp_degree);
  } else
    output_filename += "direct";
  output_filename += extra;
  return output_filename;
}

void read_data_field(const int point_count, std::vector<double>& data, const std::string file_name) {
  // reads in double csv file_name to vector data
  std::ifstream file_in;
  // file_in.open(DATA_DIR + run_information.grid + "_" + file_name + ".csv");
  file_in.open(file_name);
  std::string line;
  // double d;
  for (int i = 0; i < point_count; i++) {
    getline(file_in, line, '\n');
    // std::cout << line << std::endl;
    // d = stod(line);
    data[i] = stod(line);
    // std::cout << i << "," << line << "," << d << "," << data[i] << std::endl;
  }
  file_in.close();
}

void load_llns(const int sph_harm_comps, std::vector<double>& llns) {
  // read the LLNs
  std::ifstream fin1, fin2;
  fin1.open(DATA_DIR + std::string("lln_hs.csv"));
  fin2.open(DATA_DIR + std::string("lln_ks.csv"));
  std::string line1, line2;
  for (int i = 0; i < sph_harm_comps; i++) {
    getline(fin1, line1, ',');
    getline(fin2, line2, ',');
    llns[i] = 1 + stod(line2) - stod(line1);
  }
  fin1.close();
  fin2.close();
}

void write_state(const std::vector<double> &data, const std::string path) {
  // write sal potential
  std::ofstream write_out(path, std::ofstream::out | std::ofstream::trunc);
  for (int i = 0; i < data.size(); i++) { // write out state
    write_out << data[i] << "\n";
  }
  write_out.close();
}
