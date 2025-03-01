#ifndef H_STRUCTS_H
#define H_STRUCTS_H

#include <vector>
#include <string>

struct RunConfig {
  bool use_fast = false;
  bool use_fmm = false;
  bool use_cesaro = false; // use cesaro summed kernel or not
  std::string out_path;    // ../../run-output/ locally, on Derecho, /glade/derecho/scratch/achen/sal/
  int write_precision = 6; // number of decimal places, 6 for data visualization, 16 for error testing
  bool write_output = false;
  double radius = 1.0;
  int point_count;    // number of dynamics points
  // 2562, 10242, 40962, 163842, 655362, 2621442 for icos grid
  // full point count is 2313486 for mpas grid
  // 453657 points for mom6
  int interp_degree;      // interpolation degree
  int interp_point_count; // number of interpolation points
  int info_per_point; // how many doubles each point is, for example, storing x y z vor tracer = 5
  int time; // time 0, 1, 2, or 3
  int sph_harm_comps; // number of spherical harmonic components to use
  std::string grid; // "mpas_ocean_grid" or "icos_grid" or "cube_grid" or "ll_grid"
  bool rotate = true; // whether or not to rotate the points
  double alph = 0.01; // 3 rotation coefficients
  double beta = 0.01;
  double gamm = 0.01;
  double kernel_eps = 0.0; // kernel regularization parameter
  // epsilon for regularized kernels

  std::string initial_condition; // initial vorticity distribution
  int init_cond_param1;              // parameter for initial condition
  int init_cond_param2;           // parameter for initial condition
  bool balance_condition = true;       // balance initial condition
  // ICs: SH43, SH10

  // fast sum info
  int fast_sum_cluster_thresh; // threshold for a panel being a cluster
  double fast_sum_theta;       // well separated threshold

  // mpi info
  int mpi_P;       // total MPI ranks
  int mpi_ID;      // own MPI rank
  int point_lb; // range of assigned particles
  int point_ub;
  int point_own;
  int two_d_one_lb;
  int two_d_one_ub;
  int two_d_two_lb;
  int two_d_two_ub;
};

struct CubePanel {
  int level  = 0; // refinement level, level 0 = base cube
  bool is_leaf = true; // whether or not this is a leaf panel
  int id; // location in vector of cube panels
  int parent_panel = -1;
  int child_panel_1 = -1;
  int child_panel_2 = -1;
  int child_panel_3 = -1;
  int child_panel_4 = -1;
  int face; // face 1: x = 1, face 2: y = 1, face 3: x = -1, face 4: y = -1, face 5: z = 1, face 6: z = -1
  double min_xi;
  double mid_xi;
  double max_xi;
  double min_eta;
  double mid_eta;
  double max_eta;
  double radius;
  int point_count = 0; // number of points in this panel
  std::vector<int> points_inside;
};

struct InteractPair {
  int index_target;
  int index_source;
  int interact_type; // 0 for PP, 1 for PC, 2 for CP, 3 for CC
};

#endif
