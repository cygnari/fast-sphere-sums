#ifndef H_STRUCTS_H
#define H_STRUCTS_H

#include <vector>

struct CubePanel {
  int level  = 0; // refinement level, level 0 = base cube
  bool is_leaf = true; // whether or not this is a leaf panel
  int id; // location in vector of cube panels
  int parent_panel = -1;
  int child_panel_1 = -1;
  int child_panel_2 = -1;
  int child_panel_3 = -1;
  int child_panel_4 = -1;
  std::vector<double> vertex_1 {0, 0, 0};
  std::vector<double> vertex_2 {0, 0, 0};
  std::vector<double> vertex_3 {0, 0, 0};
  std::vector<double> vertex_4 {0, 0, 0};
  std::vector<double> vertex_1_sphere {0, 0, 0};
  std::vector<double> vertex_2_sphere {0, 0, 0};
  std::vector<double> vertex_3_sphere {0, 0, 0};
  std::vector<double> vertex_4_sphere {0, 0, 0};
  double radius;
  int point_count = 0; // number of points in this panel

};

struct IcosPanel {
  int level = 0; // refinement level, level 0 = base icosahedron
  bool is_leaf = true; // whether or not this is a leaf panel
  int id; // location in vector of icos panels
  int parent_panel = -1;
  int child_panel_1 = -1;
  int child_panel_2 = -1;
  int child_panel_3 = -1;
  int child_panel_4 = -1;
  std::vector<double> vertex_1 {0, 0, 0};
  std::vector<double> vertex_2 {0, 0, 0};
  std::vector<double> vertex_3 {0, 0, 0};
  std::vector<double> center_p {0, 0, 0};
  double radius;
  int point_count = 0; // points inside the panel
  std::vector<int> points_inside;
};

struct InteractPair {
  int index_target;
  int index_source;
  int interact_type; // 0 for PP, 1 for PC, 2 for CP, 3 for CC
};

#endif
