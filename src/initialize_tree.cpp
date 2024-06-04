#include <vector>
#include <iostream>
#include <queue>
#include <cmath>

#include "general_utils.hpp"
#include "structs.hpp"

void initialize_cube_tree(const RunConfig& run_information, std::vector<CubePanel>& cube_panels, const std::vector<double>& x,
        const std::vector<double>& y, const std::vector<double>& z) {
  // sets up fast sum cubed sphere
  cube_panels.resize(6);
  for (int i = 0; i < 6; i++) {
    cube_panels[i].id = i;
    cube_panels[i].min_xi = -M_PI / 4;
    cube_panels[i].max_xi = M_PI / 4;
    cube_panels[i].min_eta = -M_PI / 4;
    cube_panels[i].max_eta = M_PI / 4;
    cube_panels[i].face = i+1;
  }

  double xval, yval, zval;
  std::vector<double> point_on_cube;
  bool found, point_found;
  int face;

  for (int i = 0; i < run_information.point_count; i++) {
    // assign points to base cube
    xval = x[i];
    yval = y[i];
    zval = z[i];
    face = face_from_xyz(xval, yval, zval);
    // std::cout << face << std::endl;
    cube_panels[face-1].point_count += 1;
    cube_panels[face-1].points_inside.push_back(i);
  }

  double mid_xi, mid_eta, min_xi, max_xi, min_eta, max_eta, xi, eta;
  std::vector<double> xieta;
  std::vector<int> points_to_assign;
  int start, end, which_panel, point_count;

  // iterate through panels, refine where needed
  for (int i = 0; i < cube_panels.size(); i++) {
    // check if panel needs to be divided
    // std::cout << i << "," << cube_panels[i].point_count << "," << cube_panels[i].face << std::endl;
    if ((cube_panels[i].point_count >= run_information.fast_sum_cluster_thresh) and (cube_panels[i].is_leaf)) {
      // refine the panel
      cube_panels[i].is_leaf = false;
      CubePanel sub_panel1, sub_panel2, sub_panel3, sub_panel4;
      min_xi = cube_panels[i].min_xi, max_xi = cube_panels[i].max_xi, min_eta = cube_panels[i].min_eta, max_eta = cube_panels[i].max_eta;
      mid_xi = (min_xi + max_xi) * 0.5;
      mid_eta = (min_eta + max_eta) * 0.5;
      sub_panel1.min_xi = mid_xi, sub_panel1.max_xi = max_xi, sub_panel1.min_eta = mid_eta, sub_panel1.max_eta = max_eta;
      sub_panel2.min_xi = min_xi, sub_panel2.max_xi = mid_xi, sub_panel2.min_eta = mid_eta, sub_panel2.max_eta = max_eta;
      sub_panel3.min_xi = min_xi, sub_panel3.max_xi = mid_xi, sub_panel3.min_eta = min_eta, sub_panel3.max_eta = mid_eta;
      sub_panel4.min_xi = mid_xi, sub_panel4.max_xi = max_xi, sub_panel4.min_eta = min_eta, sub_panel4.max_eta = mid_eta;
      sub_panel1.level = cube_panels[i].level + 1;
      sub_panel2.level = cube_panels[i].level + 1;
      sub_panel3.level = cube_panels[i].level + 1;
      sub_panel4.level = cube_panels[i].level + 1;
      start = cube_panels.size();
      sub_panel1.id = start;
      sub_panel2.id = start+1;
      sub_panel3.id = start+2;
      sub_panel4.id = start+3;
      sub_panel1.face = cube_panels[i].face;
      sub_panel2.face = cube_panels[i].face;
      sub_panel3.face = cube_panels[i].face;
      sub_panel4.face = cube_panels[i].face;
      sub_panel1.parent_panel = cube_panels[i].id;
      sub_panel2.parent_panel = cube_panels[i].id;
      sub_panel3.parent_panel = cube_panels[i].id;
      sub_panel4.parent_panel = cube_panels[i].id;
      cube_panels[i].child_panel_1 = start;
      cube_panels[i].child_panel_2 = start+1;
      cube_panels[i].child_panel_3 = start+2;
      cube_panels[i].child_panel_4 = start+3;
      cube_panels.push_back(sub_panel1);
      cube_panels.push_back(sub_panel2);
      cube_panels.push_back(sub_panel3);
      cube_panels.push_back(sub_panel4);
      end = cube_panels.size();

      // assign points in parent panel to child panels
      points_to_assign = cube_panels[i].points_inside;
      for (int j = 0; j < points_to_assign.size(); j++) {
        xval = x[points_to_assign[j]];
        yval = y[points_to_assign[j]];
        zval = z[points_to_assign[j]];
        xieta = xieta_from_xyz(xval, yval, zval, cube_panels[i].face);
        xi = xieta[0];
        eta = xieta[1];
        min_xi = cube_panels[i].min_xi, max_xi = cube_panels[i].max_xi, min_eta = cube_panels[i].min_eta, max_eta = cube_panels[i].max_eta;
        mid_xi = (min_xi + max_xi) * 0.5;
        mid_eta = (min_eta + max_eta) * 0.5;
        if (xi < mid_xi) {
          if (eta < mid_eta) {
            which_panel = 2;
          } else {
            which_panel = 1;
          }
        } else {
          if (eta < mid_eta) {
            which_panel = 3;
          } else {
            which_panel = 0;
          }
        }
        cube_panels[start+which_panel].point_count += 1;
        cube_panels[start+which_panel].points_inside.push_back(points_to_assign[j]);
      }
      point_count = cube_panels[start].point_count + cube_panels[start+1].point_count + cube_panels[start+2].point_count + cube_panels[start+3].point_count;
      if (point_count != cube_panels[i].point_count) {
        throw std::runtime_error("Not all points assigned, line 313");
      }
    }
  }

  std::vector<double> v1, v2, v3, v4, vc;
  double d1, d2, d3, d4, d;
  // compute panel radii
  for (int i = 0; i < cube_panels.size(); i++) {
    min_xi = cube_panels[i].min_xi, max_xi = cube_panels[i].max_xi;
    min_eta = cube_panels[i].min_eta, max_eta = cube_panels[i].max_eta;
    mid_xi = (min_xi+max_xi)*0.5;
    mid_eta = (min_eta+max_eta)*0.5;
    cube_panels[i].mid_xi = mid_xi, cube_panels[i].mid_eta = mid_eta;
    v1 = xyz_from_xieta(min_xi, min_eta, cube_panels[i].face);
    v2 = xyz_from_xieta(min_xi, max_eta, cube_panels[i].face);
    v3 = xyz_from_xieta(max_xi, min_eta, cube_panels[i].face);
    v4 = xyz_from_xieta(max_xi, max_eta, cube_panels[i].face);
    vc = xyz_from_xieta(mid_xi, mid_eta, cube_panels[i].face);
    d1 = gcdist(v1, vc, run_information.radius);
    d2 = gcdist(v2, vc, run_information.radius);
    d3 = gcdist(v3, vc, run_information.radius);
    d4 = gcdist(v4, vc, run_information.radius);
    d = std::max(std::max(std::max(d1, d2), d3), d4);
    cube_panels[i].radius = d;
  }
}

void dual_tree_traversal(const RunConfig& run_information, std::vector<InteractPair>& interactions, const std::vector<CubePanel>& cube_panels) {
  std::queue<int> target_squares;
  std::queue<int> source_squares;
  // top level interactions
  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 6; j++) {
      target_squares.push(i);
      source_squares.push(j);
    }
  }

  int index_target, index_source;
  std::vector<double> c1, c2;
  double dist, separation;
  while (target_squares.size() > 0) {
    // dual tree traversal
    index_target = target_squares.front();
    index_source = source_squares.front();
    target_squares.pop();
    source_squares.pop();
    if (cube_panels[index_target].point_count == 0) continue;
    if (cube_panels[index_source].point_count == 0) continue;
    c1 = xyz_from_xieta(cube_panels[index_target].mid_xi, cube_panels[index_target].mid_eta, cube_panels[index_target].face);
    c2 = xyz_from_xieta(cube_panels[index_source].mid_xi, cube_panels[index_source].mid_eta, cube_panels[index_source].face);
    dist = gcdist(c1, c2, run_information.radius);
    separation = (cube_panels[index_target].radius + cube_panels[index_source].radius) / dist;
    if ((dist > 0) and (separation < run_information.fast_sum_theta)) {
      // well separated
      InteractPair new_interact = {index_target, index_source, 0};
      if (cube_panels[index_target].point_count > run_information.fast_sum_cluster_thresh) {
        new_interact.interact_type += 2;
      }
      if (cube_panels[index_source].point_count > run_information.fast_sum_cluster_thresh) {
        new_interact.interact_type += 1;
      }
      interactions.push_back(new_interact);
    } else {
      // not well separated
      if ((cube_panels[index_target].point_count < run_information.fast_sum_cluster_thresh) and (cube_panels[index_source].point_count < run_information.fast_sum_cluster_thresh)) {
        // both have few points
        InteractPair new_interact = {index_target, index_source, 0};
        interactions.push_back(new_interact);
      } else if (cube_panels[index_target].is_leaf and cube_panels[index_source].is_leaf) {
        // both are leaves
        InteractPair new_interact = {index_target, index_source, 0};
        interactions.push_back(new_interact);
      } else if (cube_panels[index_target].is_leaf) {
        // target is leaf, break apart source
        target_squares.push(index_target);
        target_squares.push(index_target);
        target_squares.push(index_target);
        target_squares.push(index_target);
        source_squares.push(cube_panels[index_source].child_panel_1);
        source_squares.push(cube_panels[index_source].child_panel_2);
        source_squares.push(cube_panels[index_source].child_panel_3);
        source_squares.push(cube_panels[index_source].child_panel_4);
      } else if (cube_panels[index_source].is_leaf) {
        // source is leaf, break apart target
        source_squares.push(index_source);
        source_squares.push(index_source);
        source_squares.push(index_source);
        source_squares.push(index_source);
        target_squares.push(cube_panels[index_target].child_panel_1);
        target_squares.push(cube_panels[index_target].child_panel_2);
        target_squares.push(cube_panels[index_target].child_panel_3);
        target_squares.push(cube_panels[index_target].child_panel_4);
      } else {
        // neither panel is a leaf, refine the panel with more points
        if (cube_panels[index_target].point_count >= cube_panels[index_source].point_count) {
          // target has more points, break apart target
          source_squares.push(index_source);
          source_squares.push(index_source);
          source_squares.push(index_source);
          source_squares.push(index_source);
          target_squares.push(cube_panels[index_target].child_panel_1);
          target_squares.push(cube_panels[index_target].child_panel_2);
          target_squares.push(cube_panels[index_target].child_panel_3);
          target_squares.push(cube_panels[index_target].child_panel_4);
        } else {
          // source has more points, break apart
          target_squares.push(index_target);
          target_squares.push(index_target);
          target_squares.push(index_target);
          target_squares.push(index_target);
          source_squares.push(cube_panels[index_source].child_panel_1);
          source_squares.push(cube_panels[index_source].child_panel_2);
          source_squares.push(cube_panels[index_source].child_panel_3);
          source_squares.push(cube_panels[index_source].child_panel_4);
        }
      }
    }
  }
}

void dual_tree_traversal(const RunConfig& run_information, std::vector<InteractPair>& interactions, const std::vector<CubePanel>& cube_panels_target, const std::vector<CubePanel>& cube_panels_source) {
  std::queue<int> target_squares;
  std::queue<int> source_squares;
  // top level interactions
  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 6; j++) {
      target_squares.push(i);
      source_squares.push(j);
    }
  }

  int index_target, index_source;
  std::vector<double> c1, c2;
  double dist, separation;
  while (target_squares.size() > 0) {
    // dual tree traversal
    index_target = target_squares.front();
    index_source = source_squares.front();
    target_squares.pop();
    source_squares.pop();
    if (cube_panels_target[index_target].point_count == 0) continue;
    if (cube_panels_source[index_source].point_count == 0) continue;
    c1 = xyz_from_xieta(cube_panels_target[index_target].mid_xi, cube_panels_target[index_target].mid_eta, cube_panels_target[index_target].face);
    c2 = xyz_from_xieta(cube_panels_source[index_source].mid_xi, cube_panels_source[index_source].mid_eta, cube_panels_source[index_source].face);
    dist = gcdist(c1, c2, run_information.radius);
    separation = (cube_panels_target[index_target].radius + cube_panels_source[index_source].radius) / dist;
    // separation = c1[0]*c2[0] + c1[1]*c2[1] + c1[2]*c2[2];
    if ((dist > 0) and (separation < run_information.fast_sum_theta)) {
    // if (separation < 0.25) {
      // well separated
      InteractPair new_interact = {index_target, index_source, 0};
      if (cube_panels_target[index_target].point_count > run_information.fast_sum_cluster_thresh) {
        new_interact.interact_type += 2;
      }
      if (cube_panels_source[index_source].point_count > run_information.fast_sum_cluster_thresh) {
        new_interact.interact_type += 1;
      }
      interactions.push_back(new_interact);
    } else {
      // not well separated
      if ((cube_panels_target[index_target].point_count < run_information.fast_sum_cluster_thresh) and (cube_panels_source[index_source].point_count < run_information.fast_sum_cluster_thresh)) {
        // both have few points
        InteractPair new_interact = {index_target, index_source, 0};
        interactions.push_back(new_interact);
      } else if (cube_panels_target[index_target].is_leaf and cube_panels_source[index_source].is_leaf) {
        // both are leaves
        InteractPair new_interact = {index_target, index_source, 0};
        interactions.push_back(new_interact);
      } else if (cube_panels_target[index_target].is_leaf) {
        // target is leaf, break apart source
        target_squares.push(index_target);
        target_squares.push(index_target);
        target_squares.push(index_target);
        target_squares.push(index_target);
        source_squares.push(cube_panels_source[index_source].child_panel_1);
        source_squares.push(cube_panels_source[index_source].child_panel_2);
        source_squares.push(cube_panels_source[index_source].child_panel_3);
        source_squares.push(cube_panels_source[index_source].child_panel_4);
      } else if (cube_panels_source[index_source].is_leaf) {
        // source is leaf, break apart target
        source_squares.push(index_source);
        source_squares.push(index_source);
        source_squares.push(index_source);
        source_squares.push(index_source);
        target_squares.push(cube_panels_target[index_target].child_panel_1);
        target_squares.push(cube_panels_target[index_target].child_panel_2);
        target_squares.push(cube_panels_target[index_target].child_panel_3);
        target_squares.push(cube_panels_target[index_target].child_panel_4);
      } else {
        // neither panel is a leaf, refine the panel with more points
        if (cube_panels_target[index_target].point_count >= cube_panels_source[index_source].point_count) {
          // target has more points, break apart target
          source_squares.push(index_source);
          source_squares.push(index_source);
          source_squares.push(index_source);
          source_squares.push(index_source);
          target_squares.push(cube_panels_target[index_target].child_panel_1);
          target_squares.push(cube_panels_target[index_target].child_panel_2);
          target_squares.push(cube_panels_target[index_target].child_panel_3);
          target_squares.push(cube_panels_target[index_target].child_panel_4);
        } else {
          // source has more points, break apart
          target_squares.push(index_target);
          target_squares.push(index_target);
          target_squares.push(index_target);
          target_squares.push(index_target);
          source_squares.push(cube_panels_source[index_source].child_panel_1);
          source_squares.push(cube_panels_source[index_source].child_panel_2);
          source_squares.push(cube_panels_source[index_source].child_panel_3);
          source_squares.push(cube_panels_source[index_source].child_panel_4);
        }
      }
    }
  }
}
