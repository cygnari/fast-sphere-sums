#include <vector>
#include <iostream>
#include <queue>
#include <cmath>

#include "general_utils.hpp"
#include "structs.hpp"

void initialize_icosahedron(const RunConfig& run_information, std::vector<IcosPanel>& icos_panels, const std::vector<double>& x,
        const std::vector<double>& y, const std::vector<double>& z) {
  // sets up fast sum icosahedron
  icos_panels.resize(20);
  double phi = (1 + sqrt(5)) / 2;
  std::vector<double> point_1 = project_to_sphere(0, 1, phi, run_information.radius);
  std::vector<double> point_2 = project_to_sphere(0, -1, phi, run_information.radius);
  std::vector<double> point_3 = project_to_sphere(0, 1, -phi, run_information.radius);
  std::vector<double> point_4 = project_to_sphere(0, -1, -phi, run_information.radius);
  std::vector<double> point_5 = project_to_sphere(1, phi, 0, run_information.radius);
  std::vector<double> point_6 = project_to_sphere(1, -phi, 0, run_information.radius);
  std::vector<double> point_7 = project_to_sphere(-1, phi, 0, run_information.radius);
  std::vector<double> point_8 = project_to_sphere(-1, -phi, 0, run_information.radius);
  std::vector<double> point_9 = project_to_sphere(phi, 0, 1, run_information.radius);
  std::vector<double> point_10 = project_to_sphere(phi, 0, -1, run_information.radius);
  std::vector<double> point_11 = project_to_sphere(-phi, 0, 1, run_information.radius);
  std::vector<double> point_12 = project_to_sphere(-phi, 0, -1, run_information.radius);

  icos_panels[0].vertex_1 = point_1, icos_panels[0].vertex_2 = point_2, icos_panels[0].vertex_3 = point_9;
  icos_panels[1].vertex_1 = point_1, icos_panels[1].vertex_2 = point_2, icos_panels[1].vertex_3 = point_11;
  icos_panels[2].vertex_1 = point_1, icos_panels[2].vertex_2 = point_5, icos_panels[2].vertex_3 = point_7;
  icos_panels[3].vertex_1 = point_1, icos_panels[3].vertex_2 = point_5, icos_panels[3].vertex_3 = point_9;
  icos_panels[4].vertex_1 = point_1, icos_panels[4].vertex_2 = point_7, icos_panels[4].vertex_3 = point_11;
  icos_panels[5].vertex_1 = point_2, icos_panels[5].vertex_2 = point_6, icos_panels[5].vertex_3 = point_8;
  icos_panels[6].vertex_1 = point_2, icos_panels[6].vertex_2 = point_6, icos_panels[6].vertex_3 = point_9;
  icos_panels[7].vertex_1 = point_2, icos_panels[7].vertex_2 = point_8, icos_panels[7].vertex_3 = point_11;
  icos_panels[8].vertex_1 = point_3, icos_panels[8].vertex_2 = point_4, icos_panels[8].vertex_3 = point_10;
  icos_panels[9].vertex_1 = point_3, icos_panels[9].vertex_2 = point_4, icos_panels[9].vertex_3 = point_12;
  icos_panels[10].vertex_1 = point_3, icos_panels[10].vertex_2 = point_5, icos_panels[10].vertex_3 = point_7;
  icos_panels[11].vertex_1 = point_3, icos_panels[11].vertex_2 = point_5, icos_panels[11].vertex_3 = point_10;
  icos_panels[12].vertex_1 = point_3, icos_panels[12].vertex_2 = point_7, icos_panels[12].vertex_3 = point_12;
  icos_panels[13].vertex_1 = point_4, icos_panels[13].vertex_2 = point_6, icos_panels[13].vertex_3 = point_8;
  icos_panels[14].vertex_1 = point_4, icos_panels[14].vertex_2 = point_6, icos_panels[14].vertex_3 = point_10;
  icos_panels[15].vertex_1 = point_4, icos_panels[15].vertex_2 = point_8, icos_panels[15].vertex_3 = point_12;
  icos_panels[16].vertex_1 = point_5, icos_panels[16].vertex_2 = point_9, icos_panels[16].vertex_3 = point_10;
  icos_panels[17].vertex_1 = point_6, icos_panels[17].vertex_2 = point_9, icos_panels[17].vertex_3 = point_10;
  icos_panels[18].vertex_1 = point_7, icos_panels[18].vertex_2 = point_11, icos_panels[18].vertex_3 = point_12;
  icos_panels[19].vertex_1 = point_8, icos_panels[19].vertex_2 = point_11, icos_panels[19].vertex_3 = point_12;

  for (int i = 0; i < 20; i++) {
    icos_panels[i].id = i;
  }

  double xval, yval, zval;
  std::vector<double> points (9, 0), point (3, 0);
  bool point_found;
  int status;

  for (int i = 0; i < run_information.point_count; i++) {
    // assign the points to the base icosahedron
    xval = x[i];
    yval = y[i];
    zval = z[i];
    point_found = false;
    for (int j = 0; j < 20; j++) {
      point[0] = xval;
      point[1] = yval;
      point[2] = zval;
      points[0] = icos_panels[j].vertex_1[0];
      points[1] = icos_panels[j].vertex_1[1];
      points[2] = icos_panels[j].vertex_1[2];
      points[3] = icos_panels[j].vertex_2[0];
      points[4] = icos_panels[j].vertex_2[1];
      points[5] = icos_panels[j].vertex_2[2];
      points[6] = icos_panels[j].vertex_3[0];
      points[7] = icos_panels[j].vertex_3[1];
      points[8] = icos_panels[j].vertex_3[2];

      status = linear_solve(points, point, 3, 1, 1);
      if (status > 0) {
        throw std::runtime_error("Error with barycoordinate computation in icosahedron initialize, line 77");
      }
      if ((point[0] >= 0) and (point[1] >= 0) and (point[2] >= 0)) {
        point_found = true;
        icos_panels[j].point_count++;
        icos_panels[j].points_inside.push_back(i);
      }

      if (point_found) {
        break;
      }
    }
    if (not point_found) {
      throw std::runtime_error("Not all points located");
    }
  }

  double x1, x2, x3, y1, y2, y3, z1, z2, z3, x12, x23, x31, y12, y23, y31, z12, z23, z31;
  std::vector<double> v12, v23, v31;
  int point_count;
  int start, end;
  double p_x, p_y, p_z;
  std::vector<int> points_to_assign;

  for (int i = 0; i < icos_panels.size(); i++) {
    if ((icos_panels[i].is_leaf) and (icos_panels[i].point_count > run_information.fast_sum_cluster_thresh)) {
      // refine, create 4 new sub-panels
      IcosPanel sub_panel1, sub_panel2, sub_panel3, sub_panel4;
      x1 = icos_panels[i].vertex_1[0], y1 = icos_panels[i].vertex_1[1], z1 = icos_panels[i].vertex_1[2];
      x2 = icos_panels[i].vertex_2[0], y2 = icos_panels[i].vertex_2[1], z2 = icos_panels[i].vertex_2[2];
      x3 = icos_panels[i].vertex_3[0], y3 = icos_panels[i].vertex_3[1], z3 = icos_panels[i].vertex_3[2];
      x12 = (x1 + x2) / 2, y12 = (y1 + y2) / 2, z12 = (z1 + z2) / 2;
      x23 = (x2 + x3) / 2, y23 = (y2 + y3) / 2, z23 = (z2 + z3) / 2;
      x31 = (x3 + x1) / 2, y31 = (y3 + y1) / 2, z31 = (z3 + z1) / 2;
      v12 = project_to_sphere(x12, y12, z12, run_information.radius);
      v23 = project_to_sphere(x23, y23, z23, run_information.radius);
      v31 = project_to_sphere(x31, y31, z31, run_information.radius);
      sub_panel1.vertex_1 = icos_panels[i].vertex_1, sub_panel1.vertex_2 = v31, sub_panel1.vertex_3 = v12;
      sub_panel2.vertex_1 = icos_panels[i].vertex_3, sub_panel2.vertex_2 = v23, sub_panel2.vertex_3 = v31;
      sub_panel3.vertex_1 = icos_panels[i].vertex_2, sub_panel3.vertex_2 = v12, sub_panel3.vertex_3 = v23;
      sub_panel4.vertex_1 = v12, sub_panel4.vertex_2 = v23, sub_panel4.vertex_3 = v31;
      sub_panel1.level = icos_panels[i].level + 1;
      sub_panel2.level = icos_panels[i].level + 1;
      sub_panel3.level = icos_panels[i].level + 1;
      sub_panel4.level = icos_panels[i].level + 1;
      icos_panels[i].is_leaf = false;
      start = icos_panels.size();
      sub_panel1.id = start;
      sub_panel2.id = start + 1;
      sub_panel3.id = start + 2;
      sub_panel4.id = start + 3;
      sub_panel1.parent_panel = icos_panels[i].id;
      sub_panel2.parent_panel = icos_panels[i].id;
      sub_panel3.parent_panel = icos_panels[i].id;
      sub_panel4.parent_panel = icos_panels[i].id;
      icos_panels[i].child_panel_1 = sub_panel1.id;
      icos_panels[i].child_panel_2 = sub_panel2.id;
      icos_panels[i].child_panel_3 = sub_panel3.id;
      icos_panels[i].child_panel_4 = sub_panel4.id;
      icos_panels.push_back(sub_panel1);
      icos_panels.push_back(sub_panel2);
      icos_panels.push_back(sub_panel3);
      icos_panels.push_back(sub_panel4);
      end = icos_panels.size();
      points_to_assign = icos_panels[i].points_inside;
      for (int j = 0; j < points_to_assign.size(); j++) {
        // loop over points in parent panel
        point_found = false;
        for (int k = start; k < end; k++) {
          // loop over child panels
          point[0] = x[points_to_assign[j]];
          point[1] = y[points_to_assign[j]];
          point[2] = z[points_to_assign[j]];
          points[0] = icos_panels[k].vertex_1[0];
          points[1] = icos_panels[k].vertex_1[1];
          points[2] = icos_panels[k].vertex_1[2];
          points[3] = icos_panels[k].vertex_2[0];
          points[4] = icos_panels[k].vertex_2[1];
          points[5] = icos_panels[k].vertex_2[2];
          points[6] = icos_panels[k].vertex_3[0];
          points[7] = icos_panels[k].vertex_3[1];
          points[8] = icos_panels[k].vertex_3[2];

          status = linear_solve(points, point, 3, 1, 1);
          if (status > 0) {
            throw std::runtime_error("Error with barycoordinate computation in icosahedron initialize, line 162");
          }

          if ((point[0] >= 0) and (point[1] >= 0) and (point[2] >= 0)) {
            point_found = true;
            icos_panels[k].point_count++;
            icos_panels[k].points_inside.push_back(points_to_assign[j]);
          }

          if (point_found) {
            break;
          }
        }
        if (not point_found) {
          throw std::runtime_error("Point not correctly found");
        }
      }
      point_count = icos_panels[start].point_count + icos_panels[start+1].point_count + icos_panels[start+2].point_count + icos_panels[start+3].point_count;
      if (point_count != icos_panels[i].point_count) {
        std::cout << "discrepancy: " << point_count << " vs " << icos_panels[i].point_count << std::endl;
        throw std::runtime_error("Not all points assigned");
      }
    }
  }

  double xc, yc, zc, d1, d2, d3, d;
  std::vector<double> vc;

  for (int i = 0; i < icos_panels.size(); i++) {
    x1 = icos_panels[i].vertex_1[0], y1 = icos_panels[i].vertex_1[1], z1 = icos_panels[i].vertex_1[2];
    x2 = icos_panels[i].vertex_2[0], y2 = icos_panels[i].vertex_2[1], z2 = icos_panels[i].vertex_2[2];
    x3 = icos_panels[i].vertex_3[0], y3 = icos_panels[i].vertex_3[1], z3 = icos_panels[i].vertex_3[2];
    xc = (x1 + x2 + x3) / 3.0;
    yc = (y1 + y2 + y3) / 3.0;
    zc = (z1 + z2 + z3) / 3.0;
    vc = project_to_sphere(xc, yc, zc, run_information.radius);
    xc = vc[0], yc = vc[1], zc = vc[2];
    icos_panels[i].center_p[0] = xc, icos_panels[i].center_p[1] = yc, icos_panels[i].center_p[2] = zc;
    d1 = pow(x1 - xc, 2) + pow(y1 - yc, 2) + pow(z1 - zc, 2);
    d2 = pow(x2 - xc, 2) + pow(y2 - yc, 2) + pow(z2 - zc, 2);
    d3 = pow(x3 - xc, 2) + pow(y3 - yc, 2) + pow(z3 - zc, 2);
    d = std::max(d1, std::max(d2, d3));
    icos_panels[i].radius = sqrt(d);
    if (icos_panels[i].point_count != icos_panels[i].points_inside.size()) {
      throw std::runtime_error("Error in icosahedron initialization, line 207");
    }
  }
}

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
      for (int i = 0; i < points_to_assign.size(); i++) {
        xval = x[points_to_assign[i]];
        yval = y[points_to_assign[i]];
        zval = z[points_to_assign[i]];
        xieta = xieta_from_xyz(xval, yval, zval, cube_panels[i].face);
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
        cube_panels[start+which_panel].points_inside.push_back(points_to_assign[i]);
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

void dual_tree_traversal_icos(const RunConfig& run_information, std::vector<InteractPair>& interactions, const std::vector<IcosPanel>& icos_panels) {
  // performs dual tree traversal to find {P,C} - {P,C} interactions
  std::queue<int> target_triangles;
  std::queue<int> source_triangles;

  for (int i = 0; i < 20; i++) {
    for (int j = 0; j < 20; j++) {
      target_triangles.push(i);
      source_triangles.push(j);
    }
  }

  int index_target, index_source;
  double x1, x2, y1, y2, z1, z2, dist, separation;

  while (target_triangles.size() > 0) {
    index_target = target_triangles.front();
    index_source = source_triangles.front();
    target_triangles.pop();
    source_triangles.pop();
    if (icos_panels[index_target].point_count == 0) continue;
    if (icos_panels[index_source].point_count == 0) continue;
    x1 = icos_panels[index_target].center_p[0], y1 = icos_panels[index_target].center_p[1], z1 = icos_panels[index_target].center_p[2];
    x2 = icos_panels[index_source].center_p[0], y2 = icos_panels[index_source].center_p[1], z2 = icos_panels[index_source].center_p[2];
    dist = gcdist(x1, y1, z1, x2, y2, z2, run_information.radius);
    separation = (icos_panels[index_target].radius + icos_panels[index_source].radius) / dist;
    if ((dist > 0) and (separation < run_information.fast_sum_theta)) {
      // triangles are well separated
      InteractPair new_interact = {index_target, index_source, 0};
      if (icos_panels[index_target].point_count >= run_information.fast_sum_cluster_thresh) {
        new_interact.interact_type += 2;
      }
      if (icos_panels[index_source].point_count >= run_information.fast_sum_cluster_thresh) {
        new_interact.interact_type += 1;
      }
      interactions.push_back(new_interact);
    } else {
      // not well separated
      if ((icos_panels[index_target].point_count < run_information.fast_sum_cluster_thresh) and (icos_panels[index_source].point_count < run_information.fast_sum_cluster_thresh)) {
        // both have few points
        InteractPair new_interact = {index_target, index_source, 0};
        interactions.push_back(new_interact);
      } else if (icos_panels[index_target].is_leaf and icos_panels[index_source].is_leaf) {
        // both are leaves
        InteractPair new_interact = {index_target, index_source, 0};
        interactions.push_back(new_interact);
      } else if (icos_panels[index_target].is_leaf) {
        // target is leaf, break apart source
        target_triangles.push(index_target);
        target_triangles.push(index_target);
        target_triangles.push(index_target);
        target_triangles.push(index_target);
        source_triangles.push(icos_panels[index_source].child_panel_1);
        source_triangles.push(icos_panels[index_source].child_panel_2);
        source_triangles.push(icos_panels[index_source].child_panel_3);
        source_triangles.push(icos_panels[index_source].child_panel_4);
      } else if (icos_panels[index_source].is_leaf) {
        // source is leaf, break apart target
        source_triangles.push(index_source);
        source_triangles.push(index_source);
        source_triangles.push(index_source);
        source_triangles.push(index_source);
        target_triangles.push(icos_panels[index_target].child_panel_1);
        target_triangles.push(icos_panels[index_target].child_panel_2);
        target_triangles.push(icos_panels[index_target].child_panel_3);
        target_triangles.push(icos_panels[index_target].child_panel_4);
      } else {
        // neither is leaf, break apart triangle with more points
        if (icos_panels[index_target].point_count >= icos_panels[index_source].point_count) {
          // target has more points, break apart target
          source_triangles.push(index_source);
          source_triangles.push(index_source);
          source_triangles.push(index_source);
          source_triangles.push(index_source);
          target_triangles.push(icos_panels[index_target].child_panel_1);
          target_triangles.push(icos_panels[index_target].child_panel_2);
          target_triangles.push(icos_panels[index_target].child_panel_3);
          target_triangles.push(icos_panels[index_target].child_panel_4);
        } else {
          // source has more points, break apart source
          target_triangles.push(index_target);
          target_triangles.push(index_target);
          target_triangles.push(index_target);
          target_triangles.push(index_target);
          source_triangles.push(icos_panels[index_source].child_panel_1);
          source_triangles.push(icos_panels[index_source].child_panel_2);
          source_triangles.push(icos_panels[index_source].child_panel_3);
          source_triangles.push(icos_panels[index_source].child_panel_4);
        }
      }
    }
  }
}

void dual_tree_traversal_cube(const RunConfig& run_information, std::vector<InteractPair>& interactions, const std::vector<CubePanel>& cube_panels) {
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
