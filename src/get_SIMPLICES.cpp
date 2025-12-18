#include <Rcpp.h>
#include <vector>
#include <string>
#include <algorithm>
#include <unordered_map>

using namespace Rcpp;

// Fast vertex structure
struct Vertex {
  double x, y, z;
  int id;
  std::string z_str;
};

// Fast comparison for sorting
bool compareByZ(const std::pair<double, int>& a, const std::pair<double, int>& b) {
  return a.first < b.first;
}

// [[Rcpp::export]]
List get_SIMPLICES_fast_cpp(NumericMatrix vertices_mat, 
                            NumericMatrix faces_mat, 
                            NumericMatrix edges_mat,
                            IntegerVector input_truth) {
  
  const int n_vertices = vertices_mat.nrow();
  const int n_faces = faces_mat.nrow();
  const int n_edges = edges_mat.nrow();
  
  // Pre-allocate output vectors
  NumericVector v_x(n_vertices), v_y(n_vertices);
  CharacterVector v_z(n_vertices);
  IntegerVector v_id(n_vertices);
  
  // Process vertices
  std::vector<Vertex> vertices(n_vertices);
  for (int i = 0; i < n_vertices; ++i) {
    vertices[i].x = vertices_mat(i, 0);
    vertices[i].y = vertices_mat(i, 1);
    vertices[i].z = vertices_mat(i, 2);
    vertices[i].id = input_truth[i];
    vertices[i].z_str = std::to_string(vertices[i].z);
    
    v_x[i] = vertices[i].x;
    v_y[i] = vertices[i].y;
    v_z[i] = vertices[i].z_str;
    v_id[i] = vertices[i].id;
  }
  
  // Create lookup: row index -> vertex (for edges/faces)
  // No need for hashmap since edges/faces use 1-based row indices
  
  // Process edges
  IntegerVector e_i1(n_edges), e_i2(n_edges);
  NumericVector e_ii1X(n_edges), e_ii1Y(n_edges), e_ii2X(n_edges), e_ii2Y(n_edges);
  CharacterVector e_ii1Z(n_edges), e_ii2Z(n_edges);
  CharacterVector e_label(n_edges), e_idlabel(n_edges);
  CharacterVector e_lexi_label(n_edges), e_lexi_id(n_edges);
  
  for (int i = 0; i < n_edges; ++i) {
    int idx1 = edges_mat(i, 0) - 1;  // Convert to 0-based
    int idx2 = edges_mat(i, 1) - 1;
    
    const Vertex& v1 = vertices[idx1];
    const Vertex& v2 = vertices[idx2];
    
    // Store basic info
    e_i1[i] = v1.id;
    e_i2[i] = v2.id;
    
    e_ii1X[i] = v1.x;
    e_ii1Y[i] = v1.y;
    e_ii1Z[i] = v1.z_str;
    
    e_ii2X[i] = v2.x;
    e_ii2Y[i] = v2.y;
    e_ii2Z[i] = v2.z_str;
    
    // Create label and idlabel
    std::string label = v1.z_str + " " + v2.z_str;
    std::string idlabel = std::to_string(v1.id) + " " + std::to_string(v2.id);
    
    // Sort for lexicographic order
    if (v1.z > v2.z) {
      e_lexi_label[i] = v2.z_str + " " + v1.z_str;
      e_lexi_id[i] = std::to_string(v2.id) + " " + std::to_string(v1.id);
    } else {
      e_lexi_label[i] = v1.z_str + " " + v2.z_str;
      e_lexi_id[i] = std::to_string(v1.id) + " " + std::to_string(v2.id);
    }
    
    e_label[i] = label;
    e_idlabel[i] = idlabel;
  }
  
  // Process faces
  IntegerVector f_i1(n_faces), f_i2(n_faces), f_i3(n_faces);
  NumericVector f_ii1X(n_faces), f_ii1Y(n_faces), 
  f_ii2X(n_faces), f_ii2Y(n_faces),
  f_ii3X(n_faces), f_ii3Y(n_faces);
  CharacterVector f_ii1Z(n_faces), f_ii2Z(n_faces), f_ii3Z(n_faces);
  CharacterVector f_label(n_faces), f_idlabel(n_faces);
  CharacterVector f_lexi_label(n_faces), f_lexi_id(n_faces);
  
  // Temporary arrays for sorting
  std::vector<std::pair<double, int>> z_ids(3);
  
  for (int i = 0; i < n_faces; ++i) {
    int idx1 = faces_mat(i, 0) - 1;
    int idx2 = faces_mat(i, 1) - 1;
    int idx3 = faces_mat(i, 2) - 1;
    
    const Vertex& v1 = vertices[idx1];
    const Vertex& v2 = vertices[idx2];
    const Vertex& v3 = vertices[idx3];
    
    // Store basic info
    f_i1[i] = v1.id;
    f_i2[i] = v2.id;
    f_i3[i] = v3.id;
    
    f_ii1X[i] = v1.x;
    f_ii1Y[i] = v1.y;
    f_ii1Z[i] = v1.z_str;
    
    f_ii2X[i] = v2.x;
    f_ii2Y[i] = v2.y;
    f_ii2Z[i] = v2.z_str;
    
    f_ii3X[i] = v3.x;
    f_ii3Y[i] = v3.y;
    f_ii3Z[i] = v3.z_str;
    
    // Create label and idlabel
    std::string label = v1.z_str + " " + v2.z_str + " " + v3.z_str;
    std::string idlabel = std::to_string(v1.id) + " " + 
      std::to_string(v2.id) + " " + 
      std::to_string(v3.id);
    
    // Sort for lexicographic order
    z_ids[0] = {v1.z, v1.id};
    z_ids[1] = {v2.z, v2.id};
    z_ids[2] = {v3.z, v3.id};
    
    std::sort(z_ids.begin(), z_ids.end(), compareByZ);
    
    f_lexi_label[i] = std::to_string(z_ids[0].first) + " " +
      std::to_string(z_ids[1].first) + " " +
      std::to_string(z_ids[2].first);
    f_lexi_id[i] = std::to_string(z_ids[0].second) + " " +
      std::to_string(z_ids[1].second) + " " +
      std::to_string(z_ids[2].second);
    
    f_label[i] = label;
    f_idlabel[i] = idlabel;
  }
  
  // Create data frames
  DataFrame vertices_df = DataFrame::create(
    _["X"] = v_x,
    _["Y"] = v_y,
    _["Z"] = v_z,
    _["i123"] = v_id
  );
  
  DataFrame edges_df = DataFrame::create(
    _["i1"] = e_i1,
    _["i2"] = e_i2,
    _["ii1X"] = e_ii1X,
    _["ii1Y"] = e_ii1Y,
    _["ii1Z"] = e_ii1Z,
    _["ii2X"] = e_ii2X,
    _["ii2Y"] = e_ii2Y,
    _["ii2Z"] = e_ii2Z,
    _["label"] = e_label,
    _["idlabel"] = e_idlabel,
    _["lexi_label"] = e_lexi_label,
    _["lexi_id"] = e_lexi_id
  );
  
  DataFrame faces_df = DataFrame::create(
    _["i1"] = f_i1,
    _["i2"] = f_i2,
    _["i3"] = f_i3,
    _["ii1X"] = f_ii1X,
    _["ii1Y"] = f_ii1Y,
    _["ii1Z"] = f_ii1Z,
    _["ii2X"] = f_ii2X,
    _["ii2Y"] = f_ii2Y,
    _["ii2Z"] = f_ii2Z,
    _["ii3X"] = f_ii3X,
    _["ii3Y"] = f_ii3Y,
    _["ii3Z"] = f_ii3Z,
    _["label"] = f_label,
    _["idlabel"] = f_idlabel,
    _["lexi_label"] = f_lexi_label,
    _["lexi_id"] = f_lexi_id
  );
  
  return List::create(
    _["vertices"] = vertices_df,
    _["edges"] = edges_df,
    _["faces"] = faces_df
  );
}