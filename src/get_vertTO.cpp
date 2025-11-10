#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <sstream>
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame get_vertTO_cpp(DataFrame vertex, DataFrame edge, DataFrame face) {
  int n_vertex = vertex.nrow();
  int n_edge = edge.nrow();
  int n_face = face.nrow();
  
  // Get the data columns once
  CharacterVector vertex_i123 = vertex["i123"];
  CharacterVector edge_lexi_id = edge["lexi_id"];
  CharacterVector edge_lexi_label = edge["lexi_label"];
  CharacterVector face_lexi_id = face["lexi_id"];
  CharacterVector face_lexi_label = face["lexi_label"];
  
  // Create lookup map: vertex_id -> vector of edge/face indices
  std::unordered_map<std::string, std::vector<int>> vertex_to_edges;
  std::unordered_map<std::string, std::vector<int>> vertex_to_faces;
  
  // Pre-process edges: Build lookup in O(n_edge)
  for (int j = 0; j < n_edge; j++) {
    std::string edge_str = as<std::string>(edge_lexi_id[j]);
    std::istringstream iss(edge_str);
    std::string token;
    while (iss >> token) {
      vertex_to_edges[token].push_back(j);
    }
  }
  
  // Pre-process faces: Build lookup in O(n_face)  
  for (int j = 0; j < n_face; j++) {
    std::string face_str = as<std::string>(face_lexi_id[j]);
    std::istringstream iss(face_str);
    std::string token;
    while (iss >> token) {
      vertex_to_faces[token].push_back(j);
    }
  }
  
  // Output vectors - use std::vector for efficient growth
  std::vector<std::string> vertTO_lexi_label;
  std::vector<std::string> vertTO_lexi_id;
  
  // Reserve space for efficiency
  vertTO_lexi_label.reserve(n_vertex * 12);
  vertTO_lexi_id.reserve(n_vertex * 12);
  
  // Main processing: O(n_vertex) with direct lookups!
  for (int i = 0; i < n_vertex; i++) {
    std::string v_id = as<std::string>(vertex_i123[i]);
    
    // Lookup edges for this vertex - O(1)
    auto edge_it = vertex_to_edges.find(v_id);
    if (edge_it != vertex_to_edges.end()) {
      for (int edge_idx : edge_it->second) {
        vertTO_lexi_id.push_back(as<std::string>(edge_lexi_id[edge_idx]));
        vertTO_lexi_label.push_back(as<std::string>(edge_lexi_label[edge_idx]));
      }
    }
    
    // Lookup faces for this vertex - O(1)
    auto face_it = vertex_to_faces.find(v_id);
    if (face_it != vertex_to_faces.end()) {
      for (int face_idx : face_it->second) {
        vertTO_lexi_id.push_back(as<std::string>(face_lexi_id[face_idx]));
        vertTO_lexi_label.push_back(as<std::string>(face_lexi_label[face_idx]));
      }
    }
  }
  
  // Convert back to Rcpp vectors
  return DataFrame::create(
    Named("lexi_label") = CharacterVector(vertTO_lexi_label.begin(), vertTO_lexi_label.end()),
    Named("lexi_id") = CharacterVector(vertTO_lexi_id.begin(), vertTO_lexi_id.end())
  );
}