#include <Rcpp.h>
#include <unordered_set>
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
  
  // Pre-process: Create lookup sets for edges and faces
  std::vector<std::unordered_set<std::string>> edge_sets(n_edge);
  std::vector<std::unordered_set<std::string>> face_sets(n_face);
  
  // Parse edge vertices into sets (silent)
  for (int j = 0; j < n_edge; j++) {
    std::string edge_str = as<std::string>(edge_lexi_id[j]);
    std::istringstream iss(edge_str);
    std::string token;
    while (iss >> token) {
      edge_sets[j].insert(token);
    }
  }
  
  // Parse face vertices into sets (silent)  
  for (int j = 0; j < n_face; j++) {
    std::string face_str = as<std::string>(face_lexi_id[j]);
    std::istringstream iss(face_str);
    std::string token;
    while (iss >> token) {
      face_sets[j].insert(token);
    }
  }
  
  // Output vectors
  CharacterVector vertTO_lexi_label;
  CharacterVector vertTO_lexi_id;
  
  // Main processing with O(1) lookups (silent)
  for (int i = 0; i < n_vertex; i++) {
    std::string v_id = as<std::string>(vertex_i123[i]);
    
    // Check edges - O(1) lookup per edge!
    for (int j = 0; j < n_edge; j++) {
      if (edge_sets[j].find(v_id) != edge_sets[j].end()) {
        vertTO_lexi_id.push_back(edge_lexi_id[j]);
        vertTO_lexi_label.push_back(edge_lexi_label[j]);
      }
    }
    
    // Check faces - O(1) lookup per face!
    for (int j = 0; j < n_face; j++) {
      if (face_sets[j].find(v_id) != face_sets[j].end()) {
        vertTO_lexi_id.push_back(face_lexi_id[j]);
        vertTO_lexi_label.push_back(face_lexi_label[j]);
      }
    }
  }
  
  return DataFrame::create(
    Named("lexi_label") = vertTO_lexi_label, 
    Named("lexi_id") = vertTO_lexi_id
  );
}