#include <Rcpp.h>
#include <unordered_set>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]
List get_MESH_cpp(NumericMatrix vertices, IntegerMatrix faces) {
  int n_faces = faces.nrow();
  
  // Fast edge creation using hash set
  std::unordered_set<uint64_t> edge_set;
  std::vector<std::pair<int, int>> edges;
  
  auto make_key = [](int a, int b) {
    return ((uint64_t)std::min(a, b) << 32) | (uint64_t)std::max(a, b);
  };
  
  for (int i = 0; i < n_faces; i++) {
    int v1 = faces(i, 0);
    int v2 = faces(i, 1);
    int v3 = faces(i, 2);
    
    uint64_t key12 = make_key(v1, v2);
    if (edge_set.insert(key12).second) {
      edges.push_back(std::make_pair(v1, v2));
    }
    
    uint64_t key23 = make_key(v2, v3);
    if (edge_set.insert(key23).second) {
      edges.push_back(std::make_pair(v2, v3));
    }
    
    uint64_t key31 = make_key(v3, v1);
    if (edge_set.insert(key31).second) {
      edges.push_back(std::make_pair(v3, v1));
    }
  }
  
  // Convert to data frame with only essential columns
  IntegerMatrix edges_mat(edges.size(), 2);
  for (size_t i = 0; i < edges.size(); i++) {
    edges_mat(i, 0) = edges[i].first;
    edges_mat(i, 1) = edges[i].second;
  }
  
  List edgesDF = List::create(
    Named("i1") = edges_mat(_, 0),
    Named("i2") = edges_mat(_, 1)
    // Skip length, angle, exterior, coplanar - not used!
  );
  
  return List::create(
    Named("vertices") = vertices,
    Named("faces") = faces,
    Named("edgesDF") = edgesDF
  );
}