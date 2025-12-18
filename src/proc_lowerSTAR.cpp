#include <Rcpp.h>
#include <vector>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <sstream>

using namespace Rcpp;

// Pre-process simplex data into efficient structures
struct SimplexInfo {
  int dimension;
  std::vector<int> vertices;
  std::string id;
};

SimplexInfo parseSimplex(const std::string& simplex_id) {
  SimplexInfo info;
  info.id = simplex_id;
  
  std::istringstream iss(simplex_id);
  std::string token;
  while (iss >> token) {
    info.vertices.push_back(std::stoi(token));
  }
  info.dimension = info.vertices.size() - 1;
  
  // Sort vertices for consistent comparison
  std::sort(info.vertices.begin(), info.vertices.end());
  
  return info;
}

// Fast face/coface relationship checking
bool isFaceOf(const std::vector<int>& face, const std::vector<int>& coface) {
  // Check if all vertices of face are in coface (both are sorted)
  auto face_it = face.begin();
  auto coface_it = coface.begin();
  
  while (face_it != face.end() && coface_it != coface.end()) {
    if (*face_it < *coface_it) {
      return false; // Face vertex not found in coface
    } else if (*face_it > *coface_it) {
      ++coface_it;
    } else {
      ++face_it;
      ++coface_it;
    }
  }
  
  return face_it == face.end(); // All face vertices found
}

// [[Rcpp::export]]
List process_lowerSTAR_cpp(List list_lowerSTAR, DataFrame vertex) {
  int n = list_lowerSTAR.size();
  std::unordered_set<std::string> VE_set, CR_set;
  
  IntegerVector vertex_ids = vertex["i123"];
  
  for (int i = 0; i < n; i++) {
    DataFrame df = as<DataFrame>(list_lowerSTAR[i]);
    if (df.nrows() == 0) {
      CR_set.insert(std::to_string(vertex_ids[i]));
      continue;
    }
    
    CharacterVector lexi_id = df["lexi_id"];
    
    // Parse all simplices upfront
    std::vector<SimplexInfo> simplices;
    simplices.reserve(lexi_id.size());
    
    for (int j = 0; j < lexi_id.size(); j++) {
      simplices.push_back(parseSimplex(as<std::string>(lexi_id[j])));
    }
    
    // Group simplices by dimension
    std::vector<SimplexInfo> vertices, edges, faces;
    for (const auto& simplex : simplices) {
      if (simplex.dimension == 0) {
        vertices.push_back(simplex);
      } else if (simplex.dimension == 1) {
        edges.push_back(simplex);
      } else if (simplex.dimension == 2) {
        faces.push_back(simplex);
      }
    }
    
    // First: Pair vertex with first simplex of next dimension
    // This handles the vertex-edge and vertex-face pairings
    if (!edges.empty()) {
      // Pair vertex with first edge
      VE_set.insert(std::to_string(vertex_ids[i]) + ":" + edges[0].id);
    } else if (!faces.empty()) {
      // If no edges, pair vertex with first face
      VE_set.insert(std::to_string(vertex_ids[i]) + ":" + faces[0].id);
    } else if (!vertices.empty()) {
      // If only vertices (shouldn't happen in proper filtration)
      VE_set.insert(std::to_string(vertex_ids[i]) + ":" + vertices[0].id);
    }
    
    // Now handle edge-face pairings
    if (edges.size() > 1 || faces.size() > 0) {
      // Build edge-face adjacency for O(1) lookups
      std::unordered_map<std::string, std::vector<std::string>> edge_to_faces;
      for (const auto& edge : edges) {
        for (const auto& face : faces) {
          if (isFaceOf(edge.vertices, face.vertices)) {
            edge_to_faces[edge.id].push_back(face.id);
          }
        }
      }
      
      // Greedy Morse pairing for edge-face relationships
      std::unordered_set<std::string> paired_faces;
      std::vector<std::string> VE_local, CR_local;
      
      // Start from j=1 because edges[0] is already paired with vertex
      for (size_t j = 1; j < edges.size(); j++) {
        const auto& edge = edges[j];
        auto face_it = edge_to_faces.find(edge.id);
        
        if (face_it != edge_to_faces.end() && !face_it->second.empty()) {
          // Find first unpaired face
          for (const auto& face_id : face_it->second) {
            if (paired_faces.find(face_id) == paired_faces.end()) {
              VE_local.push_back(edge.id + ":" + face_id);
              paired_faces.insert(face_id);
              break;
            }
          }
        }
        
        // If no face was paired with this edge, it's critical
        bool edge_paired = false;
        for (const auto& pair : VE_local) {
          if (pair.find(edge.id + ":") == 0) {
            edge_paired = true;
            break;
          }
        }
        // Also check if this edge was the first one paired with vertex
        if (j == 0) edge_paired = true;
        
        if (!edge_paired) {
          CR_local.push_back(edge.id);
        }
      }
      
      // Add unpaired faces to critical list
      for (const auto& face : faces) {
        if (paired_faces.find(face.id) == paired_faces.end()) {
          CR_local.push_back(face.id);
        }
      }
      
      // Also handle the case where we have only edges (no faces)
      if (faces.empty() && edges.size() > 1) {
        // All edges after the first one are critical
        for (size_t j = 1; j < edges.size(); j++) {
          CR_local.push_back(edges[j].id);
        }
      }
      
      // Update global sets
      for (const auto& pair : VE_local) VE_set.insert(pair);
      for (const auto& crit : CR_local) CR_set.insert(crit);
    }
    
    // Handle any remaining vertices (should be minimal in proper filtration)
    for (size_t j = 1; j < vertices.size(); j++) {
      CR_set.insert(vertices[j].id);
    }
  }
  
  // Convert to output
  std::vector<std::string> VE_final(VE_set.begin(), VE_set.end());
  std::vector<std::string> CR_final(CR_set.begin(), CR_set.end());
  
  return List::create(_["VE_"] = VE_final, _["CR_"] = CR_final);
}

// [[Rcpp::export]]
List proc_lowerSTAR_cpp(List list_lowerSTAR, DataFrame vertex) {
    int n = list_lowerSTAR.size();
    std::unordered_set<std::string> VE_set, CR_set;
    
    IntegerVector vertex_ids = vertex["i123"];
    
    for (int i = 0; i < n; i++) {
        DataFrame df = as<DataFrame>(list_lowerSTAR[i]);
        if (df.nrows() == 0) {
            CR_set.insert(std::to_string(vertex_ids[i]));
            continue;
        }
        
        CharacterVector lexi_id = df["lexi_id"];
        
        // Parse all simplices upfront
        std::vector<SimplexInfo> simplices;
        simplices.reserve(lexi_id.size());
        
        for (int j = 0; j < lexi_id.size(); j++) {
            simplices.push_back(parseSimplex(as<std::string>(lexi_id[j])));
        }
        
        // First pairing: vertex with first simplex
        if (!simplices.empty()) {
            VE_set.insert(std::to_string(vertex_ids[i]) + ":" + simplices[0].id);
        }
        
        if (simplices.size() <= 1) continue;
        
        // Separate faces (dimension 1) and cofaces (dimension 2)
        std::vector<SimplexInfo> faces, cofaces;
        for (size_t j = 1; j < simplices.size(); j++) {
            if (simplices[j].dimension == 1) {
                faces.push_back(simplices[j]);
            } else {
                cofaces.push_back(simplices[j]);
            }
        }
        
        if (faces.empty()) continue;
        
        // Build face-coface adjacency for O(1) lookups
        std::unordered_map<std::string, std::vector<std::string>> face_to_cofaces;
        for (const auto& face : faces) {
            for (const auto& coface : cofaces) {
                if (isFaceOf(face.vertices, coface.vertices)) {
                    face_to_cofaces[face.id].push_back(coface.id);
                }
            }
        }
        
        // Greedy Morse pairing
        std::unordered_set<std::string> paired_cofaces;
        std::vector<std::string> VE_local, CR_local;
        
        for (const auto& face : faces) {
            auto coface_it = face_to_cofaces.find(face.id);
            if (coface_it != face_to_cofaces.end() && !coface_it->second.empty()) {
                // Find first unpaired coface
                for (const auto& coface_id : coface_it->second) {
                    if (paired_cofaces.find(coface_id) == paired_cofaces.end()) {
                        VE_local.push_back(face.id + ":" + coface_id);
                        paired_cofaces.insert(coface_id);
                        break;
                    }
                }
            }
            
            // If no coface was paired with this face, it's critical
            bool face_paired = false;
            for (const auto& pair : VE_local) {
                if (pair.find(face.id + ":") == 0) {
                    face_paired = true;
                    break;
                }
            }
            if (!face_paired) {
                CR_local.push_back(face.id);
            }
        }
        
        // Add unpaired cofaces to critical list
        for (const auto& coface : cofaces) {
            if (paired_cofaces.find(coface.id) == paired_cofaces.end()) {
                CR_local.push_back(coface.id);
            }
        }
        
        // Update global sets
        for (const auto& pair : VE_local) VE_set.insert(pair);
        for (const auto& crit : CR_local) CR_set.insert(crit);
    }
    
    // Convert to output
    std::vector<std::string> VE_final(VE_set.begin(), VE_set.end());
    std::vector<std::string> CR_final(CR_set.begin(), CR_set.end());
    
    return List::create(_["VE_"] = VE_final, _["CR_"] = CR_final);
}