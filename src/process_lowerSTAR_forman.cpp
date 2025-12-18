#include <Rcpp.h>
#include <vector>
#include <string>
#include <algorithm>
#include <sstream>
#include <unordered_set>

using namespace Rcpp;

// C++11 COMPATIBLE VERSION
// [[Rcpp::export]]
List process_lowerSTAR_forman_cpp(List list_lowerSTAR, IntegerVector vertex_ids) {
  int n = list_lowerSTAR.size();
  
  std::vector<std::string> VE_all;
  std::vector<std::string> CR_all;
  
  for (int i = 0; i < n; i++) {
    int current_vertex_id = vertex_ids[i];
    
    DataFrame df = as<DataFrame>(list_lowerSTAR[i]);
    int m = df.nrows();
    
    if (m == 0) {
      CR_all.push_back(std::to_string(current_vertex_id));
      continue;
    }
    
    CharacterVector simplex_ids = df["lexi_id"];
    CharacterVector labels = df["lexi_label"];
    
    // Create vector of pairs for sorting (C++11 compatible)
    std::vector<std::pair<double, std::string> > sorted_pairs;
    
    for (int j = 0; j < m; j++) {
      std::string simplex_id = as<std::string>(simplex_ids[j]);
      std::string label = as<std::string>(labels[j]);
      
      // Extract first value from label
      double first_val = 0.0;
      std::istringstream iss(label);
      std::string token;
      if (iss >> token) {
        try {
          first_val = std::stod(token);
        } catch (...) {
          first_val = 0.0;
        }
      }
      
      sorted_pairs.push_back(std::make_pair(first_val, simplex_id));
    }
    
    // Sort by first value (C++11 compatible lambda)
    std::sort(sorted_pairs.begin(), sorted_pairs.end(),
              [](const std::pair<double, std::string>& a, 
                 const std::pair<double, std::string>& b) {
                return a.first < b.first;
              });
    
    // Extract sorted simplex IDs
    std::vector<std::string> sorted_ids;
    for (size_t j = 0; j < sorted_pairs.size(); j++) {
      sorted_ids.push_back(sorted_pairs[j].second);
    }
    
    // Track paired simplices
    std::vector<bool> paired(m, false);
    bool vertex_paired = false;
    
    // ----- RULE 1: Pair CURRENT VERTEX with an edge -----
    for (int j = 0; j < m; j++) {
      std::string simplex_id = sorted_ids[j];
      
      // Count spaces to determine vertices
      int spaces = 0;
      for (size_t k = 0; k < simplex_id.size(); k++) {
        if (simplex_id[k] == ' ') spaces++;
      }
      
      // Edge has 2 vertices (1 space)
      if (spaces == 1) {
        // Check if edge contains current vertex
        std::istringstream iss(simplex_id);
        std::string token;
        bool contains_current = false;
        
        while (iss >> token) {
          if (std::stoi(token) == current_vertex_id) {
            contains_current = true;
            break;
          }
        }
        
        if (contains_current) {
          VE_all.push_back(std::to_string(current_vertex_id) + ":" + simplex_id);
          paired[j] = true;
          vertex_paired = true;
          break;
        }
      }
    }
    
    if (!vertex_paired) {
      CR_all.push_back(std::to_string(current_vertex_id));
    }
    
    // ----- RULE 2: Pair edges with faces -----
    for (int j = 0; j < m; j++) {
      if (paired[j]) continue;
      
      std::string simplex_id = sorted_ids[j];
      int spaces = 0;
      for (size_t k = 0; k < simplex_id.size(); k++) {
        if (simplex_id[k] == ' ') spaces++;
      }
      
      // Edge (2 vertices = 1 space)
      if (spaces == 1) {
        bool edge_paired = false;
        
        // Extract edge vertices
        std::vector<int> edge_verts;
        std::istringstream iss_edge(simplex_id);
        std::string token;
        while (iss_edge >> token) {
          edge_verts.push_back(std::stoi(token));
        }
        
        // Look for a face that contains this edge
        for (int k = j + 1; k < m; k++) {
          if (paired[k]) continue;
          
          std::string face_id = sorted_ids[k];
          int face_spaces = 0;
          for (size_t l = 0; l < face_id.size(); l++) {
            if (face_id[l] == ' ') face_spaces++;
          }
          
          // Face (3 vertices = 2 spaces)
          if (face_spaces == 2) {
            // Extract face vertices
            std::vector<int> face_verts;
            std::istringstream iss_face(face_id);
            while (iss_face >> token) {
              face_verts.push_back(std::stoi(token));
            }
            
            // Check if face contains all edge vertices
            bool contains_edge = true;
            for (size_t ev = 0; ev < edge_verts.size(); ev++) {
              bool found = false;
              for (size_t fv = 0; fv < face_verts.size(); fv++) {
                if (edge_verts[ev] == face_verts[fv]) {
                  found = true;
                  break;
                }
              }
              if (!found) {
                contains_edge = false;
                break;
              }
            }
            
            if (contains_edge) {
              VE_all.push_back(simplex_id + ":" + face_id);
              paired[j] = true;
              paired[k] = true;
              edge_paired = true;
              break;
            }
          }
        }
        
        if (!edge_paired) {
          CR_all.push_back(simplex_id);
        }
      }
    }
    
    // ----- RULE 3: Unpaired faces are critical -----
    for (int j = 0; j < m; j++) {
      if (paired[j]) continue;
      
      std::string simplex_id = sorted_ids[j];
      int spaces = 0;
      for (size_t k = 0; k < simplex_id.size(); k++) {
        if (simplex_id[k] == ' ') spaces++;
      }
      
      // Face (3 vertices = 2 spaces)
      if (spaces == 2) {
        CR_all.push_back(simplex_id);
      }
    }
  }
  
  // Remove duplicates
  std::unordered_set<std::string> cr_set;
  for (size_t i = 0; i < CR_all.size(); i++) {
    cr_set.insert(CR_all[i]);
  }
  
  CR_all.clear();
  CR_all.insert(CR_all.end(), cr_set.begin(), cr_set.end());
  
  // Convert to R vectors
  CharacterVector VE_final(VE_all.size());
  CharacterVector CR_final(CR_all.size());
  
  for (size_t i = 0; i < VE_all.size(); i++) {
    VE_final[i] = VE_all[i];
  }
  
  for (size_t i = 0; i < CR_all.size(); i++) {
    CR_final[i] = CR_all[i];
  }
  
  // Statistics
  int ve_count = 0;
  for (size_t i = 0; i < VE_all.size(); i++) {
    size_t colon_pos = VE_all[i].find(':');
    std::string first_part = VE_all[i].substr(0, colon_pos);
    int spaces = 0;
    for (size_t j = 0; j < first_part.size(); j++) {
      if (first_part[j] == ' ') spaces++;
    }
    if (spaces == 0) ve_count++;
  }
  
  Rcout << "Forman gradient (C++11):\n";
  Rcout << "  Vertex->Edge pairs: " << ve_count << "\n";
  Rcout << "  Total gradient pairs: " << VE_final.size() << "\n";
  Rcout << "  Critical simplices: " << CR_final.size() << "\n";
  
  return List::create(
    Named("VE_") = VE_final,
    Named("CR_") = CR_final
  );
}