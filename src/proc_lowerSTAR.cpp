#include <Rcpp.h>
#include <vector>
#include <string>
#include <unordered_set>
#include <algorithm>

using namespace Rcpp;

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
    std::string first_id = as<std::string>(lexi_id[0]);
    VE_set.insert(std::to_string(vertex_ids[i]) + ":" + first_id);
    
    std::vector<std::string> faceV, cofaceV;
    for (int j = 1; j < lexi_id.size(); j++) {
      std::string current = as<std::string>(lexi_id[j]);
      (std::count(current.begin(), current.end(), ' ') == 1 ? faceV : cofaceV).push_back(current);
    }
    
    if (faceV.empty()) continue;
    
    std::vector<std::string> PAIR;
    std::unordered_set<std::string> FACEC;
    std::vector<std::string> CRIT;
    
    for (size_t ii = 0; ii < faceV.size(); ii++) {
      const std::string& face = faceV[ii];
      std::vector<std::string> matches;
      for (const auto& coface : cofaceV) {
        if (coface.find(face) != std::string::npos) matches.push_back(coface);
      }
      
      if (ii > 0) {
        for (const auto& p : PAIR) {
          if (p.find(face) != std::string::npos && !matches.empty()) {
            CRIT.push_back(matches[0]);
            goto next_face;
          }
        }
      }
      
      if (!matches.empty()) {
        std::string coface_used = matches[0];
        PAIR.push_back(face + ":" + coface_used);
        FACEC.insert(face);
        cofaceV.erase(std::remove(cofaceV.begin(), cofaceV.end(), coface_used), cofaceV.end());
      } else {
        CRIT.push_back(face);
      }
      next_face:;
    }
    
    for (const auto& c : cofaceV) CRIT.push_back(c);
    for (const auto& f : faceV) if (!FACEC.count(f)) CRIT.push_back(f);
    for (const auto& p : PAIR) VE_set.insert(p);
    for (const auto& c : CRIT) CR_set.insert(c);
  }
  
  std::vector<std::string> VE_final(VE_set.begin(), VE_set.end());
  std::vector<std::string> CR_final(CR_set.begin(), CR_set.end());
  std::sort(CR_final.begin(), CR_final.end());
  
  return List::create(_["VE_"] = VE_final, _["CR_"] = CR_final);
}