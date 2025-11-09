#include <Rcpp.h>
#include <regex>
#include <vector>

using namespace Rcpp;

// [[Rcpp::export]]
List get_vOUT(DataFrame v, CharacterVector ve_e_p) {
  
  CharacterVector vertex_i123 = v["i123"];
  
  int n = vertex_i123.size();
  List vOUT;
  NumericVector region(n);  // 'region' vector to store the region for each iteration

  for (int i = 0; i < n; i++) {
    Rcout << "_____" << i + 1 << "_____\r";

    std::string pattern = "\\b" + as<std::string>(vertex_i123[i]) + "\\b";
    std::regex reg(pattern);

    std::vector<std::string> vpe_;
    for (int j = 0; j < ve_e_p.size(); j++) {
      if (std::regex_search(as<std::string>(ve_e_p[j]), reg)) {
        vpe_.push_back(as<std::string>(ve_e_p[j]));
      }
    }

    std::vector<std::string> vpe_flat;
    for (std::vector<std::string>::size_type j = 0; j < vpe_.size(); j++) {
      std::regex reg("[0-9]+");
      std::sregex_iterator it(vpe_[j].begin(), vpe_[j].end(), reg);
      std::sregex_iterator end;
      
      while (it != end) {
        vpe_flat.push_back(it->str());
        ++it;
      }
    }

    std::vector<std::string> vtofind;
    for (std::vector<std::string>::size_type j = 0; j < vpe_flat.size(); j++) {
      if (!std::regex_search(vpe_flat[j], reg)) {
        vtofind.push_back(vpe_flat[j]);
      }
    }

    std::vector<std::string> vFOUND;
    if (!vtofind.empty()) {
      while (!vtofind.empty()) {
        std::string find = vtofind[0];
        vtofind.erase(vtofind.begin());
        
        std::vector<std::string> vpe__;
        for (int j = 0; j < ve_e_p.size(); j++) {
          std::regex reg("\\b" + find + "\\b");
          
          if (std::regex_search(as<std::string>(ve_e_p[j]), reg)) {
            vpe__.push_back(as<std::string>(ve_e_p[j]));
          }
        }
        std::vector<std::string> vpe__flat;
        for (std::vector<std::string>::size_type j = 0; j < vpe__.size(); j++) {
          std::regex reg("[0-9]+");
          std::sregex_iterator it(vpe__[j].begin(), vpe__[j].end(), reg);
          std::sregex_iterator end;
          
          while (it != end) {
            vpe__flat.push_back(it->str());
            ++it;
          }
        }
        for (std::vector<std::string>::size_type j = 0; j < vpe__flat.size(); j++) {
          if (!std::regex_search(vpe__flat[j], reg)) {
            vFOUND.push_back(vpe__flat[j]);
          }
        }
      }
    } else {
      // Continue to the next iteration if vtofind is empty
      continue;
    }

    // Adding minima
    std::sort(vFOUND.begin(), vFOUND.end());
    auto last = std::unique(vFOUND.begin(), vFOUND.end());
    vFOUND.erase(last, vFOUND.end());
    
    std::vector<double> vFOUND_double;
    std::transform(vFOUND.begin(), vFOUND.end(), std::back_inserter(vFOUND_double), [](const std::string& str) {
      return std::stod(str);
    });

    // Dequeuing visited and found vertices
    std::vector<std::string> ve_e_p_filtered;
    for (int j = 0; j < ve_e_p.size(); j++) {
      if (!std::regex_search(as<std::string>(ve_e_p[j]), reg)) {
        ve_e_p_filtered.push_back(as<std::string>(ve_e_p[j]));
      }
    }

    // Storing
    DataFrame df = DataFrame::create(Named("vFOUND") = vFOUND_double, Named("region") = i + 1);
    vOUT.push_back(df);
  }

  return vOUT;
}
