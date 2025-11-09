#include <regex> 
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame get_vertTO_cpp(DataFrame vertex, DataFrame edge, DataFrame face) {
  int n = vertex.nrow();

  CharacterVector vertex_i123 = vertex["i123"];
  CharacterVector edge_lexi_id = edge["lexi_id"];
  CharacterVector edge_lexi_label = edge["lexi_label"];
  CharacterVector face_lexi_id = face["lexi_id"];
  CharacterVector face_lexi_label = face["lexi_label"];

  CharacterVector vertTO_lexi_label;
  CharacterVector vertTO_lexi_id;

  for (int i = 0; i < n; i++) {
    // Use std::regex to match the pattern
    std::string pattern = "\\b" + as<std::string>(vertex_i123[i]) + "\\b";
    std::regex regexPattern(pattern);

    for (int j = 0; j < edge.nrow(); j++) {
      if (std::regex_search(as<std::string>(edge_lexi_id[j]), regexPattern)) {
        vertTO_lexi_id.push_back(edge_lexi_id[j]);
        vertTO_lexi_label.push_back(edge_lexi_label[j]);
      }
    }

    for (int j = 0; j < face.nrow(); j++) {
      if (std::regex_search(as<std::string>(face_lexi_id[j]), regexPattern)) {
        vertTO_lexi_id.push_back(face_lexi_id[j]);
        vertTO_lexi_label.push_back(face_lexi_label[j]);
      }
    }
  }

  DataFrame vertTO = DataFrame::create(Named("lexi_label") = vertTO_lexi_label, 
                                       Named("lexi_id") = vertTO_lexi_id);

  return vertTO;
}

