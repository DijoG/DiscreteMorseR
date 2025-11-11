#include <Rcpp.h>
#include <vector>
#include <string>
#include <unordered_map>
#include <cstring>

using namespace Rcpp;

// Fast string view approach
class StringView {
private:
  const char* start_;
  size_t length_;
public:
  StringView(const char* start, size_t length) : start_(start), length_(length) {}
  
  std::string str() const { return std::string(start_, length_); }
  
  // For unordered_map compatibility
  bool operator==(const StringView& other) const {
    return length_ == other.length_ && 
           std::strncmp(start_, other.start_, length_) == 0;
  }
};

// Hash for StringView
struct StringViewHash {
  std::size_t operator()(const StringView& sv) const {
    // Simple hash - for better performance consider FNV-1a
    std::size_t h = 0;
    for (size_t i = 0; i < sv.str().length(); ++i) {
      h = h * 31 + sv.str()[i];
    }
    return h;
  }
};

// [[Rcpp::export]]
List get_PRECOMPUTEDvert_cpp(CharacterVector lexi_ids, CharacterVector lexi_labels) {
  int n = lexi_ids.size();
  
  // Use pre-allocation for better performance
  std::unordered_map<std::string, std::vector<int>> vertex_map;
  vertex_map.reserve(n * 3); // Reserve space to avoid rehashing
  
  NumericVector first_verts(n);
  
  for (int i = 0; i < n; i++) {
    // MAX SPEED: First vertex extraction
    const char* label = lexi_labels[i];
    first_verts[i] = std::atof(label); // atof automatically stops at space
    
    // MAX SPEED: Manual tokenization without string copies where possible
    const char* simplex = lexi_ids[i];
    const char* token_start = simplex;
    const char* p = simplex;
    
    while (true) {
      if (*p == ' ' || *p == '\0') {
        if (token_start != p) { // Non-empty token
          // Create string only when necessary
          std::string vertex_id(token_start, p - token_start);
          vertex_map[vertex_id].push_back(i + 1);
        }
        if (*p == '\0') break;
        token_start = p + 1;
      }
      p++;
    }
  }
  
  // Efficient conversion to R data structures
  int map_size = vertex_map.size();
  List vertex_index(map_size);
  CharacterVector names(map_size);
  int idx = 0;
  
  for (const auto& pair : vertex_map) {
    names[idx] = pair.first;
    vertex_index[idx] = IntegerVector(pair.second.begin(), pair.second.end());
    idx++;
  }
  
  vertex_index.attr("names") = names;
  
  return List::create(
    _["vertex_index"] = vertex_index,
    _["first_verts_z"] = first_verts
  );
}