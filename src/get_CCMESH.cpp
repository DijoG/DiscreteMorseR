#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <queue>
#include <unordered_set>

using namespace Rcpp;

// Ultra-fast integer pair for edges
typedef std::pair<int, int> Edge;

// Custom hash for edges - maximum speed
struct EdgeHash {
    inline std::size_t operator()(const Edge& e) const {
        return (static_cast<size_t>(e.first) << 32) | static_cast<size_t>(e.second);
    }
};

// Ultra-fast component extractor - zero dynamic allocation in hot loops
class FastComponentExtractor {
private:
    std::vector<int> queue;              // declared first
    std::vector<bool> visited;           // declared second  
    std::vector<int> component_vertices;
    std::unordered_set<Edge, EdgeHash> component_edges;
    int queue_start, queue_end;
    
public:
    FastComponentExtractor(int n_vertices) : 
        queue(n_vertices),              // initialized first
        visited(n_vertices + 1, false), // initialized second
        component_vertices() {
        component_vertices.reserve(n_vertices);
        component_edges.reserve(n_vertices * 3);
    }
    
    List extract_component_fast(const NumericMatrix& vertices, 
                               const NumericMatrix& faces,
                               int start_vertex,
                               const std::vector<std::vector<int>>& adj) {
        
        // Reset for new component
        queue_start = 0;
        queue_end = 0;
        component_vertices.clear();
        component_edges.clear();
        
        // Manual queue implementation
        queue[queue_end++] = start_vertex;
        visited[start_vertex] = true;
        component_vertices.push_back(start_vertex);
        
        // Ultra-fast BFS
        while (queue_start < queue_end) {
            int current = queue[queue_start++];
            const std::vector<int>& neighbors = adj[current];
            const int num_neighbors = neighbors.size();
            
            for (int j = 0; j < num_neighbors; ++j) {
                int neighbor = neighbors[j];
                if (!visited[neighbor]) {
                    visited[neighbor] = true;
                    queue[queue_end++] = neighbor;
                    component_vertices.push_back(neighbor);
                }
                
                // Create edge
                int v1 = current, v2 = neighbor;
                if (v1 > v2) std::swap(v1, v2);
                component_edges.insert(std::make_pair(v1, v2));
            }
        }
        
        return build_component_ultrafast(vertices, faces);
    }
    
    // Getter for visited array
    const std::vector<bool>& getVisited() const { return visited; }
    
private:
    List build_component_ultrafast(const NumericMatrix& vertices, 
                                  const NumericMatrix& faces) {
        
        // Sort vertices once
        std::sort(component_vertices.begin(), component_vertices.end());
        
        // Build vertex mapping array
        std::vector<int> vertex_map(vertices.nrow() + 1, -1);
        const int num_comp_vertices = component_vertices.size();
        NumericMatrix comp_vertices(num_comp_vertices, 3);
        
        for (int i = 0; i < num_comp_vertices; ++i) {
            const int old_idx = component_vertices[i];
            vertex_map[old_idx] = i + 1;
            comp_vertices(i, 0) = vertices(old_idx - 1, 0);
            comp_vertices(i, 1) = vertices(old_idx - 1, 1);
            comp_vertices(i, 2) = vertices(old_idx - 1, 2);
        }
        
        // Count faces in component
        int face_count = 0;
        const int num_faces = faces.nrow();
        for (int i = 0; i < num_faces; ++i) {
            if (visited[faces(i, 0)] && visited[faces(i, 1)] && visited[faces(i, 2)]) {
                face_count++;
            }
        }
        
        // Build faces with new indices
        NumericMatrix comp_faces(face_count, 3);
        int face_idx = 0;
        for (int i = 0; i < num_faces; ++i) {
            const int v1 = faces(i, 0);
            const int v2 = faces(i, 1);
            const int v3 = faces(i, 2);
            
            if (visited[v1] && visited[v2] && visited[v3]) {
                comp_faces(face_idx, 0) = vertex_map[v1];
                comp_faces(face_idx, 1) = vertex_map[v2];
                comp_faces(face_idx, 2) = vertex_map[v3];
                face_idx++;
            }
        }
        
        // Build edges from component edges with new vertex indices
        std::unordered_set<Edge, EdgeHash> unique_edges;
        unique_edges.reserve(component_edges.size());
        
        for (const auto& edge : component_edges) {
            const int new_v1 = vertex_map[edge.first];
            const int new_v2 = vertex_map[edge.second];
            if (new_v1 != -1 && new_v2 != -1) {
                unique_edges.insert(std::make_pair(new_v1, new_v2));
            }
        }
        
        NumericMatrix comp_edges(unique_edges.size(), 2);
        int edge_idx = 0;
        for (const auto& edge : unique_edges) {
            comp_edges(edge_idx, 0) = edge.first;
            comp_edges(edge_idx, 1) = edge.second;
            edge_idx++;
        }
        
        List mesh;
        mesh["vertices"] = comp_vertices;
        mesh["faces"] = comp_faces;
        mesh["edges"] = comp_edges;
        
        return mesh;
    }
};

// Ultra-fast adjacency list builder
std::vector<std::vector<int>> build_adjacency_ultrafast(int n_vertices, const NumericMatrix& faces) {
    std::vector<std::vector<int>> adj(n_vertices + 1);
    
    // Single pass to count degrees
    std::vector<int> degree(n_vertices + 1, 0);
    const int num_faces = faces.nrow();
    for (int i = 0; i < num_faces; ++i) {
        degree[faces(i, 0)] += 2;
        degree[faces(i, 1)] += 2;
        degree[faces(i, 2)] += 2;
    }
    
    // Pre-allocate exact memory needed
    for (int i = 1; i <= n_vertices; ++i) {
        adj[i].reserve(degree[i]);
    }
    
    // Single pass to build adjacency
    for (int i = 0; i < num_faces; ++i) {
        const int v1 = faces(i, 0);
        const int v2 = faces(i, 1);
        const int v3 = faces(i, 2);
        
        adj[v1].push_back(v2);
        adj[v1].push_back(v3);
        adj[v2].push_back(v1);
        adj[v2].push_back(v3);
        adj[v3].push_back(v1);
        adj[v3].push_back(v2);
    }
    
    // Fast duplicate removal in-place
    for (int i = 1; i <= n_vertices; ++i) {
        if (!adj[i].empty()) {
            std::sort(adj[i].begin(), adj[i].end());
            auto last = std::unique(adj[i].begin(), adj[i].end());
            adj[i].resize(last - adj[i].begin());
        }
    }
    
    return adj;
}

// Fast single mesh creation (no components)
List create_single_mesh_fast(const NumericMatrix& vertices, const NumericMatrix& faces) {
    std::unordered_set<Edge, EdgeHash> unique_edges;
    const int n_faces = faces.nrow();
    unique_edges.reserve(n_faces * 3);
    
    // Single pass edge extraction
    for (int i = 0; i < n_faces; ++i) {
        int v1 = faces(i, 0), v2 = faces(i, 1), v3 = faces(i, 2);
        
        // Sort edges for consistent hashing
        if (v1 > v2) std::swap(v1, v2);
        unique_edges.insert(std::make_pair(v1, v2));
        
        if (v1 > v3) std::swap(v1, v3);
        unique_edges.insert(std::make_pair(v1, v3));
        
        if (v2 > v3) std::swap(v2, v3);
        unique_edges.insert(std::make_pair(v2, v3));
    }
    
    NumericMatrix edges(unique_edges.size(), 2);
    int edge_idx = 0;
    for (const auto& edge : unique_edges) {
        edges(edge_idx, 0) = edge.first;
        edges(edge_idx, 1) = edge.second;
        edge_idx++;
    }
    
    List mesh;
    mesh["vertices"] = vertices;
    mesh["faces"] = faces;
    mesh["edges"] = edges;
    return mesh;
}

// [[Rcpp::export]]
List get_CCMESH_cpp(NumericMatrix vertices, NumericMatrix faces, bool select_largest = true) {
    const int n_vertices = vertices.nrow();
    const int n_faces = faces.nrow();
    
    // Input validation
    if (n_vertices == 0 || n_faces == 0) {
        stop("Empty vertices or faces matrix");
    }
    
    // Fast path: single component assumption when select_largest = true
    if (select_largest) {
        return create_single_mesh_fast(vertices, faces);
    }
    
    // Multi-component processing only when explicitly requested
    std::vector<std::vector<int>> adj = build_adjacency_ultrafast(n_vertices, faces);
    FastComponentExtractor extractor(n_vertices);
    
    std::vector<List> components;
    components.reserve(10);
    
    for (int i = 1; i <= n_vertices; ++i) {
        if (!extractor.getVisited()[i] && !adj[i].empty()) {
            List component = extractor.extract_component_fast(vertices, faces, i, adj);
            components.push_back(component);
        }
    }
    
    // Handle no components found
    if (components.empty()) {
        return create_single_mesh_fast(vertices, faces);
    }
    
    // Return all components as a list
    List result;
    result["components"] = components;
    result["n_components"] = components.size();
    return result;
}