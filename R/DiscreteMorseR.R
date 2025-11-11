#' Discrete Morse Theory Analysis
#' 
#' @description Compute Morse vector field and critical simplices from 3D meshes.
#' @name DiscreteMorseR
#' @docType package
#' @useDynLib DiscreteMorseR, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @import dplyr purrr stringr data.table gtools readr future furrr
NULL

#' Add decimal formatting
#' 
#' @param x Numeric vector
#' @param k Number of decimal places
#' @return Formatted character vector
#' @keywords internal
add_DECIMAL <- function(x, k) format(round(x, k), nsmall = k)

#' Compute lexicographically sorted simplex labels
#'
#' @param df Data frame with label and idlabel columns
#' @return Data table with lexicographically sorted labels
#' @keywords internal
get_lexIDLAB <- function(df) {
  # Trim labels
  LA = stringr::str_trim(df$label, "both")
  IDLA = stringr::str_trim(df$idlabel, "both")
  
  # Split and get sorting indices
  LA_split = stringr::str_split(LA, " ")
  m1 = purrr::map(LA_split, ~get_MIXEDSORT_cpp(add_DECIMAL(as.numeric(.), k = 9)))
  
  # Apply sorting to both label types
  mlab = purrr::map2(LA_split, m1, ~paste(.x[.y], collapse = " "))
  mid = purrr::map2(stringr::str_split(IDLA, " "), m1, ~paste(.x[.y], collapse = " "))
  
  # Create result data.table
  result = data.table::data.table(
    lexi_label = unlist(mlab),
    lexi_id = unlist(mid)
  )
  
  return(result)
}

#' Extract and process simplices from mesh
#'
#' @param mesh Mesh object from MeshesOperations
#' @param txt_dirout Directory for output files (optional)
#' @return List of processed simplices
#' @keywords internal
get_SIMPLICES <- function(mesh, txt_dirout = "") {
  
  # Vertices (0-simplex)
  mesh_ver = mesh$vertices %>%
    data.table::data.table() %>%
    data.table::setnames(., c("X", "Y", "Z")) %>%
    dplyr::mutate(i123 = 1:dplyr::n(), .before = 1) %>%
    dplyr::distinct() %>% 
    dplyr::mutate(Z = as.character(Z))
  
  if (txt_dirout != "") {
    readr::write_tsv(mesh_ver, stringr::str_c(txt_dirout, "/vertices.txt"))
  }
  
  # Edges (1-simplex)
  ldf = list(mesh$edgesDF %>%
               data.table::data.table(), mesh_ver)
  
  mesh_edge = ldf %>%
    purrr::reduce(dplyr::left_join, by = c("i1" = "i123")) %>%
    dplyr::rename(ii1X = X, ii1Y = Y, ii1Z = Z)
  
  ldf = list(mesh_edge, mesh_ver)
  
  mesh_edge = ldf %>%
    purrr::reduce(dplyr::left_join, by = c("i2" = "i123")) %>%
    dplyr::rename(ii2X = X, ii2Y = Y, ii2Z = Z) %>%
    dplyr::mutate(
      label = stringr::str_c(as.character(ii1Z), as.character(ii2Z), sep = " "),
      idlabel = stringr::str_c(as.character(i1), as.character(i2), sep = " ")
    )
  
  out = get_lexIDLAB(mesh_edge)
  
  # Use cbind!
  mesh_edge_final = cbind(mesh_edge, out)
  
  if (txt_dirout != "") {
    readr::write_tsv(mesh_edge_final, stringr::str_c(txt_dirout, "/edges.txt"))
  }
  
  # Faces (2-simplex)
  ldfa = list(mesh$faces %>%
                data.table::data.table() %>%
                data.table::setnames(., c("i1", "i2", "i3")), mesh_ver)
  
  mesh_f = ldfa %>%
    purrr::reduce(dplyr::left_join, by = c("i1" = "i123")) %>%
    dplyr::rename(ii1X = X, ii1Y = Y, ii1Z = Z)
  
  ldfa = list(mesh_f, mesh_ver)
  
  mesh_f = ldfa %>%
    purrr::reduce(dplyr::left_join, by = c("i2" = "i123")) %>%
    dplyr::rename(ii2X = X, ii2Y = Y, ii2Z = Z)
  
  ldfa = list(mesh_f, mesh_ver)
  
  mesh_face = ldfa %>%
    purrr::reduce(dplyr::left_join, by = c("i3" = "i123")) %>%
    dplyr::rename(ii3X = X, ii3Y = Y, ii3Z = Z) %>%
    dplyr::mutate(
      label = stringr::str_c(as.character(ii1Z), as.character(ii2Z), as.character(ii3Z), sep = " "),
      idlabel = stringr::str_c(as.character(i1), as.character(i2), as.character(i3), sep = " ")
    )
  
  out = get_lexIDLAB(mesh_face)
  
  # Use cbind!
  mesh_face_final = cbind(mesh_face, out)
  
  if (txt_dirout != "") {
    readr::write_tsv(mesh_face_final, stringr::str_c(txt_dirout, "/faces.txt"))
  }
  
  return(list(
    vertices = mesh_ver,
    edges = mesh_edge_final,
    faces = mesh_face_final
  ))
}

#' Ultra-fast mesh preparation 
#' @param vertices Vertex data from alphahull 
#' @param faces Face data from alphahull
#' @return Mesh data expected by compute_morse_complex()
#' @export
get_MESH <- function(vertices, faces) {
  
  vertices = as.matrix(vertices)
  faces = as.matrix(faces)
  
  # Handle alpha hull format
  if(nrow(vertices) == 4) vertices = t(vertices)
  if(nrow(faces) == 3) faces = t(faces)
  if(ncol(vertices) == 4) vertices = vertices[, 1:3]
  if(min(faces) == 0) faces = faces + 1
  
  result = get_MESH_cpp(vertices, faces)
  
  mesh = list(
    vertices = as.matrix(result$vertices),
    faces = as.matrix(result$faces),
    edgesDF = as.data.frame(result$edgesDF)
  )
  
  colnames(mesh$vertices) = c("X", "Y", "Z")
  colnames(mesh$faces) = c("i1", "i2", "i3")
  
  message("Minimal mesh: ", nrow(mesh$vertices), " vertices, ",
          nrow(mesh$faces), " faces, ", nrow(mesh$edgesDF), " edges")
  
  return(mesh)
}

#' Compute lower star filtration for vertices
#'
#' @param vertex Vertex data
#' @param edge Edge data  
#' @param face Face data
#' @param dirout Output directory (optional)
#' @param cores Number of cores (for consistency)
#' @return List of lower star sets
#' @keywords internal
get_lowerSTAR <- function(vertex, edge, face, dirout = NULL, cores = 1) {
  
  # Pre-compute connections
  all_connections = get_vertTO_cpp(vertex, edge, face)
  lexi_labels = all_connections$lexi_label
  lexi_ids = all_connections$lexi_id
  
  # Pre-compute everything upfront
  first_verts = sapply(strsplit(lexi_labels, " "), function(x) as.numeric(x[1]))
  simplex_vertices = strsplit(lexi_ids, " ")
  
  # Build reverse index: vertex_id -> simplex_indices
  vertex_to_simplices = list()
  for (i in seq_along(simplex_vertices)) {
    for (v_id in simplex_vertices[[i]]) {
      if (is.null(vertex_to_simplices[[v_id]])) {
        vertex_to_simplices[[v_id]] = integer(0)
      }
      vertex_to_simplices[[v_id]] = c(vertex_to_simplices[[v_id]], i)
    }
  }
  
  lowerSTAR = vector("list", nrow(vertex))
  
  for (i in 1:nrow(vertex)) {
    v_id = as.character(vertex$i123[i])
    v_z = as.numeric(vertex$Z[i])
    
    simplex_indices = vertex_to_simplices[[v_id]]
    
    if (!is.null(simplex_indices) && length(simplex_indices) > 0) {
      valid_indices = simplex_indices[first_verts[simplex_indices] <= v_z]
      
      if (length(valid_indices) > 0) {
        valid_simplices = all_connections[valid_indices, ]
        
        order_idx = gtools::mixedorder(valid_simplices$lexi_label, decreasing = FALSE)
        sorted_simplices = valid_simplices[order_idx, ]
        
        # Create data.frame WITHOUT factors
        result = data.frame(
          id = rep(as.integer(v_id), nrow(sorted_simplices)),
          lexi_label = sorted_simplices$lexi_label,
          lexi_id = sorted_simplices$lexi_id,
          stringsAsFactors = FALSE
        )
        
        lowerSTAR[[i]] = result
        
        if (!is.null(dirout)) {
          output_lines = apply(result, 1, function(row) {
            paste(row, collapse = "\t")
          })
          
          # Append to file with cat() for exact format matching
          cat(output_lines, file = file.path(dirout, "lowerSTAR.txt"), 
              sep = "\n", append = TRUE)
        }
      }
    }
  }
  
  lowerSTAR = lowerSTAR[!sapply(lowerSTAR, is.null)]
  return(lowerSTAR)
}

#' Process lower star to get Morse pairings
#'
#' @param list_lowerSTAR Lower star data
#' @param vertex Vertex data
#' @return List with vector field and critical simplices
#' @keywords internal
proc_lowerSTAR <- function(list_lowerSTAR, vertex) {
  
  VE_ = c()
  CR_ = c()
  
  for (i in seq_along(list_lowerSTAR)) {
    if (nrow(list_lowerSTAR[[i]]) == 0) {
      CR_ = c(CR_, as.character(vertex[i, "i123"]))
    } else {
      PAIR_v_e_id = stringr::str_c(vertex[i, "i123"], ":", list_lowerSTAR[[i]]$lexi_id[1])
      
      faceV = list_lowerSTAR[[i]]$lexi_id[-1][purrr::map_dbl(purrr::map(list_lowerSTAR[[i]]$lexi_id[-1], ~stringr::str_split_1(., " ")), length) == 2]
      cofaceV = list_lowerSTAR[[i]]$lexi_id[-1][purrr::map_dbl(purrr::map(list_lowerSTAR[[i]]$lexi_id[-1], ~stringr::str_split_1(., " ")), length) == 3]
      
      if (length(faceV) == 0) {
        VE_ = c(VE_, PAIR_v_e_id)
        next
      }
      
      PAIR = c()
      FACEC = c()
      CRIT = c()
      for (ii in seq_along(faceV)) {
        cofaces = cofaceV[purrr::map_lgl(cofaceV, ~stringr::str_detect(., stringr::str_c("^\\b", faceV[ii], "\\b")))]
        
        if (ii > 1 & any(stringr::str_detect(PAIR, faceV[ii]))) {
          CRIT = c(CRIT, cofaces[1])
        } else {
          if (length(cofaces) > 0) {
            PAIR_c = cofaces[1]
            PAIR_  = stringr::str_c(faceV[ii], ":", PAIR_c)
            
            cofaceV = cofaceV[-which(cofaceV == PAIR_c)]
            PAIR    = c(PAIR, PAIR_)
            FACEC   = c(FACEC, faceV[ii])
          } else {
            CRIT = c(CRIT, faceV[ii])
            next
          }
          CRIT = c(cofaceV, setdiff(faceV, FACEC)) %>% 
            gtools::mixedsort()
        }
      }
      
      VE_ = c(VE_, PAIR_v_e_id, PAIR)
      CR_ = c(CR_, CRIT)
    }
  }
  
  return(list(VE_ = VE_ %>% unique(),
              CR_ = CR_ %>% na.omit() %>% unique() %>% gtools::mixedsort()))
}

#' Compute lower star in parallel 
#'
#' @param vertex Vertex data
#' @param edge Edge data
#' @param face Face data
#' @param output_dir Output directory
#' @param cores Number of cores (default: available cores-1)
#' @param batch_size Number of vertices per batch (default minimum: 5000, allowed maximum: 20000)
#' @return Combined lower star results
#' @export
compute_lowerSTAR_parallel <- function(vertex, edge, face, output_dir = NULL, 
                                       cores = NULL, batch_size = 5000) {
  
  if (!requireNamespace("furrr", quietly = TRUE)) {
    stop("Package 'furrr' required for parallel processing")
  }
  
  if (is.null(cores)) {
    cores = parallel::detectCores() - 1
  }
  
  n_vertex = nrow(vertex)
  
  # Batch sizing
  if (is.null(batch_size)) {
    batch_size = ceiling(n_vertex / cores)
    batch_size = max(batch_size, 5000)
    batch_size = min(batch_size, 20000)
  }
  
  batches = split(1:n_vertex, ceiling(seq_along(1:n_vertex) / batch_size))
  total_batches = length(batches)
  
  message("🚀 Rocket-fast parallel: ", n_vertex, " vertices, ", 
          total_batches, " batches, ", cores, " cores")
  
  # PRE-COMPUTE connections to avoid passing large data
  message("Pre-computing vertex connections...")
  all_connections = get_vertTO_cpp(vertex, edge, face)
  
  # Set up parallel processing with increased memory
  options(future.globals.maxSize = 1024 * 1024 * 1024)  # 1GB
  future::plan(future::multisession, workers = cores)
  
  # Optimized batch processor - only passes vertex indices, not large data
  process_batch = function(batch_indices, batch_num) {
    message("Batch ", batch_num, "/", total_batches, " (vertices ", 
            min(batch_indices), "-", max(batch_indices), ")")
    
    batch_results = list()
    for (i in batch_indices) {
      v_id = vertex$i123[i]
      v_z = as.numeric(vertex$Z[i])
      
      # Filter from pre-computed connections
      vertTO = all_connections[grepl(paste0("\\b", v_id, "\\b"), all_connections$lexi_id), ]
      
      if (nrow(vertTO) > 0) {
        vertTO = vertTO[order(gtools::mixedorder(vertTO$lexi_label, decreasing = FALSE)), ]
        
        first_verts = sapply(strsplit(vertTO$lexi_label, " "), function(x) as.numeric(x[1]))
        valid_simplices = first_verts <= v_z
        
        if (any(valid_simplices)) {
          result = cbind(data.frame(id = v_id), vertTO[valid_simplices, ])
          batch_results[[i]] = result
          
          if (!is.null(output_dir)) {
            write.table(result, 
                        file = file.path(output_dir, "lowerSTAR.txt"),
                        sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
          }
        }
      }
    }
    
    batch_results = batch_results[!sapply(batch_results, is.null)]
    message("✅ Batch ", batch_num, "/", total_batches, " complete - ", 
            length(batch_results), " lower star sets")
    
    return(batch_results)
  }
  
  # Add batch numbers
  batch_list <- lapply(seq_along(batches), function(i) {
    list(indices = batches[[i]], number = i)
  })
  
  # Process batches in parallel
  results <- furrr::future_map(batch_list, function(batch) {
    process_batch(batch$indices, batch$number)
  }, .options = furrr::furrr_options(seed = TRUE))
  
  # Combine results
  final_results <- unlist(results, recursive = FALSE)
  
  message("🎉 Parallel computation complete: ", length(final_results), " total lower star sets")
  
  return(final_results)
}

#' Compute Morse complex from mesh
#'
#' @param mesh Mesh object from MeshesOperations
#' @param output_dir Optional directory to save results. If provided, writes:
#'   - vertices.txt, edges.txt, faces.txt: Mesh simplices
#'   - vector_field.txt: Gradient vector field pairs  
#'   - critical_simplices.txt: Critical simplices
#'   - lowerSTAR.txt: Lower star filtration data
#' @param parallel Whether to use parallel processing (default: TRUE)
#' @param cores Number of cores for parallel processing (default: 4)
#' @param batch_size Number of vertices per batch in parallel processing
#' @return List with Morse vector field and critical simplices
#' @export
compute_morse_complex <- function(mesh, output_dir = NULL, parallel = TRUE, 
                                  cores = 4, batch_size = NULL) {
  
  message("Step 1: Computing simplices")
  simplices = get_SIMPLICES(mesh, output_dir)
  
  message("Step 2: Computing lower star filtration")
  if (parallel) {
    lower_star = compute_lowerSTAR_parallel(
      simplices$vertices, simplices$edges, simplices$faces, 
      output_dir, cores = cores, batch_size = batch_size
    )
  } else {
    lower_star = get_lowerSTAR(simplices$vertices, simplices$edges, simplices$faces, output_dir)
  }
  
  message("Step 3: Processing Morse pairings")
  morse_complex = proc_lowerSTAR(lower_star, simplices$vertices)
  
  # Write final results to files if output_dir provided
  if (!is.null(output_dir)) {
    message("Step 4: Writing final results to files")
    
    # Write vector field (one pair per line)
    writeLines(morse_complex$VE_, file.path(output_dir, "vector_field.txt"))
    
    # Write critical simplices (one simplex per line)
    writeLines(morse_complex$CR_, file.path(output_dir, "critical_simplices.txt"))
    
    # Write simplices as tab-separated text files (no row limits)
    write.table(simplices$vertices, file.path(output_dir, "vertices.txt"), 
                sep = "\t", row.names = FALSE, quote = FALSE)
    write.table(simplices$edges, file.path(output_dir, "edges.txt"), 
                sep = "\t", row.names = FALSE, quote = FALSE)
    write.table(simplices$faces, file.path(output_dir, "faces.txt"), 
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    message("Final results written to:")
    message("  - ", file.path(output_dir, "vector_field.txt"))
    message("  - ", file.path(output_dir, "critical_simplices.txt"))
    message("  - ", file.path(output_dir, "vertices.txt"))
    message("  - ", file.path(output_dir, "edges.txt")) 
    message("  - ", file.path(output_dir, "faces.txt"))
    message("  - ", file.path(output_dir, "lowerSTAR.txt"))
  }
  
  message("Discrete Morse Theory analysis complete!")
  return(list(
    vector_field = morse_complex$VE_,    # Gradient vector field
    critical = morse_complex$CR_,        # Critical simplices
    simplices = simplices
  ))
}

#' Get center coordinates of a simplex
#' 
#' @param simplex Simplex identifier string
#' @param vertices Vertex coordinates
#' @return Numeric vector of coordinates
#' @keywords internal
get_simplex_center <- function(simplex, vertices) {
  vert_ids = as.numeric(strsplit(simplex, " ")[[1]])
  vert_coords = vertices[vertices$i123 %in% vert_ids, c("X", "Y", "Z")]
  
  if (nrow(vert_coords) > 0) {
    c(mean(as.numeric(vert_coords$X)),
      mean(as.numeric(vert_coords$Y)),
      mean(as.numeric(vert_coords$Z)))
  } else {
    c(NA, NA, NA)
  }
}

#' Plot gradient vector field
#' 
#' @param vector_field Morse vector field
#' @param vertices Vertex coordinates
#' @param color Arrow color
#' @param length_scale Scale factor for arrow length
#' @param alpha Arrow transparency
#' @keywords internal
plot_gradient_field <- function(vector_field, vertices, color, length_scale, alpha) {
  
  for (pair in vector_field) {
    parts = strsplit(pair, ":")[[1]]
    if (length(parts) == 2) {
      from_simplex = parts[1]
      to_simplex = parts[2]
      
      from_coords = get_simplex_center(from_simplex, vertices)
      to_coords = get_simplex_center(to_simplex, vertices)
      
      if (!any(is.na(from_coords)) && !any(is.na(to_coords))) {
        arrow_vec = (to_coords - from_coords) * length_scale
        scaled_to = from_coords + arrow_vec
        
        rgl::arrows3d(from_coords, scaled_to, 
                      color = color, alpha = alpha, lwd = 2)
      }
    }
  }
}

#' Plot critical points
#' 
#' @param critical Critical simplices
#' @param vertices Vertex coordinates
#' @param colors Colors for different simplex types
#' @param size Point size
#' @keywords internal
plot_critical_points <- function(critical, vertices, colors, size) {
  
  for (crit in critical) {
    parts = strsplit(crit, " ")[[1]]
    
    if (length(parts) == 1) {
      # Critical vertex
      vert_id = as.numeric(parts[1])
      vert_data = vertices[vertices$i123 == vert_id, ]
      if (nrow(vert_data) > 0) {
        rgl::points3d(as.numeric(vert_data$X), 
                      as.numeric(vert_data$Y), 
                      as.numeric(vert_data$Z), 
                      color = colors[1], size = size)
      }
      
    } else if (length(parts) == 2) {
      # Critical edge (plot midpoint)
      edge_verts = as.numeric(parts)
      v1 = vertices[vertices$i123 == edge_verts[1], ]
      v2 = vertices[vertices$i123 == edge_verts[2], ]
      if (nrow(v1) > 0 && nrow(v2) > 0) {
        midpoint = c(
          mean(c(as.numeric(v1$X), as.numeric(v2$X))),
          mean(c(as.numeric(v1$Y), as.numeric(v2$Y))),
          mean(c(as.numeric(v1$Z), as.numeric(v2$Z)))
        )
        rgl::points3d(midpoint[1], midpoint[2], midpoint[3], 
                      color = colors[2], size = size)
      }
      
    } else if (length(parts) == 3) {
      # Critical face (plot centroid)
      face_verts = as.numeric(parts)
      v1 = vertices[vertices$i123 == face_verts[1], ]
      v2 = vertices[vertices$i123 == face_verts[2], ]
      v3 = vertices[vertices$i123 == face_verts[3], ]
      if (nrow(v1) > 0 && nrow(v2) > 0 && nrow(v3) > 0) {
        centroid = c(
          mean(c(as.numeric(v1$X), as.numeric(v2$X), as.numeric(v3$X))),
          mean(c(as.numeric(v1$Y), as.numeric(v2$Y), as.numeric(v3$Y))),
          mean(c(as.numeric(v1$Z), as.numeric(v2$Z), as.numeric(v3$Z)))
        )
        rgl::points3d(centroid[1], centroid[2], centroid[3], 
                      color = colors[3], size = size)
      }
    }
  }
}

#' Visualize Morse complex gradient and critical points
#'
#' @param mesh Mesh object from MeshesOperations
#' @param morse_complex Output from compute_morse_complex()
#' @param vertex_coords Optional: Specific vertex coordinates to highlight
#' @param gradient_color Color for gradient arrows (default: "pink")
#' @param critical_colors Colors for critical points: c(vertex, edge, face)
#' @param arrow_length Scale factor for arrow length (default: 1)
#' @param point_size Size of critical points (default: 8)
#' @param alpha Transparency for mesh and arrows (default: 0.3)
#' @return rgl scene
#' @export
visualize_morse_gradient <- function(mesh, morse_complex, vertex_coords = NULL,
                                     gradient_color = "pink", 
                                     critical_colors = c("firebrick2", "darkslategray1", "forestgreen"),
                                     arrow_length = 1, point_size = 8, alpha = 0.3) {
  
  required_packages = c("rgl", "MeshesOperations")
  missing_packages = required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  
  if (length(missing_packages) > 0) {
    stop("Required packages missing: ", paste(missing_packages, collapse = ", "), 
         ". Install with: install.packages(c('", paste(missing_packages, collapse = "', '"), "'))")
  }
  
  vertices = morse_complex$simplices$vertices
  vector_field = morse_complex$vector_field
  critical = morse_complex$critical
  
  rgl::open3d()
  rgl::bg3d("white")
  
  rgl::shade3d(MeshesOperations::toRGL(mesh), color = "grey80", alpha = alpha)
  
  if (length(vector_field) > 0) {
    plot_gradient_field(vector_field, vertices, gradient_color, arrow_length, alpha)
  }
  
  if (length(critical) > 0) {
    plot_critical_points(critical, vertices, critical_colors, point_size)
  }
  
  if (!is.null(vertex_coords)) {
    rgl::points3d(vertex_coords$X, vertex_coords$Y, vertex_coords$Z, 
                  color = "darkgoldenrod1", size = point_size * 1.5)
  }
  
  rgl::legend3d("topright", 
                legend = c("Gradient", "Critical Vertex", "Critical Edge", "Critical Face"),
                col = c(gradient_color, critical_colors), 
                pch = c(NA, 16, 16, 16),
                lty = c(1, NA, NA, NA),
                cex = 1.2)
  
  message("Morse complex visualization complete!")
  message("Critical points: ", length(critical))
  message("Gradient arrows: ", length(vector_field))
}

#' Create summary visualization of Morse complex
#'
#' @param mesh Mesh object
#' @param morse_complex Output from compute_morse_complex()
#' @param output_file Optional file to save screenshot
#' @return Multi-panel rgl visualization
#' @export
visualize_morse_summary <- function(mesh, morse_complex, output_file = NULL) {
  
  if (!requireNamespace("rgl", quietly = TRUE)) {
    stop("Package 'rgl' required for visualization")
  }
  
  if (!requireNamespace("MeshesOperations", quietly = TRUE)) {
    stop("Package 'MeshesOperations' required for mesh conversion")
  }
  
  rgl::open3d()
  rgl::par3d(mfrow = c(2, 2))
  
  rgl::next3d()
  rgl::shade3d(MeshesOperations::toRGL(mesh), color = "grey80", alpha = 0.3)
  rgl::title3d("Mesh")
  
  rgl::next3d()
  rgl::shade3d(MeshesOperations::toRGL(mesh), color = "grey80", alpha = 0.1)
  plot_critical_points(morse_complex$critical, morse_complex$simplices$vertices, 
                       c("firebrick2", "darkslategray1", "forestgreen"), 8)
  rgl::title3d("Critical Points")
  
  rgl::next3d()
  rgl::shade3d(MeshesOperations::toRGL(mesh), color = "grey80", alpha = 0.1)
  plot_gradient_field(morse_complex$vector_field, morse_complex$simplices$vertices, 
                      "pink", 1, 0.7)
  rgl::title3d("Gradient Field")
  
  rgl::next3d()
  rgl::shade3d(MeshesOperations::toRGL(mesh), color = "grey80", alpha = 0.1)
  plot_gradient_field(morse_complex$vector_field, morse_complex$simplices$vertices, 
                      "pink", 1, 0.5)
  plot_critical_points(morse_complex$critical, morse_complex$simplices$vertices, 
                       c("firebrick2", "darkslategray1", "forestgreen"), 6)
  rgl::title3d("Combined View")
  
  rgl::par3d(mfrow = c(1, 1))
  
  if (!is.null(output_file)) {
    rgl::rgl.snapshot(output_file)
    message("Screenshot saved to: ", output_file)
  }
}
