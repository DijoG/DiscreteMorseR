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
    readr::write_tsv(mesh_ver, stringr::str_c(txt_dirout, "/trees_mesh_ver.txt"))
  }
  
  # Edges (1-simplex)
  ldf = list(mesh$edgesDF %>%
               data.table::data.table() %>%
               dplyr::select(-dplyr::any_of(c("angle", "exterior", "coplanar"))), mesh_ver)
  
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
    readr::write_tsv(mesh_edge_final, stringr::str_c(txt_dirout, "/trees_mesh_edge.txt"))
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
    readr::write_tsv(mesh_face_final, stringr::str_c(txt_dirout, "/trees_mesh_face.txt"))
  }
  
  return(list(
    vertices = mesh_ver,
    edges = mesh_edge_final,
    faces = mesh_face_final
  ))
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
  
  lowerSTAR = list()
  for (i in 1:nrow(vertex)) {
    vertTO = get_vertTO_cpp(vertex[i,], edge, face) %>%
      data.table::data.table() 
    
    if (nrow(vertTO) > 0) {
      vertTO = vertTO %>% 
        dplyr::arrange(order(gtools::mixedorder(lexi_label, decreasing = FALSE)))
      
      tf = purrr::map(vertTO$lexi_label, ~ stringr::str_split_1(., " ")[1] %>% as.double) <= as.double(vertex[i,]$Z)
      
      lowerSTAR[[i]] = vertTO[tf,] %>% 
        dplyr::mutate(id = vertex$i123[i], .before = 1)
      
      if (!is.null(dirout)) {
        data.table::fwrite(lowerSTAR[[i]], 
                           file = stringr::str_c(dirout, "/lowerSTAR.txt"),
                           sep = "\t", append = TRUE, col.names = FALSE)
      }
    }
  }
  
  if (is.null(dirout)) {
    return(lowerSTAR)
  }
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

#' Compute lower star in parallel (simple and reliable)
#'
#' @param vertex Vertex data
#' @param edge Edge data
#' @param face Face data
#' @param output_dir Output directory
#' @param cores Number of cores (default: available cores-1)
#' @param batch_size Number of vertices per batch
#' @return Combined lower star results
#' @export
compute_lowerSTAR_parallel <- function(vertex, edge, face, output_dir = NULL, 
                                       cores = NULL, batch_size = 1000) {
  
  if (!requireNamespace("furrr", quietly = TRUE)) {
    stop("Package 'furrr' required for parallel processing")
  }
  
  # Auto-detect cores if not specified
  if (is.null(cores)) {
    cores = parallel::detectCores() - 1
  }
  
  n_vertex = nrow(vertex)
  
  # Calculate batches
  if (is.null(batch_size)) {
    batch_size = ceiling(n_vertex / (cores * 3))
    batch_size = max(batch_size, 1000)
    batch_size = min(batch_size, 50000)
  }
  
  batches = split(1:n_vertex, ceiling(seq_along(1:n_vertex) / batch_size))
  
  message("Parallel computation: ", n_vertex, " vertices, ", 
          length(batches), " batches, ", cores, " cores")
  
  # Set up parallel processing
  future::plan(future::multisession, workers = cores)
  
  # Simple batch processor using the original get_lowerSTAR
  process_batch = function(batch_indices) {
    batch_vertex = vertex[batch_indices, ]
    get_lowerSTAR(batch_vertex, edge, face, output_dir, cores = 1)
  }
  
  # Process batches in parallel
  results = furrr::future_map(batches, process_batch, 
                              .options = furrr::furrr_options(seed = TRUE))
  
  # Combine results
  final_results = unlist(results, recursive = FALSE)
  
  message("Parallel computation complete: ", length(final_results), " lower star sets")
  
  return(final_results)
}

#' Compute Morse complex from mesh
#'
#' @param mesh Mesh object from MeshesOperations
#' @param output_dir Optional directory to save results
#' @param parallel Whether to use parallel processing (default: TRUE)
#' @param cores Number of cores for parallel processing (default: 4)
#' @return List with Morse vector field and critical simplices
#' @export
compute_morse_complex <- function(mesh, output_dir = NULL, 
                                  parallel = TRUE, 
                                  cores = 4) {
  
  message("Step 1: Computing simplices")
  simplices = get_SIMPLICES(mesh, output_dir)
  
  message("Step 2: Computing lower star filtration")
  if (parallel) {
    lower_star = compute_lowerSTAR_parallel(
      simplices$vertices, simplices$edges, simplices$faces, 
      output_dir, cores = cores
    )
  } else {
    lower_star = get_lowerSTAR(simplices$vertices, simplices$edges, simplices$faces, output_dir)
  }
  
  message("Step 3: Processing Morse pairings")
  morse_complex = proc_lowerSTAR(lower_star, simplices$vertices)
  
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
  
  if (!requireNamespace("rgl", quietly = TRUE)) {
    stop("Package 'rgl' required for visualization. Install with: install.packages('rgl')")
  }
  
  if (!requireNamespace("MeshesOperations", quietly = TRUE)) {
    stop("Package 'MeshesOperations' required for mesh conversion")
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
