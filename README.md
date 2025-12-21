# DiscreteMorseR 🚀

[![Parallel](https://img.shields.io/badge/Parallel-20+%20cores-green.svg)]()
[![C++](https://img.shields.io/badge/C++-Optimized-blue.svg)]()
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

The DiscreteMorseR package delivers ultra-fast C++ backend Morse gradient field and critical simplices (0-simplices: vertices, 1-simplices: edges, 2-simplices: faces) parallel computation. Perfect for LiDAR data, computational topology, and Morse theory applications.

# Installation

```r
devtools::install_github("DijoG/DiscreteMorseR")
devtools::install_github("DijoG/ahull3D")
library(DiscreteMorseR)
```
Quick test of C++ functions
```r
DiscreteMorseR::get_MIXEDSORT_cpp(c("2", "1", "12 45", "25 256", "11 8", "256 23"))
DiscreteMorseR::add_DECIMAL(215.2585589, 3)
```
# Dependencies

  - `lidR` - LiDAR data processing
  - `tidyverse` - Data manipulation
  - `data.table` - Fast data operations
  - `clustermq` - Parallelization backend 

# Usage

## Data Preparation
```r
library(lidR);library(tidyverse);library(data.table)

# Load LiDAR data
trees <- lidR::readLAS("D:/Gergo/DiscreteMorseR/lasref/12tree_exampleN.las")
lidR::plot(trees, pal = "grey98")

# Create matrix input for alpha hull
lasdf <- 
  ilas@data[, c("X", "Y", "Z", "pid")] %>%
  as.data.frame() %>%
  distinct(X, Y, Z, .keep_all = TRUE) %>%
  as.matrix()

# Generate alpha hull
a <- ahull3D::ahull3D(lasdf[,1:3], 
  input_truth = lasdff[,4], 
  alpha = .1) 

# Extract largest connected component mesh
mesh <- DiscreteMorseR::get_CCMESH(a)
```
## Morse Complex Analysis

**Real-world computation on tree point cloud (226,267 vertices):**
```r
tictoc::tic()
morse_complex <- DiscreteMorseR::compute_MORSE_complex(
  mesh, 
  output_dir = "D:/Gergo/DiscreteMorseR/12_output",
  cores = 12,
  batch_size = 5000  # Increase for large datasets
)  
tictoc::toc()
# ~2.5 minutes for typical TLS tree point clouds
```
<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/DMR/DMRprocv.png" width="800">

**🚀 Performance Highlights:**
- ✅ **226,267 vertices** processed in parallel  
- ✅ **12 cores** utilized (~92% CPU efficiency)
- ✅ **100% completion rate** - all lower star sets computed
- ✅ **Complete Morse analysis** in ~2.5 minutes
- ✅ **Automatic file export** of all results

## Analyze results 
```r
crit_types <- sapply(strsplit(morse_complex$critical, " "), length)
table(crit_types)
# 1 = vertices (0-simplices, minima), 2 = edges (1-simplices), 3 = faces (2-simplices)
crit_types
     1      2      3 
225137    115   1005   

```
## Visualization
```r
# Critical simplices only  
p <- DiscreteMorseR::visualize_MORSE_2d(
  morse_complex, 
  projection = "XZ",
  point_alpha = .6,
  point_size = .8,
  plot_gradient = FALSE,
  max_points = 30000
)
print(p)
```
<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/DMR/DMR_xz.png" width="800">

```r
# Multi-panel: all projections
pp <- DiscreteMorseR::visualize_MORSE_2d(
  morse_complex, 
  projection = "XY",
  point_alpha = .6,
  point_size = .8,
  plot_gradient = FALSE,
  max_points = 30000
)
print(pp)
```
<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/DMR/DMR_xy.png" width="800" height="800">

```r
# Multi-panel: all projections
ppp <- DiscreteMorseR::visualize_MORSE_2d_panel(
  morse_complex, 
  point_alpha = .6,
  point_size = .8,
  plot_gradient = FALSE,
  max_points = 30000
)
print(ppp)
```
<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/DMR/DMR_3.png" width="800" height="800">

## Save Visualization
```r
DiscreteMorseR::save_MORSE_2d(
  morse_complex,
  filename = "D:/Gergo/DiscreteMorseR/png/DMR_xz.png",
  projection = "XZ",
  point_alpha = .6,
  point_size = .8,
  plot_gradient = F,
  max_points = 30000,
  width = 6,
  height = 5
)
```