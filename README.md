# DiscreteMorseR 🚀

The DiscreteMorseR package delivers ultra-fast C++ backend Morse gradient field and critical simplices (0-simplices: vertices, 1-simplices: edges, 2-simplices: faces) parallel computation. 

# Installation

```r
devtools::install_github("DijoG/DiscreteMorseR")
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
  - `AlphaHull3D` - Alpha hull generation
  - `clustermq` - Parallelization backend 

# Usage

## Data Preparation
```r
library(lidR);library(tidyverse);library(data.table)

# Load LiDAR data
ilas <- lidR::readLAS("D:/Gergo/DiscreteMorseR/lasref/12tree_exampleN.las")
trees <- lidR::filter_poi(ilas, treeid %in% c(2:7))
plot(trees, pal = "grey98")

# Create matrix input for alpha hull
lasdf <- 
  ilas@data[,1:3] %>%
  as.data.frame() %>%
  distinct() %>%
  as.matrix()

# Generate alpha hull
a <- AlphaHull3D::ahull3d(lasdf, alpha = .1) 

# Extract largest connected component mesh
mesh <- DiscreteMorseR::get_MESH(vertices = a$vb, faces = a$it)
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
# ~3.5 minutes for typical TLS tree point clouds
```
<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/DMR/DMRproc.png" width="650">

**🚀 Performance Highlights:**
- ✅ **226,267 vertices** processed in parallel  
- ✅ **12 cores** utilized (~38% CPU efficiency)
- ✅ **100% completion rate** - all lower star sets computed
- ✅ **Complete Morse analysis** in ~3.5 minutes
- ✅ **Automatic file export** of all results

# Analyze results
```r
crit_types <- sapply(strsplit(morse_complex$critical, " "), length)
table(crit_types)
# 1 = vertices, 2 = edges, 3 = faces
```
## Visualization
```r
# Gradient field only
p <- DiscreteMorseR::visualize_MORSE_2d(
  morse_complex, 
  projection = "XZ",
  point_alpha = .6,
  point_size = .8,
  plot_critical = FALSE,
  max_points = 30000
)
print(p)
```
<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/DMR/DMR01.png" width="650">

```r
# Critical simplices only  
pp <- DiscreteMorseR::visualize_MORSE_2d(
  morse_complex, 
  projection = "XZ",
  point_alpha = .6,
  point_size = .8,
  plot_gradient = FALSE,
  max_points = 30000
)
print(pp)
```
<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/DMR/DMR02.png" width="650">

```r
# Multi-panel: all projections
ppp <- DiscreteMorseR::visualize_MORSE_2d_panel(
  morse_complex,
  point_alpha = .5,
  point_size = .8,
  plot_gradient = FALSE,
  max_points = 30000
)
print(ppp)
```
<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/DMR/DMR03.png" width="650">