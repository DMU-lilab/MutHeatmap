# MutHeatmap
Mutation heatmap with multiple mutation types

Plot a heat map with multiple mutation types within a gene

The gene heat map that we usually see is a grid of genes, a type of mutation,
but in fact the same gene often has multiple mutation types in the same patient,
so the traditional heat map drawing tool can't satisfy us. For the problem, this
function uses image to draw the preliminary heat map, and then uses points to
add the second mutation in the form of a square. The third mutation is added in
turn. At the same time, the square position is slightly moved and accompanied by
a slight reduction in size to achieve a better display effect, and up to four
mutations can be represented on one heat map grid.

### Before run MutHeatmap, one should known:

1. If you want to draw a heat map with a histogram, because the image size is too large, use the pdf
function and give it a large enough width and length;

2. The default is the type of mutation annotated by annovar;

3. Because there are too many factors affecting the alignment of the heat map and the histogram during
drawing, it is difficult to adjust the corresponding mar, mex, and oma parameters to achieve better
results. Therefore, it is recommended to quickly draw a rough estimate and then use inkscape or
adobe for layout alignment.

4. If the point of the mutation type in the heat map is too small, the width and length of the pdf file
should be reduced.

### installation
#### Option 1(Recommended)
Using devtools
```
install.packages("devtools")
library(devtools)
install_github("DMU-lilab/MutHeatmap")
```
### Option 2
Manual installation
#### Checkout the latest release of MutHeatmap from GitHub
```git clone https://github.com/DMU-lilab/MutHeatmap.git```
#### Install R dependencies (in R)
 ```install.packages("data.table") # version > 1.10.4```

#### Install the FindPeek package
From the command line and in the directory where FindPeek github was cloned.
```R CMD INSTALL MutHeatmap ```

### Usage
```?MutHeatmap```to access the documentation pages for function MutHeatmap
