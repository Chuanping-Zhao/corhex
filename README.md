# corhex

corhex is an R package for generating correlation hexagon plots. It helps you quickly visualize correlations between multiple groups within the R environment, particularly suited for large-scale datasets.

# Installation

```
if(!require(devtools)){install.packages("devtools")}
if(!require(spatstat.geom)){install.packages("spatstat.geom")}#dependency package
# Install corHex package from GitHub
devtools::install_github("Chuanping-Zhao/corhex")
```
# Example Data

```
library(corhex)
data("demo", package = "corhex")
```
```
output>
demo
# A tibble: 847 × 7
   Index      replicate1 replicate2 replicate3 replicate4 replicate5 replicate6
   <chr>           <dbl>      <dbl>      <dbl>      <dbl>      <dbl>      <dbl>
 1 protein_1       38860      66779      83308      51912      45539      53776
 2 protein_2      302630     345410     265360     240090     325560     316900
 3 protein_3     9404800   10450000   10737000    9861500   12372000   10057000
 4 protein_4      290580     196550     568370     411740     583660     250190
 5 protein_5      148370     124730     134680     103870     117820     110700
 6 protein_6    36017000   34625000   39045000   30246000   33735000   26250000
 7 protein_7     2092400    1997200    2525800    1582500    2855600    1745800
 8 protein_8      251060     236500     431380     130350     179880      49719
 9 protein_9     1559300    1428300     959220    1322900    1302300     780380
10 protein_10     144310      69745      71703      69709      66113      55467
# ℹ 837 more rows
# ℹ Use `print(n = ...)` to see more rows
```





# Using the cor_hex Function

```
result=cor_hex(dt=demo,
               id.col="Index",
               cor.method=c("pearson", "kendall", "spearman")[3],
               savefile="outputfile",
               singleplotsize=c(3,2.5),
               facetplotsize=c(3*3,2.5*3),
               bin=50,#Number of bins for hex plot. Default is 50.
               logL=FALSE,#Whether to log-transform the data (log10). Default is TRUE.
               pointcolor=c("A","B","C","D","E")[3],#Color scheme for density point plots. Options include "A", "B", "C", "D", "E". Default is "C".
               pointsize=0.5,#Size of points in density point plots. Default is 0.5
               kde2d.n=50,#Grid size for kernel density estimation using `MASS::kde2d`. Default is 50.
               plottype=c("hex","point")[2])#Type of plot to generate: "hex" for hexagon plots, "point" for density point plots. Default is "point"
```

result:

![facetwrapPlot_cor](https://github.com/user-attachments/assets/33b3ea49-0427-44c0-8869-3e4a69889e7f)

# Using the corpairs Function

```
corpairs(dt=demo,
         id.col="Index",
         cor.method=c("pearson", "kendall", "spearman")[3],
         savefile="outputfile_corpairs",
         corpairplotsize=c(3*3,2.5*3),#width height
         bin=50,
         logL=FALSE,
         plottype=c("hex","point")[2],
         pointcolor=c("A","B","C","D","E")[3],
         pointsize=1,
         kde2d.n=50)
```

plotting result:

![pairsCorplot](https://github.com/user-attachments/assets/7ba7574b-54e3-42ab-bab3-09deec941108)