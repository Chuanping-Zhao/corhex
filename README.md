# corhex

corhex is an R package for generating correlation hexagon plots. It helps you quickly visualize correlations between multiple groups within the R environment, particularly suited for large-scale datasets.

# Installation

```
if(!require(devtools)){install.packages("devtools")}
# Install corHex package from GitHub
devtools::install_github("Chuanping-Zhao/corhex")
```
# Example Data

dtplot:

![image](https://github.com/user-attachments/assets/883b168d-2fa6-4a70-8a72-b15601be5c09)




# Using the cor_hex Function

```
corhex::cor_hex(
dt=dtplot,
id.col="protein",
cor.method=c("pearson", "kendall", "spearman")[1],
savefile="outputfile",
singleplotsize=c(3,2.5),#width height
facetplotsize=c(3*3,2.5*3),#width height
                bin=50
)
```

# result

![image](https://github.com/user-attachments/assets/9a75bb03-9bed-454a-a193-e2e0cf107946)



