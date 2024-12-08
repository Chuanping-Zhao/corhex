# corhex

corhex is an R package for generating correlation hexagon plots. It helps you quickly visualize correlations between multiple groups within the R environment, particularly suited for large-scale datasets.

# Installation

```
# Install devtools
install.packages("devtools")

# Install corHex package from GitHub
devtools::install_github("Chuanping-Zhao/corhex")
```
# Example Data

dtplot:

![image](https://github.com/user-attachments/assets/57e4b70d-a923-4aef-ac8a-1302e3d50f7e)



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

![image](https://github.com/user-attachments/assets/39a24ccb-7a14-411d-a84a-6e9c3d162e78)

