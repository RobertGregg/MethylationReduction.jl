# MethylationReduction

[![Build Status](https://github.com/RobertGregg/MethylationReduction.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/RobertGregg/MethylationReduction.jl/actions/workflows/CI.yml?query=branch%3Amaster)


The goal of this repository is to investigate dimension reduction algorithms applied to DNA methylation data. These data are large (850,000+ variables) and tend to overwhelm many conventional dimension reduction algorithms (for example, some methods require a correlation matrix to be formed which is infeasible at this scale). Possible methods we plan to investigate include:

- PCA
- Non-negative matrix factorization
- Random projections
- Autoencoders

## PCA

Julia provides a [MultivariateStats.jl](https://github.com/JuliaStats/MultivariateStats.jl) which we can use to do conventional PCA, kernal PCA, etc.

## NMF

Currently, this package implements a GPU compatible NMF algorithm based on Hieracial Alternating Least Squares (HALS). The goal is to decompose the methylation data into two matrices $W$ and $H$:

$$X \approx WH$$

Columns in `W` can be thought of as latent variables which are linear combinations of methylation sites. `H` holds the coefficients that map methylation sites to the latent variables.

Here is an example using the package:

```julia
using MethylationReduction

X = rand(100,100)
sort!(X,dims=1) #to give some structure to X
k = 10
nmf = NMFCache(X, k; α=0.2);
solveNMF(nmf, verbose=false, maxiter=200)
```

The `NMFCache` structure holds two matrices `W` and `H` as well as an estimation of the reconstructed matrix `WHᵀ`. The variable `k` is the number of latent factors we are reducing to. The α value controls the sparsity of `H`, with larger values leading to sparser solutions. The objective function used is: 

$$ \left \lVert X - WH \right \rVert_F + \alpha \left \lVert H \right \rVert_1$$

We can view these matrices using the dot notation similar to python. Additionally, we can observe the error in the reconstruction:

```julia
nmf.WHᵀ
residual(nmf)
```

To run this algorithm on the GPU, you simply need to provide a GPU array, such as those provided by ![CUDA.jl](https://cuda.juliagpu.org/stable/).
```julia
using CUDA

X = CUDA.rand(100,100) #creates a CuMatrix
sort!(X,dims=1) #to give some structure to X
k = 10
nmf = NMFCache(X, k; α=0.2);
solveNMF(nmf, verbose=false, maxiter=200)
```

# Connecting to R

There are a few R libraries that can interface between R and Julia. Here we will demonstrate the `JuliaConnectoR` library.

```R
library(JuliaConnectoR)

#Import the nmf functionality
MethylationReduction <- juliaImport("MethylationReduction")

#Create a random data matrix in R
m <- matrix(runif(1000*100) , ncol = 100)

#Run the NMF algorithm
nmf <- MethylationReduction$NMFCache(m, 20, α=0.1)
MethylationReduction$solveNMF(nmf, verbose=F)

#Extract the latent variables
W <- nmf$W

```
