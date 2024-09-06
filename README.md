# nestedIDR
**Description**: The package fits the idr model for hierarchical structured data (e.g ChIP-seq data from different labs with multiple replicates) to measure the reproducibility and replicability of findings from replicate experiments in multi-source studies. Using a nested copula mixture model that characterizes the interdependence between replication experiments both across and within sources, the hierarchical idr quantifies reproducibility and replicability of each candidate simultaneously in a coherent framework.

This project reports R code for implementing the methods and reproducing the results in Ranalli et al. "A statistical framework for measuring reproducibility and replicability of high-throughput experiments from multiple sources".

The following libraries shall be installed in R before using code from
this project:

library(mixtools)

library(idr)
prova 
