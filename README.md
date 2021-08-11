# State Transitions and Growth (STAG) model <img src="https://user-images.githubusercontent.com/43134255/128992022-60ca394f-741a-45a3-8759-818d8b457cd9.png" width="100" height="100">


Modeling glioblastoma heterogeneity as a dynamic network of cell states

## Requirements
Matlab (R2019a)  
CVX: Matlab Software for Disciplined Convex Programming

## Input files
STAG requires an input matrix where each row represents a single cell, and the columns are genes (opt.), barcode, timepoint, treatment and state (opt.) according to:  

Genes: Normalized expression values, if states should be assigned within the pipeline.  
Barcode: The barcode retrieved for each individual cell.  
Timepoint: The timepoint the cell was sampled.  
Treatment: Single digit (1,2,...,n). If no treatment, 1 for all cells.  
State: Single digit (1,2,...,n) if this has been assigned before. Can also be assigned within the pipeline.  

## Provided data
The file "expressionMatrices_20210622.mat" includes barcoded data for the three GBM cell lines U3017MG, U3065MG and U3071MG. 

## Pipeline
Follow the script main.m

## Output
A-matrix with transition and growth rates
