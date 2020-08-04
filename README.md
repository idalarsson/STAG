# State Transitions and Growth (STAG) model
Genome-scale modeling of temporal dynamics in glioblastoma

## Requirements
Matlab (R2019a)  
CVX: Matlab Software for Disciplined Convex Programming

## Input files
STAG requires an input matrix where each row represents a single cell, and the columns are genes (opt.), barcode, timepoint, treatment and state (opt.) according to  

Genes: Normalized expression values, if states should be assigned within the pipeline.
Barcode: The barcode retrieved for each individual cell  
Timepoint: The timepoint the cell was sampled  
Treatment: Single digit (1,2,...,n). If no treatment, 1 for all cells.  
State: Single digit (1,2,...,n) if this has been assigned before. Can also be assigned within the pipeline.  

## Pipeline
Follow the script main.m

## Output
A-matrix with transition and growth rates
