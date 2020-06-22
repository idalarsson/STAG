# Barcoding
Genome-scale modeling of phenotypic switching in glioblastoma

## Requirements
Matlab (R2019a)  
CVX: Matlab Software for Disciplined Convex Programming

## Input files
The scripts require an input matrix where each row represents a single cell, and the columns are genes, barcode, timepoint, treatment and (opt.) state according to  

Genes: Normalized expression values  
Barcode: The barcode retrieved for each individual cell  
Timepoint: The timepoint the cell was sampled  
Treatment: Single digit (1,2,...,n). If no treatment, 1 for all cells.  
State: Single digit (1,2,...,n) if this has been assigned before. Can also be assigned within the code.  

## Pipeline
Follow the script main.m

## Output files

