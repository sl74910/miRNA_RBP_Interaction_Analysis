# miRNA_RBP_Interaction_Analysis  

## Overview  
This repository contains scripts used to analyze complex interactions between microRNAs (miRNAs) and RNA-binding proteins (RBPs), with a focus on their shared mRNA targets. The analysis combines a partial information decomposition (PID) method with additional data merging and filtering processes, revealing the cooperative and competitive relationships between miRNAs and RBPs. 

## Files in the Repository  
- **RBP_miRNA_Merger.cpp**: A C++ script that merges and analyzes miRNA-RBP interactions based on mRNA targets. This script performs filtering, data merging, and further analysis of overlapping and adjacent interactions between miRNAs and RBPs.
- **PIDSyR.R**: The script implementing the PID method to analyze miRNA-RBP interactions.
- **plot.R**: Script for visualizing the interaction patterns between miRNAs and RBPs.

## Methodology  
The RBP_miRNA_Merger.cpp script begins by filtering out interactions where the RBP and mRNA genes are identical. Following this, it merges miRNA-mRNA and RBP-mRNA interaction datasets, producing intermediate files that facilitate further analysis. The merged dataset enables the identification of potential competitive or synergistic interactions between miRNAs and RBPs that target the same mRNA.

To dissect these interactions, the Partial Information Decomposition (PID) method is applied. This technique breaks down the total information that miRNAs and RBPs provide about their shared mRNA target into distinct components, offering a nuanced view of how these molecules collaborate or compete in regulating gene expression.The analysis leverages both sequence data and gene expression profiles, integrating multiple cancer-related datasets to investigate the universality and specificity of miRNA-RBP interactions across different cancer types.

## Results  
The analysis reveals that miRNAs and RBPs often exhibit synergistic regulatory effects, working together to inhibit the expression of target mRNAs. The results suggest that these interactions play a significant role in gene expression regulation at the post-transcriptional level, particularly in cancer.



