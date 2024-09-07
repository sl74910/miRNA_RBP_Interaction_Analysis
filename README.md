# miRNA_RBP_Interaction_Analysis  

## Overview  
This repository contains scripts used to analyze complex interactions between microRNAs (miRNAs) and RNA-binding proteins (RBPs), with a focus on their shared mRNA targets. The analysis combines a partial information decomposition (PID) method with additional data merging and filtering processes, revealing the cooperative and competitive relationships between miRNAs and RBPs. The project provides insights into gene expression regulation, particularly within the context of cancer biology.

## Files in the Repository  
- **RBP_miRNA_Merger.cpp**: A C++ script that merges and analyzes miRNA-RBP interactions based on mRNA targets. This script performs filtering, data merging, and further analysis of overlapping and adjacent interactions between miRNAs and RBPs.
- **PIDSyR.R**: The script implementing the PID method to analyze miRNA-RBP interactions.
- **plot.R**: Script for visualizing the interaction patterns between miRNAs and RBPs.

## Methodology  
The RBP_miRNA_Merger.cpp script first filters out interactions where the RBP and mRNA genes are the same. It then merges miRNA-mRNA and RBP-mRNA interaction datasets, generating intermediate files. This merged data enables the detection of potential competitive or synergistic interactions between miRNAs and RBPs targeting the same mRNA.The PID method is used to decompose the total information that miRNAs and RBPs provide about their shared target mRNA into distinct components. This allows for a detailed understanding of how these molecules interact to regulate gene expression. The analysis is based on sequence data and gene expression profiles, integrating multiple cancer datasets to explore the universality and specificity of miRNA-RBP interactions.

## Results  
The analysis reveals that miRNAs and RBPs often exhibit synergistic regulatory effects, working together to inhibit the expression of target mRNAs. The results suggest that these interactions play a significant role in gene expression regulation at the post-transcriptional level, particularly in cancer.



