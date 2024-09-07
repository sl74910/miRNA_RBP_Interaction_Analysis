# miRNA_RBP_Interaction_Analysis  

## Overview  
This repository contains the R scripts used for analyzing the complex interactions between microRNAs (miRNAs) and RNA-binding proteins (RBPs) in the context of shared mRNA targets. The analysis employs a partial information decomposition (PID) approach to quantify the cooperative and competitive relationships between miRNAs and RBPs, providing insights into their regulatory roles.

## Files in the Repository  
- **PIDSyR.R**: The script implementing the PID method to analyze miRNA-RBP interactions.
- **plot.R**: Script for visualizing the interaction patterns between miRNAs and RBPs.

## Methodology  
The PID method is used to decompose the total information that miRNAs and RBPs provide about their shared target mRNA into distinct components. This allows for a detailed understanding of how these molecules interact to regulate gene expression. The analysis is based on sequence data and gene expression profiles, integrating multiple cancer datasets to explore the universality and specificity of miRNA-RBP interactions.

## Results  
The analysis reveals that miRNAs and RBPs often exhibit synergistic regulatory effects, working together to inhibit the expression of target mRNAs. The results suggest that these interactions play a significant role in gene expression regulation at the post-transcriptional level, particularly in cancer.



