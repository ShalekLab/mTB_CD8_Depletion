## Code for scRNA-seq analysis included in Winchell et al 2023


Preprocessing, clustering and QC code for all analysis is in:
0_Preprocessing.ipynb


Code to generate analyses for each figure:

Figure 5: General clustering and macrophage and cDC subclustering. For iterative clustering code see preprocessing notebook.
A-D: Figure5_all_cell_clustering.ipynb
E-H: Figure5_macrophage_cell_clustering.ipynb
J-M: Figure5_cDC_clustering.ipynb
N-I: Functional_enrichment.ipynb


Figure 6: Figure6_tnk_clustering.ipynb

Figure 7:
A. Differential expression within Cytotoxic NK/T cells between anti-CD8a depleted samples and control: pseudobulk_differential_expression.ipynb

B. Violin plots of lineage genes between groups
Figure6_tnk_clustering.ipynb

C. Gene ontology analysis: Tcell_functional_enrichment.ipynb

Figure 8: 
Interaction analysis details: 
1. Generate pseudobulk counts for clusters of interest
2. Generate csv containing percent of cells expressing each gene for each cluster
3. Run DESeq2 to get differential expressed genes between conditions then run NicheNet on those lists - Nichenet_pseudobulk_interaction_analysis.ipynb
4. Calculate interaction potentials and annotate NicheNET list based on these significant potentials - Figure8_Receptor_Ligand_Interactions.ipynb


A,B: Figure8_Receptor_Ligand_interactions.ipynb
C. Figure8_S10.ipynb 
 



Supplemental:
S4: See preprocessing notebook for details
A. Preprocessing notebook
B. Figure5_all_cell_and_myeloid_clustering.ipynb
C. 0_Preprocessing.ipynb
D. 0_Preprocessing.ipynb
E. 0_Preprocessing.ipynb
F.Figure6_tnk_clusteing.ipynb

S5: See figure 8 description for interaction analysis details (A-E).
F-G: Figure8_S10.ipynb



Pseudobulk differential expression: pseudobulk_differential_expression.ipynb


Differential abundance: DA_with_edgeR.ipynb
