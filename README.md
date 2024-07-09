### All code used during my MRes/PhD. Some are to be inputed straight into the terminal in GENOMICS ENGLAND others ar scripts, Jupyter workbooks or R studio codes. 

Binary_enhancer_specificty_heatmap.R : R studio code used to create a heatmap using binary data for ATAC signals in enhancer regions of in house Endothelial cell ATAC and Encode ATAC. 


Cohort_variant_search.sh : This is a script used in the GEL terminal to automatically find desired variants within the cohort.


DE_NOVO_SEARCH.txt : This is the code used in the GEL terminal for the DE novo variant filtering process. 

ECs_combidning_bedfiles.sh : This is a script that combines bed files form the inhouse Endothelial cell ATAC seqeuncing. The output was used to generate a binary heatmap of LEC specific enhancer regions.The code for the heatmap is found in Binary_enhancer_specificty_heatmap.R

    
Encode_CRES_SPEC.sh : This is a script that combines ENCODE broadpeak files and adds a column with the corresponding cell name. This is then intersected with the defined enhancer regions. This was used to generate a file with the replicate ATAC peaks combined which was then used for a binary heatmap. The code for the heatmap is found in Binary_enhancer_specificty_heatmap.R


Enhancer_specificity_heatmaps.ipynb : A visual studio code workbook which includes all  code used to generate Enhancer specificity heatmaps based of ATAC signal strength. The submit_dbaCount.sh file is used in this workbook. The code for the heatmap is found in Binary_enhancer_specificty_heatmap.R

GEL_relatedeness.txt : This is the code used in the GEL terminal to find the kinship coefficients and determine relatedness between participants. 

Graphtyper.txt : This is a script used in the GEL terminal to look for variants within regions of a gene in all participants in the cohort. 

Motifbreakr.R : This is an R studio code to find predicted TF motifs that may be altered owing to SNPs. 

Submit_dbaCount.sh : This is a script used within the Enhancer_specificity_heatmaps.ipynb to generate the dbObj.counted file needed for the heatmap.



