#motifbreakR for intronic variants

#Install MotifBreakR
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")

BiocManager::install("motifbreakR", force = TRUE)

library("motifbreakR")

# your SNPs from bed file
setwd('/Users/tamaravujic/Desktop/ERG/ERG_SNP_TF_MOTIF_analysis/rs148479581')
snps.bed.file <- "rs148479581.bed"

# Read in Single Nucleotide Variants
library(BSgenome)
available.SNPs()

#install and load SNPs from rsID
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
install.packages("/Users/tamaravujic/Desktop/ERG/ERG_SNP_analysis_TV/MotifBreakR/SNPlocs.Hsapiens.dbSNP155.GRCh38", repos = NULL, type = "source")

library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(BSgenome.Hsapiens.UCSC.hg38) 

#collect SNPs - checking for rsIDs in dbsNP 155
#snps.mb <- snps.from.file(file = snps.bed.file,
#                               dbSNP = SNPlocs.Hsapiens.dbSNP155.GRCh38,
#                               search.genome = BSgenome.Hsapiens.UCSC.hg38,
#                              format = "bed", check.unnamed.for.rsid = TRUE)

#collect SNPs - NOT checking for rsIDs in dbsNP 155
snps.mb <- snps.from.file(file = snps.bed.file,
                          search.genome = BSgenome.Hsapiens.UCSC.hg38,
                          format = "bed", check.unnamed.for.rsid = FALSE)

snps.mb

#find broken motifs
BiocManager::install("MotifDb", force = TRUE)
library(MotifDb)
MotifDb


#We have leveraged the MotifList introduced by MotifDb to include an additional set of motifs that have been gathered to this package:
data(motifbreakR_motif)
motifbreakR_motif

#access all data sets for motifBreakR analysis 
data(hocomoco)
data(homer)
data(encodemotif)
data(factorbook)
homer
hocomoco
encodemotif
factorbook

#run motifbreakR - database hocomoco
results.hocomoco <- motifbreakR(snpList = snps.mb[1:1], hocomoco, filterp = TRUE,
                       threshold = 1e-4,
                       method = "ic",
                       bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                       BPPARAM = BiocParallel::SerialParam())

#run motifbreakR - database homer
results.homer <- motifbreakR(snpList = snps.mb[1:1], homer, filterp = TRUE,
                                threshold = 1e-4,
                                method = "ic",
                                bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                                BPPARAM = BiocParallel::SerialParam())

#run motifbreakR - database encodemotif
results.encodemotif <- motifbreakR(snpList = snps.mb[1:1], encodemotif, filterp = TRUE,
                             threshold = 1e-4,
                             method = "ic",
                             bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                             BPPARAM = BiocParallel::SerialParam())

#run motifbreakR - database factorbook
results.factorbook <- motifbreakR(snpList = snps.mb[1:1], factorbook, filterp = TRUE,
                                   threshold = 1e-4,
                                   method = "ic",
                                   bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                                   BPPARAM = BiocParallel::SerialParam())

################################################################################
#view results for variant 
################################################################################

rs2836373.hocomoco <- results.hocomoco[results.hocomoco$SNP_id == "chr21:37850535:C:T"]
rs2836373.hocomoco

rs2836373.encodemotif <- results.encodemotif[results.encodemotif$SNP_id == "chr21:37850535:C:T"]
rs2836373.encodemotif

rs2836373.homer <- results.homer[results.homer$SNP_id == "chr21:37850535:C:T"]
rs2836373.homer

rs2836373.factorbook <- results.factorbook[results.factorbook$SNP_id == "chr21:37850535:C:T"]
rs2836373.factorbook

#Calculate P value
rs2836373.hocomoco <- calculatePvalue(rs2836373.hocomoco, granularity = 1e-6)
rs2836373.hocomoco


rs2836373.encodemotif <- calculatePvalue(rs2836373.encodemotif, granularity = 1e-6)
rs2836373.encodemotif

rs2836373.homer <- calculatePvalue(rs2836373.homer, granularity = 1e-6)
rs2836373.homer

rs2836373.factorbook <- calculatePvalue(rs2836373.factorbook, granularity = 1e-6)
rs2836373.factorbook

#plot variant rs2836373
install.packages("Cairo")
plotMB(results = results.hocomoco, rsid = "chr21:37850535:C:T", effect = "strong", altAllele = "T")

plotMB(results = results.encodemotif, rsid = "chr21:37850535:C:T", effect = "strong", altAllele = "T")

plotMB(results = results.factorbook, rsid = "chr21:37850535:C:T", effect = "strong", altAllele = "T")

plotMB(results = results.homer, rsid = "chr21:37850535:C:T", effect = "strong", altAllele = "T")



################################################################################
# a txt file with the session info and package versions - wil be needed when time to write 
writeLines(capture.output(sessionInfo()), "_sessionInfo.txt")

