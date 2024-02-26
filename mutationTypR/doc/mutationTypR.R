## ----style, echo = FALSE, results = 'asis'------------------------------------
library(BiocStyle)

## ----echo = FALSE-------------------------------------------------------------
library(knitr)

## ----echo = FALSE-------------------------------------------------------------
library(mutationTypR)

## -----------------------------------------------------------------------------
library(VariantAnnotation)
library(GenomeInfoDb)
library(BSgenome.Hsapiens.UCSC.hg19)

## -----------------------------------------------------------------------------
vcffile <- system.file("extdata", "chr22.vcf.gz", package = "VariantAnnotation")
vcf_data <- VariantAnnotation::readVcf(vcffile, genome = "hg19")
vcf_data <- GenomeInfoDb::renameSeqlevels(vcf_data, 'chr22')
Hsapiens <- BSgenome.Hsapiens.UCSC.hg19

## -----------------------------------------------------------------------------
mutation_types <- extractMutationTypes(vcf_data, Hsapiens, 3)

## -----------------------------------------------------------------------------
summary_df <- summarizeAndPlotMutations(vcf_data, Hsapiens, context_length = 3, plot = TRUE, top_n = 20)
print(summary_df)

