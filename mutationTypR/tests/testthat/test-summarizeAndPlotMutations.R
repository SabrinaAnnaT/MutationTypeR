library(testthat)
library(VariantAnnotation)
library(ggplot2)
library(BSgenome.Hsapiens.UCSC.hg19)
library(mutationTypR)

context("Testing summarizeAndPlotMutations Function")

test_that('summarizeAndPlotMutations returns correct data frame structure', {
  vcffile <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
  vcf <- readVcf(vcffile, "hg19")
  vcf <- GenomeInfoDb::renameSeqlevels(vcf, 'chr22')

  summary_df <- summarizeAndPlotMutations(vcf, Hsapiens, context_length = 3, plot = FALSE, top_n = 20)

  # Check if the result is a data frame
  expect_true(is.data.frame(summary_df))

  # Check if the data frame has the correct columns
  expected_cols <- c("MutationWithContext", "Frequency", "MutationType")
  expect_equal(names(summary_df), expected_cols)

  # Check if the top_n filtering is applied correctly
  expect_true(nrow(summary_df) <= 20)
})

test_that('summarizeAndPlotMutations handles incorrect inputs', {
  # Assuming that incorrect inputs should result in an error
  expect_error(summarizeAndPlotMutations("invalid_vcf", Hsapiens, 3),
               "Error message for invalid vcf input",
               info = "Should give an error with invalid vcf input.")

  expect_error(summarizeAndPlotMutations(vcf, "invalid_genome", 3),
               "Error message for invalid genome object",
               info = "Should give an error with invalid genome object.")
})

