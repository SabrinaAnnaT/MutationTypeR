library(testthat)
library(VariantAnnotation)
library(BSgenome.Hsapiens.UCSC.hg19)

context("Testing extractMutationTypes Function")

# Test case for correct parameters
test_that('Running extractMutationTypes with correct parameters', {
  vcffile <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
  vcf <- readVcf(vcffile, "hg19")
  vcf <- GenomeInfoDb::renameSeqlevels(vcf, 'chr22')
  extracted_mutations <- extractMutationTypes(vcf, Hsapiens, 3)
  expect_length(extracted_mutations, length(vcf))
  for(i in seq_along(extracted_mutations)) {
    expect_match(extracted_mutations[i], '[.>.]')
  }
})

#Test case for incorrect VCF object type
test_that('Running extractMutationTypes with incorrect VCF object type', {
  expect_error(extractMutationTypes("incorrect_object", Hsapiens, 3))
})

#Test case for incorrect reference genome object
test_that('Running extractMutationTypes with incorrect reference genome object', {
  vcffile <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
  vcf <- readVcf(vcffile, "hg19")
  vcf <- GenomeInfoDb::renameSeqlevels(vcf, 'chr22')

  expect_error(extractMutationTypes(vcf, "incorrect_genome", 3))
})

#Test case for invalid context length (non-integer/string)
test_that('Running extractMutationTypes with invalid context length', {
  vcffile <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
  vcf <- readVcf(vcffile, "hg19")
  vcf <- GenomeInfoDb::renameSeqlevels(vcf, 'chr22')

  expect_error(extractMutationTypes(vcf, Hsapiens, "invalid_context_length"))
})

#Test case for even context length
test_that('Running extractMutationTypes with even context length', {
  vcffile <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
  vcf <- readVcf(vcffile, "hg19")
  vcf <- GenomeInfoDb::renameSeqlevels(vcf, 'chr22')
  expect_error(extractMutationTypes(vcf, Hsapiens, 4))
})

#Test case for checking the format of mutation types
test_that('Checking the format of mutation types', {
  vcffile <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
  vcf <- readVcf(vcffile, "hg19")
  vcf <- GenomeInfoDb::renameSeqlevels(vcf, 'chr22')

  mutation_types <- extractMutationTypes(vcf, Hsapiens, 3)

  for(i in seq_along(mutation_types)) {
    expect_match(mutation_types[i], pattern = "[[A,C,G,T]>[A,C,G,T]]")
  }
})

