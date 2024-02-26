#' Extract and Validate Mutation Types from VCF Data
#'
#' This function extracts and validates mutation types from a VCF file, considering a specified
#' context length around each mutation site. It ensures that only single nucleotide variants
#' (SNVs) are considered and reports mutations with reference bases C or T only, avoiding
#' redundancy by converting mutations with reference bases A or G to their reverse complements.
#'
#' @param vcf_input File path to the VCF file or a CollapsedVCF object.
#' @param ref_genome Reference genome object.
#' @param context_length Integer specifying the total length of the genomic context,
#' including the mutated nucleotide.
#' @return A vector of mutation types.
#' @importFrom VariantAnnotation readVcf
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings DNAString
#' @importFrom Biostrings reverseComplement
#' @examples
#' vcf<- system.file("extdata", "chr22.vcf.gz", package = "VariantAnnotation")
#' vcf_data <- readVcf(vcf, genome = "hg19")
#' vcf_data <- renameSeqlevels(vcf_data, 'chr22')
#' Hsapiens <- BSgenome.Hsapiens.UCSC.hg19
#' extractMutationTypes(vcf_data, Hsapiens, 3)
#'
#' @export
extractMutationTypes <- function(vcf_input, ref_genome, context_length) {
  # Validate and read VCF
  if (is(vcf_input, "CollapsedVCF")) {
    vcf_data <- vcf_input
  } else if (is.character(vcf_input) && file.exists(vcf_input)) {
    vcf_data <- readVcf(vcf_input, genome = ref_genome)
  } else {
    stop("vcf_input must be a valid file path or a CollapsedVCF object.")
  }

  # Extract mutation data
  mutation_data <- data.frame(
    chrom = seqnames(vcf_data),
    pos = start(vcf_data),
    ref_base = as.character(ref(vcf_data)),
    alt_base = as.character(unlist(alt(vcf_data))),
    stringsAsFactors = FALSE
  )

  # Filter out invalid nucleotides and non-SNVs, and ensure positions are numeric
  valid_nucleotides <- c("A", "C", "G", "T")
  mutation_data <- mutation_data[!is.na(mutation_data$ref_base) & !is.na(mutation_data$alt_base), ]
  mutation_data <- mutation_data[mutation_data$ref_base %in% valid_nucleotides & mutation_data$alt_base %in% valid_nucleotides, ]
  mutation_data <- mutation_data[nchar(mutation_data$ref_base) == 1 & nchar(mutation_data$alt_base) == 1, ]
  if (!all(is.numeric(mutation_data$pos))) {
    stop("Positions must be numeric.")
  }

  # Adjust context length and extract sequences
  flanking_bases <- (context_length - 1) / 2
  mutation_ranges <- GRanges(
    seqnames = mutation_data$chrom,
    ranges = IRanges(start = mutation_data$pos - flanking_bases, end = mutation_data$pos + flanking_bases)
  )

  # Get sequences surrounding mutations
  surrounding_seq <- getSeq(ref_genome, mutation_ranges)
  mutation_types <- character(length(surrounding_seq))

  # Generate mutation types
  for (idx in seq_along(surrounding_seq)) {
    context_seq <- surrounding_seq[idx]
    ref_base <- mutation_data$ref_base[idx]
    alt_base <- mutation_data$alt_base[idx]

    # Handle reverse complement for mutations with ref base A or G
    if (ref_base %in% c("A", "G")) {
      context_seq <- reverseComplement(context_seq)
      ref_base <- as.character(reverseComplement(DNAString(ref_base)))
      alt_base <- as.character(reverseComplement(DNAString(alt_base)))
    }

    # Ensure all mutation types have C or T as the mutated reference base and construct the mutation type string
    mutation_types[idx] <- paste0(
      substr(as.character(context_seq), 1, flanking_bases),
      "[", ref_base, ">", alt_base, "]",
      substr(as.character(context_seq), flanking_bases + 2, context_length)
    )
  }

  return(mutation_types)
}

#' Summarize and Plot Mutation Types with Context
#'
#' Summarizes the mutation types extracted from a VCF file and plots
#' the frequency of each mutation type categorized by reference and mutated nucleotides.
#' It displays a refined plot showing only the most common mutation types for legibility,
#' while returning a complete frequency table.
#'
#' @param vcf_file File path to the VCF file or a CollapsedVCF object.
#' @param ref_genome Reference genome object.
#' @param context_length Integer specifying the length of the genomic context (default: 3).
#' @param plot Logical indicating whether to plot the mutation counts (default: TRUE).
#' @param top_n Integer specifying the number of most common mutation types to display in the plot (default: 20).
#' @return A data frame containing the counts of each mutation type with context.
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 theme_minimal
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 coord_flip
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 labs
#' @examples
#' summarizeAndPlotMutations(vcf, Hsapiens,3)
#'
#' @export
summarizeAndPlotMutations <- function(vcf_file, ref_genome, context_length = 3, plot = TRUE, top_n = 20) {
  mutation_types <- extractMutationTypes(vcf_file, ref_genome, context_length)

  # Create a summary table of all mutation types with context
  mutation_summary <- table(mutation_types)
  mutation_summary_df <- as.data.frame(mutation_summary, stringsAsFactors = FALSE)
  names(mutation_summary_df) <- c("MutationWithContext", "Frequency")

  # Create a copy of the data frame to use for plotting
  plot_df <- mutation_summary_df

  # Order by frequency for plotting and select top_n most common types
  plot_df <- plot_df[order(-plot_df$Frequency), ]
  if(nrow(plot_df) > top_n) {
    plot_df <- head(plot_df, top_n)
  }

  # Extract the mutation type for grouping
  plot_df$MutationType <- gsub(".*\\[(.*)\\].*", "\\1", plot_df$MutationWithContext)

  # Define colors f
  colors <- c("C>T" = "mediumblue", "T>C" = "skyblue", "C>A" = "salmon", "T>A" = "orange", "C>G" = "forestgreen", "T>G" = "plum")

  # Prepare the plot
  p <- ggplot(plot_df, aes(x = MutationWithContext, y = Frequency, fill = MutationType)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = colors) +
    coord_flip() +
    theme_minimal() +
    theme(
      axis.text.y = element_text(angle = 0, hjust = 1),
      plot.title = element_text(size = 16, face = "bold"),
      legend.position = "bottom"
    ) +
    labs(
      title = "Top Mutation Context Frequencies",
      x = "Mutation Context",
      y = "Frequency",
      fill = "Mutation Type"
    ) +
    facet_wrap(~ MutationType, scales = "free", ncol = 1) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

  # Display the plot
  if (plot) {
    print(p)
  }

  # Return the full summary table
  return(mutation_summary_df)
}
