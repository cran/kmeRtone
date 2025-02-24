# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Count k-mers from string(s) using a simple hash table.
#' 
#' Count only observed k-mers. Biostrings::oligonucleotideFrequency reports all
#' possible k-mers. For k > 12, the memory for creating empty k-mer counts
#' spiked and crashed the R session.
#'
#' @param sequences Sequence strings.
#' @param k Size of k-mer.
#' @return A vector of k-mer counts. The counts of multiple sequences are
#'    combined, similar to Biostrings::oligonucleotideFrequency simplify.as
#'    "collapsed".
#'
#' @export
countKmers <- function(sequences, k) {
    .Call(`_kmeRtone_countKmers`, sequences, k)
}

#' Locate a middle sequence pattern and count its sequence context.
#'
#' This function searches for a specified middle pattern within a given sequence. 
#' It then counts the occurrences of specific context patterns within a defined window
#' size around the middle pattern. The function returns a map where keys are the 
#' counts of context patterns found and values are the frequencies of these counts.
#'
#' @param sequence A string representing the sequence to be analyzed.
#' @param mid_pattern A string representing the middle pattern to search for within the sequence.
#' @param window An integer specifying the size of the surrounding window around the middle pattern.
#' @param context_patterns A vector of strings representing the context patterns to search for within the window.
#' @return A std::unordered_map<int,int> where keys are the counts of context patterns found 
#'         and values are the frequencies of these counts.
#'
#'
#' @export
countMidPatternContext2 <- function(sequence, mid_pattern, window, context_patterns) {
    .Call(`_kmeRtone_countMidPatternContext2`, sequence, mid_pattern, window, context_patterns)
}

#' Count Relevant K-mers with Specified Middle Pattern from Sequence String(s)
#'
#' This function scans through each sequence in the provided vector, locating a specified middle pattern.
#' For each occurrence of the middle pattern, the function extracts and counts the surrounding k-mers. 
#' The k-mers are identified based on the given k-mer size and centered around the middle pattern.
#'
#' @param sequences A vector of strings, each representing a sequence to be analyzed.
#' @param k An integer specifying the size of the k-mers to be extracted and counted.
#' @param mid_pattern A string representing the middle pattern to search for within each sequence.
#' @return A std::unordered_map with k-mers as keys and their counts as values.
#'
#' @export
countMidPatternKmers <- function(sequences, k, mid_pattern) {
    .Call(`_kmeRtone_countMidPatternKmers`, sequences, k, mid_pattern)
}

#' Ccount sequence context of given point positions.
#'
#' @param sequence A sequence to slide.
#' @param points Middle point positions.
#' @param len Length of the middle point.
#' @param window Size of a surrounding window.
#' @param context_patterns Context patterns to search for.
#' @return A named vector of frequency of counts.
#'
#' @export
countPointContext2 <- function(sequence, points, len, window, context_patterns) {
    .Call(`_kmeRtone_countPointContext2`, sequence, points, len, window, context_patterns)
}

#' Count k-mers in given ranges of a sequence.
#'
#' Slide and update the cummulated table count.
#'
#' @param sequence A sequence to count.
#' @param starts Start positions.
#' @param ends End positions.
#' @param k K-mer size.
#' @return A k-mer-named vector of count.
#'
#' @export
countRangedKmers <- function(sequence, starts, ends, k) {
    .Call(`_kmeRtone_countRangedKmers`, sequence, starts, ends, k)
}

#' Count sequence content in a sliding window of a sequence.
#'
#' @param sequence A sequence to slide.
#' @param window Size of a window.
#' @param pattern A pattern to search for.
#' @return A numeric vector of count.
#'
#' @export
countSlidingWindow <- function(sequence, window, pattern) {
    .Call(`_kmeRtone_countSlidingWindow`, sequence, window, pattern)
}

#' Count sequence content in a sliding window of a sequence.
#'
#' @param sequence A sequence to slide.
#' @param window Size of a window.
#' @param patterns Patterns of the same size to search for.
#' @return Named vector of frequency of count.
#'
#' @export
countSlidingWindow2 <- function(sequence, window, patterns) {
    .Call(`_kmeRtone_countSlidingWindow2`, sequence, window, patterns)
}

#' Count sequence content in a given sequence.
#'
#' stringi has no function that search within substring without memory copy it.
#' This function has two versions. One without the need to memory copy denoted
#' as ***. The only downside to this is std::string::find cannot stop searching
#' past end of substring. I manage to at least stop it as soon as possible. If
#' the pattern is long and rare, it won't stop until it find post-substring
#' pattern. The other version is memory copy substring but as this operation is
#' in the loop, the memory is still within comfortable range. c++17 has
#' std::string_view that solve this but still new and not widely available. Use
#' count_substring_regex to avoid memory copy.
#'
#' @param sequence A sequence to map.
#' @param start Start positions.
#' @param end End positions.
#' @param pattern A pattern to search for.
#' @return A numeric vector of count.
#'
#' @export
count_substring_fixed <- function(sequence, start, end, pattern) {
    .Call(`_kmeRtone_count_substring_fixed`, sequence, start, end, pattern)
}

#' Count sequence content in a given sequence.
#'
#' stringi has no function that search within substring without memory creating
#' it. This function solve that. Unlike count_substring_fixed, this function
#' does not need to memory copy substring.
#'
#' @param sequence A sequence to map.
#' @param start Start positions.
#' @param end End positions.
#' @param pattern A regex pattern to search for within start and end positions.
#' @return A numeric vector of count.
#'
#' @export
count_substring_regex <- function(sequence, start, end, pattern) {
    .Call(`_kmeRtone_count_substring_regex`, sequence, start, end, pattern)
}

#' Simulate a population given ranges of chromosome sequence to mutate.
#'
#' @param chrom_seq A chromosome sequence.
#' @param starts Start positions.
#' @param ends End positions.
#' @param strand Strand type: "+" or "-".
#' @param snv_df A table of SNV frequency. Columns: position, base, count.
#' @param pop_size Size of population.
#' @param top_kmers Extreme k-mers i.e. highly susceptible k-mers.
#' @param central_pattern K-mer central pattern.
#' @param k K-mer size.
#' @return A count matrix with 4 rows for total top k-mers and susceptible
#'    k-mers in sense and antisense. Columns correspond to population
#'    individuals.
#'
simulatePopulation <- function(chrom_seq, starts, ends, strand, snv_df, pop_size, top_kmers, central_pattern, k) {
    .Call(`_kmeRtone_simulatePopulation`, chrom_seq, starts, ends, strand, snv_df, pop_size, top_kmers, central_pattern, k)
}

