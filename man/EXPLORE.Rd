% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EXPLORE.R
\name{EXPLORE}
\alias{EXPLORE}
\title{Function generates various exploratory analyses.}
\usage{
EXPLORE(
  case.coor.path,
  genome.name,
  strand.sensitive,
  k,
  case.pattern,
  output.path,
  case,
  genome,
  control,
  genome.path,
  single.case.len,
  rm.dup,
  case.coor.1st.idx,
  coor.load.limit,
  genome.load.limit,
  genome.fasta.style,
  genome.ncbi.db,
  use.UCSC.chr.name,
  verbose
)
}
\arguments{
\item{case.coor.path}{Path to case coordinates.}

\item{genome.name}{Genome name (e.g., hg19, hg38).}

\item{strand.sensitive}{Boolean indicating if strand sensitivity is considered.}

\item{k}{K-mer size.}

\item{case.pattern}{String patterns to consider in the analysis.}

\item{output.path}{Output directory path for exploration plots.}

\item{case}{Coordinate class object or similar structure for case data.}

\item{genome}{Genome class object or similar structure.}

\item{control}{Control class object or similar structure.}

\item{genome.path}{Path to genome fasta files.}

\item{single.case.len}{Length of single cases.}

\item{rm.dup}{Boolean indicating if duplicates should be removed.}

\item{case.coor.1st.idx}{Indexing of case coordinates.}

\item{coor.load.limit}{Maximum number of coordinates to load.}

\item{genome.load.limit}{Maximum number of genome data to load.}

\item{genome.fasta.style}{Fasta file style for genome data.}

\item{genome.ncbi.db}{NCBI database for genome data.}

\item{use.UCSC.chr.name}{Boolean indicating if UCSC chromosome naming is used.}

\item{verbose}{Boolean indicating if verbose output is enabled.}
}
\value{
Output directory containing exploration plots.
}
\description{
Function generates various exploratory analyses.
}
