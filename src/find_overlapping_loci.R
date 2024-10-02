#!/usr/bin/env Rscript

library(data.table)
library(GenomicRanges)
library(rtracklayer)

#############
# FUNCTIONS #
#############

decode_gff_dt <- function(gff_dt) {
  my_dt <- copy(gff_dt)

  col_order <- names(my_dt)
  char_cols <- col_order[
    my_dt[, sapply(.SD, class)] == "character"
  ]

  my_dt[
    , (char_cols) := lapply(.SD, utils::URLdecode),
    .SDcols = char_cols
  ]

  return(my_dt)
}

fread_gff <- function(
    gff_file) {
  my_gff <- import(gff_file)

  my_dt <- data.table(data.frame(my_gff))
  factor_cols <- names(my_dt)[sapply(my_dt, is.factor)]
  my_dt[, (factor_cols) := lapply(.SD, as.character), .SDcols = factor_cols]

  return(my_dt)
}

generate_seqinfo_from_fai <- function(fai_file, genome_label) {
  fai_data <- fread(fai_file, header = FALSE)
  return(
    fai_data[
      , Seqinfo(
        seqnames = V1, seqlengths = V2, genome = rep(genome_label, .N)
      )
    ]
  )
}

make_granges_from_captus_gff <- function(
    gff_file, seqinfo_object) {
  gff_dt <- fread_gff(gff_file = gff_file)
  gff_dt_decoded <- decode_gff_dt(gff_dt)
  gff_granges <- makeGRangesFromDataFrame(
    gff_dt_decoded,
    keep.extra.columns = TRUE,
    seqinfo = seqinfo_object,
    na.rm = TRUE
  )
  return(gff_granges)
}

###########
# GLOBALS #
###########



# the query gff
query_gff_file <- "output/010_captus/qos.peakall/min1000000/03_extractions/qos.1000000__captus-ext/01_coding_NUC/NUC_contigs.gff"

# the reference, i.e. usually mega353
ref_gff_file <- "output/010_captus/qos.mega353/min1000000/03_extractions/qos.1000000__captus-ext/01_coding_NUC/NUC_contigs.gff"

fai_file <- "test/file.fasta.fai"
genome_label <- "qos"

# second test with pzijinensis
query_gff_file <- "output/010_captus/pzijinensis.peakall/min1000000/03_extractions/pzijinensis.1000000__captus-ext/01_coding_NUC/NUC_contigs.gff"
ref_gff_file <- "output/010_captus/pzijinensis.mega353/min1000000/03_extractions/pzijinensis.1000000__captus-ext/01_coding_NUC/NUC_contigs.gff"
fai_file <- "test/pzijinensis.1000000.fasta.fai"
genome_label <- "pzijinensis"


########
# MAIN #
########

seqinfo_object <- generate_seqinfo_from_fai(
  fai_file,
  genome_label
)

query_gff <- make_granges_from_captus_gff(
  query_gff_file,
  seqinfo_object = seqinfo_object
)

ref_gff <- make_granges_from_captus_gff(
  ref_gff_file,
  seqinfo_object = seqinfo_object
)

# Filter for 'protein_match' records
query_loci <- query_gff[query_gff$type == "protein_match:NUC"]
ref_loci <- ref_gff[ref_gff$type == "protein_match:NUC"]

# Find overlapping protein_match records
overlaps <- findOverlaps(
  query = query_loci,
  subject = ref_loci, 
  ignore.strand = TRUE
)

proximate_loci <- findOverlaps(
  query = query_loci,
  subject = ref_loci, 
  ignore.strand = TRUE,
  maxgap= 10e3
)


