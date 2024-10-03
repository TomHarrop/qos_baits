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
  gff_granges <- makeGRangesListFromDataFrame(
    gff_dt_decoded,
    keep.extra.columns = TRUE,
    seqinfo = seqinfo_object,
    na.rm = TRUE,
    split.field = "Name"
  )
  return(gff_granges)
}

process_hits <- function(overlaps, query_loci, ref_loci) {
  # all we care about is overlap so we want the extent of each hit
  query_hits <- as.data.table(range(query_loci[queryHits(overlaps)]))
  ref_hits <- as.data.table(range(ref_loci[subjectHits(overlaps)]))

  names(query_hits) <- paste("query_hit", names(query_hits), sep = "_")
  names(ref_hits) <- paste("ref_hit", names(ref_hits), sep = "_")

  return(
    cbind(ref_hits, query_hits)
  )
}


###########
# GLOBALS #
###########

# the query gff
query_gff_file <- "output/010_captus/qos.peakall/min1000000/03_extractions/qos.1000000__captus-ext/01_coding_NUC/NUC_contigs.gff"

# the reference, i.e. usually mega353
ref_gff_file <- "output/010_captus/qos.mega353/min1000000/03_extractions/qos.1000000__captus-ext/01_coding_NUC/NUC_contigs.gff"

# genome information
fai_file <- "test/file.fasta.fai"
genome_label <- "qos"

# second test with pzijinensis
query_gff_file <- "output/010_captus/pzijinensis.peakall/min1000000/03_extractions/pzijinensis.1000000__captus-ext/06_assembly_annotated/pzijinensis.1000000_hit_contigs.gff"
ref_gff_file <- "output/010_captus/pzijinensis.mega353/min1000000/03_extractions/pzijinensis.1000000__captus-ext/06_assembly_annotated/pzijinensis.1000000_hit_contigs.gff"
fai_file <- "test/pzijinensis.1000000.fasta.fai"
genome_label <- "pzijinensis"
ref_name <- "mega353"
query_name <- "peakall"
overlapping_loci_file <- "test/overlapping.csv"
proximate_loci_file <- "test/nearby.csv"

# order of output columns
output_cols <- c(
  "genome_name",
  "ref_name",
  "query_name",
  "ref_hit_group_name",
  "ref_hit_seqnames",
  "ref_hit_start",
  "ref_hit_end",
  "ref_hit_strand",
  "query_hit_group_name",
  "query_hit_seqnames",
  "query_hit_start",
  "query_hit_end",
  "query_hit_strand"
)


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
query_loci <- query_gff[sapply(query_gff, function(x) any(x$type == "protein_match:NUC"))]
ref_loci <- ref_gff[sapply(ref_gff, function(x) any(x$type == "protein_match:NUC"))]

# Not for GRangesList
# query_loci <- query_gff[query_gff$type == "protein_match:NUC"]
# ref_loci <- ref_gff[ref_gff$type == "protein_match:NUC"]

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
  maxgap = 10e3
)

# process the overlap information
full_hits_dt <- process_hits(overlaps, query_loci, ref_loci)
setorder(full_hits_dt, ref_hit_seqnames, ref_hit_start, ref_hit_width)
full_hits_dt[
  , c("ref_name", "query_name", "genome_name")
  := .(ref_name, query_name, genome_label)
]

proximate_hits_dt <- process_hits(proximate_loci, query_loci, ref_loci)
setorder(proximate_hits_dt, ref_hit_seqnames, ref_hit_start, ref_hit_width)
proximate_hits_dt[
  , c("ref_name", "query_name", "genome_name")
  := .(ref_name, query_name, genome_label)
]

# write the output
fwrite(full_hits_dt[, ..output_cols],
       overlapping_loci_file)

fwrite(proximate_hits_dt[, ..output_cols],
       overlapping_loci_file)

sessionInfo()
