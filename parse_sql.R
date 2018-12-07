#!/usr/bin/Rscript
suppressMessages(library(RSQLite))
suppressMessages(library(glue))
suppressMessages(library(stringr))
suppressMessages(library(DataCombine))

# Setting Help
'usage: parse_sql.R [--input=<file> --start=<int> --end=<int> --fofn=<file> --regex=<chr> --type=<chr> --prefix=<chr> --outdir=<chr>]
options:
  -i, --input=<file>    sqlite database outputed
  -s, --start=<int>     OPTIONAL: retrieve elements from this start position.
  -e, --end=<int>       OPTIONAL: retrieve elements until this end position.
  -f, --fofn=<file>     OPTIONAL: retrieve elements based on ids from fofn file.
  -r, --regex=<chr>     OPTIONAL: retrieve elements based on searching a pattern in specific column of GFF file. Example: feature|resistance will search for elements that \
                        have the pattern resistance in the feature column. GFF Columns: chr|source|feature|attributes.
  -t, --type=<chr>      Type of FASTA to output: nt|aa|both [default: both]
  -p, --prefix=<chr>    Output prefix
  -d, --outdir=<chr>    Output directory. Default: Current directory

  At least one of the OPTIONAL parameters must be supplied. The others are required.
  The script will treat the OPTIONAL arguments as follows:
      * Start position alone will retrieve all elements which start position is greater than the value given;
      * End position alone will retrieve all elements which end position is less than the value given;
      * Start and end together will retrieve all elements located between that genomic range;
      * The fofn file containing one element id (the first attribute of the last gff column) per line will retrieve the elements that have this id.' -> doc

# Parse parameters
suppressMessages(library(docopt))
opt <- docopt(doc)

if (is.null(opt$start) && is.null(opt$end) && is.null(opt$fofn) && is.null(opt$regex)){
  stop("At least one of the OPTIONAL parameters (start|end|fofn) must be supplied\n", call.=FALSE)
}

if (is.null(opt$input)){
  stop("At least one argument must be supplied (input file)\n", call.=FALSE)
}

if (is.null(opt$outdir)) {
  outdir <- "."
} else { outdir <- opt$outdir }

## Creating functions
### To retrieve attribute information
getAttributeField <- function (x, field, attrsep = ";") { 
  s = strsplit(x, split = attrsep, fixed = TRUE) 
  sapply(s, function(atts) { 
    a = strsplit(atts, split = "=", fixed = TRUE) 
    m = match(field, sapply(a, "[", 1)) 
    if (!is.na(m)) { rv = a[[m]][2] 
    } 
    else { 
      rv = as.character(NA) 
    } 
    return(rv) 
  }) 
}
### To extract FASTA based on IDs
get_from_ids <- function(y, x) {
    dbGetQuery(con, paste("SELECT * FROM ", x, " WHERE ID='", y, "'", sep=""))
}
###
get_gff_from_ids <- function(y) {
  dbGetQuery(con, paste("SELECT * FROM FinalGFF WHERE attributes LIKE '%", y, "%'", sep=""))
}

## Loading SQL database driver
drv <- dbDriver("SQLite")
dbname <- file.path(opt$input)
con <- dbConnect(drv, dbname=dbname)

## Getting Data Out
#head(dbListTables(con)) # Lists tables in database
#dbListFields(con, "FinalGFF") # Lists table fields
#dbListFields(con, "NucleotideFasta") # Lists table fields
#dbListFields(con, "ProteinFasta") # Lists table fields

### Searching elements that have the pattern resistance in feature column
if (!is.null(opt$regex)) {
  #### Subset GFF
  list    <- strsplit(opt$regex, "|", fixed = TRUE)
  string  <- glue("SELECT * FROM FinalGFF WHERE {list[[1]][1]} LIKE '%{list[[1]][2]}%'")
  out     <- suppressWarnings(dbGetQuery(con, string))
  outname <- glue("{outdir}/subseted.gff")
  write.table(out, file = outname, quote = FALSE, sep = "\t", 
              col.names = TRUE, row.names = FALSE)
}

### Searching elements by genomic region
#### Only with start position
if (!is.null(opt$start) && is.null(opt$end)) {
  ##### Subset GFF
  start <- opt$start
  string <- glue("SELECT * FROM FinalGFF WHERE start > {start}")
  out     <- suppressWarnings(dbGetQuery(con, string))
  outname <- glue("{outdir}/subseted.gff")
  write.table(out, file = outname, quote = FALSE, sep = "\t", 
              col.names = TRUE, row.names = FALSE)
} else if (is.null(opt$start) && !is.null(opt$end)) {
  ##### Subset GFF
  end <- opt$end
  string <- glue("SELECT * FROM FinalGFF WHERE end < {end}")
  out     <- suppressWarnings(dbGetQuery(con, string))
  outname <- glue("{outdir}/subseted.gff")
  write.table(out, file = outname, quote = FALSE, sep = "\t", 
              col.names = TRUE, row.names = FALSE)
} else if (!is.null(opt$start) && !is.null(opt$end)) {
  ##### Subset GFF
  start <- opt$start
  end <- opt$end
  string <- glue("SELECT * FROM FinalGFF WHERE start > {start} AND end < {end}")
  out     <- suppressWarnings(dbGetQuery(con, string))
  outname <- glue("{outdir}/subseted.gff")
  write.table(out, file = outname, quote = FALSE, sep = "\t", 
              col.names = TRUE, row.names = FALSE)
}

### Searching elements by ids
if (!is.null(opt$fofn)) {
  fil <- file("ids.fofn")
  ids <- readLines(fil, n = -1)
  ids <- toupper(ids)
  suppressWarnings(res <- lapply(ids, get_gff_from_ids))
  out <- suppressWarnings(data.frame(t(sapply(res,c))))
  df <- apply(out,2,as.character)
  outname <- glue("{outdir}/subseted.gff")
  write.table(df, file = outname, quote = FALSE, sep = "\t", 
              col.names = TRUE, row.names = FALSE)
}

### Subset FASTA
if (opt$type == 'both') {
  # Get aa FASTA
  ids <- getAttributeField(as.character(out$attributes), "id", ";")
  ids <- toupper(ids)
  suppressWarnings(res <- lapply(ids, get_from_ids, x = "ProteinFasta"))
  out <- suppressWarnings(data.frame(t(sapply(res,c))))
  out_fasta <- 
    glue('>{out$ID} {out$Comment}\n{out$Sequence}')
  outname <- glue("{outdir}/subseted_aa.fasta")
  write(out_fasta, file = outname, sep = "\n")
  # Get nt FASTA
  suppressWarnings(res <- lapply(ids, get_from_ids, x = "NucleotideFasta"))
  out <- suppressWarnings(data.frame(t(sapply(res,c))))
  out_fasta <- 
    glue('>{out$ID} {out$Comment}\n{out$Sequence}')
  outname <- glue("{outdir}/subseted_nt.fasta")
  write(out_fasta, file = outname, sep = "\n")
} else if (opt$type == nt) {
  # Get nt FASTA
  ids <- getAttributeField(out$attributes, "id", ";")
  ids <- toupper(ids)
  suppressWarnings(res <- lapply(ids, get_from_ids, x = "NucleotideFasta"))
  out <- suppressWarnings(data.frame(t(sapply(res,c))))
  out_fasta <- 
    glue('>{out$ID} {out$Comment}\n{out$Sequence}')
  outname <- glue("{outdir}/subseted_nt.fasta")
  write(out_fasta, file = outname, sep = "\n")
} else if (opt$type == aa) {
  # Get aa FASTA
  ids <- getAttributeField(out$attributes, "id", ";")
  ids <- toupper(ids)
  suppressWarnings(res <- lapply(ids, get_from_ids, x = "ProteinFasta"))
  out <- suppressWarnings(data.frame(t(sapply(res,c))))
  out_fasta <- 
    glue('>{out$ID} {out$Comment}\n{out$Sequence}')
  outname <- glue("{outdir}/subseted_aa.fasta")
  write(out_fasta, file = outname, sep = "\n")
}