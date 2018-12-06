#!/usr/bin/Rscript
suppressMessages(library(RSQLite))
suppressMessages(library(glue))
suppressMessages(library(stringr))
suppressMessages(library(DataCombine))

# Setting Help
'usage: parse_sql.R [--input=<file> --gff=<file> --out=<chr> --database=<chr> --type=<chr> --scoverage=<int>]
options:
  -i, --input=<file>    sqlite database outputed
  -s, --start=<int>     OPTIONAL: retrieve elements from this start position.
  -e, --end=<int>       OPTIONAL: retrieve elements until this end position.
  -f, --fofn=<file>     OPTIONAL: retrieve elements based on ids from fofn file.
  -t, --type=<chr>      Type of FASTA to output: nt|aa|both [default: both]
  -p, --prefix=<chr>    Output prefix
  -d, --outdir=<chr>    Output directory. Default: Current directory

  At least one of the OPTIONAL parameters must be supplied. The others are required.
  How the script will treat the OPTIONAL arguments:
      * Start position alone will retrieve all elements which start position is
  greater than the value given;
      * End position alone will retrieve all elements which end position is less
  than the value given;
      * Start and end together will retrieve all elements located between that 
  genomic range;
      * The fofn file containing one element id (the first attribute of the last gff
  column) per line will retrieve the elements that have this id.' -> doc

# Parse parameters
suppressMessages(library(docopt))
opt <- docopt(doc)

if (is.null(opt$input)){
  stop("At least one argument must be supplied (input file)\n", call.=FALSE)
}

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
get_from_ids <- function(y) {
    dbGetQuery(con, paste("SELECT * FROM ProteinFasta WHERE ID='", y, "'", sep=""))
  }
### Function to extract elements by region
region = function(s, e) {
  dbGetQuery(con, paste("SELECT * FROM FinalGFF WHERE start > ", s, " AND end < ", e, sep=""))
}

## Loading SQL database driver
drv <- dbDriver("SQLite")
dbname <- file.path("~/Downloads/canu_pacbio_only/sqlDB/archaea_cynthia.sqlite")
con <- dbConnect(drv, dbname=dbname)

## Getting Data Out
head(dbListTables(con)) # Lists tables in database
dbListFields(con, "FinalGFF") # Lists table fields
dbListFields(con, "NucleotideFasta") # Lists table fields
dbListFields(con, "ProteinFasta") # Lists table fields

### Searching elements that have the pattern resistance in feature column
out <- 
  dbGetQuery(con, "SELECT * FROM FinalGFF WHERE feature LIKE '%resistance%'")

### Searching elements by region
suppressWarnings(res <- region(100000, 200000))

### Get FASTA based on ids
ids <- c("cgdabima_00001", "cgdabima_00003", "cgdabima_00005")
ids <- toupper(ids)

suppressWarnings(res <- lapply(ids, get_from_ids))
out <- suppressWarnings(data.frame(t(sapply(res,c))))
out_fasta <- 
  glue('>{out$ID} {out$Comment}\n{out$Sequence}')

### Get FASTA based on region
suppressWarnings(res <- region(100000, 200000))
ids <- getAttributeField(res$attributes, "id", ";")
ids <- toupper(ids)
suppressWarnings(res <- lapply(ids, get_from_ids))
out <- suppressWarnings(data.frame(t(sapply(res,c))))
out_fasta <- 
  glue('>{out$ID} {out$Comment}\n{out$Sequence}')
