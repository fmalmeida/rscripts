#!/usr/bin/Rscript
'usage: gff2sql.R [--input=<file> --out=<chr> --nucleotide=<file> --aminoacid=<file>]

options:
  -i, --input=<file>         GFF file to transform in SQL
  -o, --out=<chr>            SQL database name to output [default: out.sql]
  -n, --nucleotide=<file>    Takes in the nucleotide FASTA.
  -a, --aminoacid=<file>     Takes in the protein FASTA' -> doc
library("docopt")
opt <- docopt(doc)
suppressMessages(library(RSQLite))

## Loading SQL database driver
drv <- dbDriver("SQLite")
dbname <- file.path(".", opt$out)
con <- dbConnect(drv, dbname=dbname)

## Loading GFF file
data <- read.delim(opt$input, header = FALSE, stringsAsFactors = FALSE)

### Give data a header
names(data) <- c("chr", "source", "feature", "start", "end", 
                 "score", "strand", "frame", "attributes")

## Create SQL table to store GFF data
suppressWarnings(dbGetQuery(con, "CREATE Table FinalGFF (chr TEXT, source TEXT, feature TEXT, 
           start INTEGER, end INTEGER, score INTEGER, strand TEXT, 
           frame INTEGER, attributes TEXT)"))

## Create sql rule
sql <- "INSERT INTO FinalGFF VALUES ($chr, $source, $feature, 
$start, $end, $score, $strand, $frame, $attributes)"

## Open db
suppressWarnings(dbBegin(con))

## Send rule
res <- suppressWarnings(dbSendQuery(con,sql))

# Insert data based on rule
suppressWarnings(dbBind(res, data))
suppressWarnings(dbFetch(res))
suppressWarnings(dbClearResult(res))
### Close db
suppressWarnings(dbCommit(con))

## Loading Protein fasta
genes_aa <- read.delim(opt$aminoacid, header = FALSE)
names(genes_aa) <- c("ID", "Comment", "Sequence")

## Create SQL table to store Protein FASTA
suppressWarnings(dbGetQuery(con, "CREATE Table ProteinFasta (ID TEXT, Comment TEXT, Sequence TEXT)"))

## Create sql rule
sql <- "INSERT INTO ProteinFasta VALUES ($ID, $Comment, $Sequence)"

## Open db
suppressWarnings(dbBegin(con))

## Send rule
res <- suppressWarnings(dbSendQuery(con,sql))

## Insert data based on rule
suppressWarnings(dbBind(res, genes_aa))
suppressWarnings(dbFetch(res))
suppressWarnings(dbClearResult(res))

## Close db
suppressWarnings(dbCommit(con))

## Loading Nucleotide fasta
genes_nt <- read.delim(opt$nucleotide, header = FALSE)
names(genes_nt) <- c("ID", "Comment", "Sequence")

## Create SQL table to store Nucleotide FASTA
suppressWarnings(dbGetQuery(con, "CREATE Table NucleotideFasta (ID TEXT, Comment TEXT, Sequence TEXT)"))

## Create sql rule
sql <- "INSERT INTO NucleotideFasta VALUES ($ID, $Comment, $Sequence)"

## Open db
suppressWarnings(dbBegin(con))

## Send rule
res <- suppressWarnings(dbSendQuery(con,sql))

## Insert data based on rule
suppressWarnings(dbBind(res, genes_nt))
suppressWarnings(dbFetch(res))
suppressWarnings(dbClearResult(res))

## Close db
suppressWarnings(dbCommit(con))