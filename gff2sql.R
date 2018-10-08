#!/usr/bin/Rscript

suppressMessages(library(RSQLite))
suppressMessages(library(optparse))

## Setting parameters
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="GFF file to be converted into a sql database", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="myNewDb.sqlite",
              help="SQL database output name [default= %default]", metavar="character"),
  make_option(c("-n", "--nucleotide"), type = "character", default=NULL,
              help="Nucleotide FASTA containing annotated genes [default= %default]", metavar="character"),
  make_option(c("-a", "--aminoacid"), type = "character", default=NULL,
              help="Protein FASTA containing annotated genes [default= %default]", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)\n", call.=FALSE)
}

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