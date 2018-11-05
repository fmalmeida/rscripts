#!/usr/bin/Rscript
# Setting Help
'usage: render_report.R [--input=<file> --prefix=<chr> --fasta=<chr> --gff=<file> --gffDir=<chr> --gbk=<chr> --command=<chr> --config=<file> --DBice=<true_or_false> --iceSummary=<file> --iceGff=<file> --isResfams=<true_or_false> --resfamsSummary=<file> --resfamsGff=<file>]

options:
  --input=<file>                  Rmd file to render
  --prefix=<chr>                  Prefix used
  --fasta=<chr>                   Fasta File Used
  --gff=<file>                    GFF file produced
  --gffDir=<chr>                  GFF file directory
  --gbk=<chr>                     Genbank File produced
  --command=<chr>                 Command executed
  --config=<file>                 Configuration file used
  --DBice=<boolean>    [Default: TRUE]
  --iceSummary=<file>            ices summary table
  --iceGff=<file>                ices gff file
  --isResfams=<true_or_false>    [Default: TRUE]
  --resfamsSummary=<file>        resfams summary table
  --resfamsGff=<file>            resfams gff file' -> doc

# Parse parameters
suppressMessages(library(docopt))
opt <- docopt(doc)

rmarkdown::render(opt$input, 
                  params = list( 
  prefix = opt$prefix, 
  fasta = opt$fasta, 
  gff = opt$gff, 
  gbk = opt$gbk, 
  command = opt$command, 
  config = opt$config, 
  gff_dir = opt$gffDir,
  ice_summary = opt$iceSummary, 
  ice_gff = opt$iceGff, 
  resfams_summary = opt$resfamsSummary, 
  resfams_gff = opt$resfamsGff
  ))
