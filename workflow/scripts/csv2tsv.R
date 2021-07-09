#!/usr/bin/env Rscript

library(readr)

args = commandArgs(trailingOnly=TRUE)
infile = args[[1]]
outfile = args[[2]]

d = read_csv(infile, guess_max = 1e+6)
d

write_tsv(d, outfile)

