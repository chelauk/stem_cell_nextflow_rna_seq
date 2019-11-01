#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(QoRTs)

decoder <- read.table(args, header = T, stringsAsFactors = FALSE)
res <- read.qc.results.data(infile.dir = './', decoder = decoder, autodetectMissingSamples = T)
makeMultiPlot.colorBySample(res, outfile.dir = "./", plot.device.name = "pdf",  chromosome.name.style = "ENSEMBL", sequencing.type = "RNA", insertSize.plot.xlim = c(0,600))
