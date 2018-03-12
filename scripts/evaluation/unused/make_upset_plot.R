#!/usr/bin/env Rscript

suppressMessages(library(optparse))
suppressMessages(library(UpSetR))

UPSETR_BASE_OUT = "/TL/deep-external01/nobackup/pebert/cloudshare/mpiinf/phd/chapter_projects/crossspecies/supplement"

cmd_option_list = list(make_option(c("-i", "--input"), type="character", default=NULL),
                       make_option(c("-c", "--combinations"), type="character", default=NULL),
                       make_option(c("-p", "--path"), type="character", default=UPSETR_BASE_OUT),
                       make_option(c("-o", "--output"), type="character", default=NULL))

opt_parser = OptionParser(option_list=cmd_option_list)
cmd_opts = parse_args(opt_parser)

if (is.null(cmd_opts$input) | is.null(cmd_opts$combinations)) {
    stop('Input file required')
}

con = file(cmd_opts$input, open = "r")

plotList = list()
rows = 0
while (TRUE) {
    oneLine = readLines(con, n = 1, warn = FALSE)
    row = unlist(strsplit(oneLine, "\t"))
    if (length(row) < 1) {
        break
    }
    row.label = row[1]
    row.data = row[2:length(row)]
    plotList[[row.label]] = row.data
    rows = rows + 1
  }

close(con)

if (cmd_opts$combinations == "human") {
    plot.comb = list(list("ESC_ML", "ESC_Orth", "CD4_ML", "CD4_Orth", "Hepa_ML", "Hepa_Orth"),
                     list("ESC_ML", "CD4_ML", "Hepa_ML"),
                     list("ESC_Orth", "CD4_Orth", "Hepa_Orth"),
                     list("ESC_ML", "ESC_Orth"),
                     list("CD4_ML", "CD4_Orth"),
                     list("Hepa_ML", "Hepa_Orth"),
                     list("Hepa_ML"), list("CD4_ML"), list("ESC_ML"))
} else {
    plot.comb = list(list("ESC_ML", "ESC_Orth", "CD4_ML", "CD4_Orth", "Liver_ML", "Liver_Orth"),
                     list("ESC_ML", "CD4_ML", "Liver_ML"),
                     list("ESC_Orth", "CD4_Orth", "Liver_Orth"),
                     list("ESC_ML", "ESC_Orth"),
                     list("CD4_ML", "CD4_Orth"),
                     list("Liver_ML", "Liver_Orth"),
                     list("Liver_ML"), list("CD4_ML"), list("ESC_ML"))
}


plot.data = fromList(plotList)

png(paste(cmd_opts$path, cmd_opts$output, sep='/'), res=300, width=8, height=5, units="in")

upset(plot.data, order.by = "freq", text.scale = 1.5,
      nsets=nrows, scale.intersections = 'identity',
      intersections = plot.comb)

dev.off()

quit(save="no", status=0)
