#!/usr/bin/env Rscript

suppressMessages(library(optparse))
suppressMessages(library(GenomicRanges))
suppressMessages(library(LOLA))

LOLA_CORE_PATH = "/TL/deep/fhgfs/projects/pebert/thesis/refdata/lola/LOLACore"
LOLA_EXT_PATH = "/TL/deep/fhgfs/projects/pebert/thesis/refdata/lola/LOLAExt"
LOLA_CUSTOM_PATH = "/TL/deep/fhgfs/projects/pebert/thesis/refdata/lola/LOLACustom"


cmd_option_list = list(make_option(c("-a", "--assembly"), type="character", default=NULL),
                       make_option(c("-i", "--input"), type="character", default=NULL),
                       make_option(c("-u", "--universe"), type="character", default=NULL),
                       make_option(c("-c", "--cores"), type="numeric", default=2),
                       make_option(c("-d", "--database"), type="character", default="custom"),
                       make_option(c("-o", "--output"), type="character", default=NULL)
                       )
opt_parser = OptionParser(option_list=cmd_option_list)
cmd_opts = parse_args(opt_parser)

if (is.null(cmd_opts$assembly) | is.null(cmd_opts$input))
{
    stop('Assembly and input file have to be specified')
}

load_paths = list(core=LOLA_CORE_PATH, ext=LOLA_EXT_PATH, custom=LOLA_CUSTOM_PATH)
regdb_load_path = paste(load_paths[[cmd_opts$database]], cmd_opts$assembly, sep='/')

regiondb = loadRegionDB(regdb_load_path)

#col_idx = which(regiondb$regionAnno$collection == cmd_opts$universe, arr.ind=TRUE)
#
#cells = unique(regiondb$regionAnno[col_idx, ]$cellType)
#
#print(cells)
#
#quit(save="no", status=1)

user_roi = readBed(cmd_opts$input)

universe = unlist(regiondb$regionGRL[which(regiondb$regionAnno$collection == cmd_opts$universe)])
universe = disjoin(universe)

results = runLOLA(user_roi, universe, regiondb, cores=cmd_opts$cores)

writeCombinedEnrichment(results, outFolder= "/home/pebert/temp/creepiest/cons_genes/mmu_lola.tsv", includeSplits=FALSE)

quit(save="no", status=0)
