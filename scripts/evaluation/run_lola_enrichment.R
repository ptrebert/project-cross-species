#!/usr/bin/env Rscript

suppressMessages(library(optparse))
suppressMessages(library(GenomicRanges))
suppressMessages(library(LOLA))

LOLA_CORE_PATH = "/TL/deep/fhgfs/projects/pebert/thesis/refdata/lola/LOLACore"
LOLA_EXT_PATH = "/TL/deep/fhgfs/projects/pebert/thesis/refdata/lola/LOLAExt"
LOLA_CUSTOM_PATH = "/TL/deep/fhgfs/projects/pebert/thesis/refdata/lola/LOLACustom"

LOLA_BASE_OUT = "/TL/deep/fhgfs/projects/pebert/thesis/projects/cross_species/processing/norm/task_lola"
LOLA_SUPP_OUT = "/TL/deep-external01/nobackup/pebert/cloudshare/mpiinf/phd/chapter_projects/crossspecies/supplement/lola"

cmd_option_list = list(make_option(c("-a", "--assembly"), type="character", default=NULL),
                       make_option(c("-i", "--input"), type="character", default=NULL),
                       make_option(c("-u", "--universe"), type="character", default='promoter'),
                       make_option(c("-c", "--cores"), type="numeric", default=2),
                       make_option(c("-d", "--database"), type="character", default="custom"),
                       make_option(c("-o", "--output"), type="character", default=LOLA_BASE_OUT),
                       make_option(c("-r", "--redefineUserSet"), type="logical", default=FALSE),
                       make_option(c("-s", "--suffix"), type="character", default=NULL))

opt_parser = OptionParser(option_list=cmd_option_list)
cmd_opts = parse_args(opt_parser)

if (is.null(cmd_opts$assembly) | is.null(cmd_opts$input) | is.null(cmd_opts$suffix))
{
    stop('Assembly, input file and folder suffix have to be specified')
}
load_paths = list(core=LOLA_CORE_PATH, ext=LOLA_EXT_PATH, custom=LOLA_CUSTOM_PATH)
regdb_load_path = paste(load_paths[[cmd_opts$database]], cmd_opts$assembly, sep='/')

regiondb = loadRegionDB(regdb_load_path)

if (cmd_opts$universe == 'promoter') {
    col_idx = which(regiondb$regionAnno$collection == 'genes' &
                    grepl('_prom.bed', regiondb$regionAnno$filename, fixed=TRUE),
                    arr.ind=TRUE)
} else if (cmd_opts$universe == 'body') {
    col_idx = which(regiondb$regionAnno$collection == 'genes' &
                    grepl('_body.bed', regiondb$regionAnno$filename, fixed=TRUE),
                    arr.ind=TRUE)
} else {
    stop('Unknown universe requested')
}

#quit(save="no", status=1)

user_roi = readBed(cmd_opts$input)

universe = unlist(regiondb$regionGRL[col_idx])
universe = disjoin(universe)

results = runLOLA(user_roi, universe, regiondb, cores=cmd_opts$cores, redefineUserSets=cmd_opts$redefineUserSet)

out_folder = paste(cmd_opts$output, paste(cmd_opts$assembly, cmd_opts$universe,
                                          cmd_opts$database, cmd_opts$suffix, sep='_'), sep='/')

writeCombinedEnrichment(results, outFolder=out_folder, includeSplits=FALSE)

out_file_default = paste(out_folder, 'allEnrichments.tsv', sep='/')

out_file_name = paste(cmd_opts$assembly, cmd_opts$universe, cmd_opts$database, cmd_opts$suffix, 'LOLA', sep='_')
out_file_name = paste(out_file_name, 'tsv', sep='.')

out_file_path = paste(LOLA_SUPP_OUT, out_file_name, sep='/')

foo = file.copy(out_file_default, out_file_path, overwrite=TRUE, recursive=FALSE)

quit(save="no", status=0)
