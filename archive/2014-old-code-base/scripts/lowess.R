#!/TL/opt/bin/Rscript

# R script to perform lowess smoothing of signal tracks
# expects a line-based file of the form
# chrom x-values y-values
# and writes the result in bedgraph format back to the same file
# chrom start end value


options(warn=2)
sink(type="message")
options(scipen=1000000)
library(multicore)

smooth.signal <- function(regionstring, fraction, iter, res)
{
	tmp <- unlist(strsplit(regionstring, split="\t", fixed=TRUE))
	chrom <- tmp[1]
	xvals <- as.numeric(unlist(strsplit(tmp[2], split=" ", fixed=TRUE)))
	yvals <- as.numeric(unlist(strsplit(tmp[3], split=" ", fixed=TRUE)))
	yvals.smooth <- unlist((lowess(xvals, yvals, f=fraction, iter=iter))$y)
	yvals.smooth[yvals.smooth <= 0] <- 0.00001
	result.regions <- c()
	start <- xvals[1]
	data_pos <- res %/% 2
	for (i in c(1:(length(xvals) / res)))
	{
		end <- start + res
		value <- round(yvals.smooth[data_pos], 3)
		region <- paste(chrom, start, end, value, sep="\t")
		result.regions <- c(result.regions, region)
		start <- end
		data_pos <- data_pos + res
	}
	outstring <- paste(result.regions, collapse="\n")
	return (outstring)
}

args <- commandArgs(trailingOnly=TRUE)
FILENAME <- args[1]
RESOLUTION <- as.numeric(args[2])
NCPU <- as.numeric(args[3])
FRACTION <- as.numeric(args[4])
ITERATION <- as.numeric(args[5])
con.in <- file(FILENAME, open="rt")
rough.regions <- readLines(con.in, n=-1, ok=TRUE)
close(con.in)
smooth.regions <- mclapply(rough.regions, smooth.signal, fraction=FRACTION, iter=ITERATION, res=RESOLUTION, mc.cores=NCPU, mc.preschedule=TRUE, mc.silent=TRUE, mc.cleanup=TRUE)
#smooth.regions <- lapply(rough.regions, smooth.signal, fraction=FRACTION, iter=ITERATION, res=RESOLUTION)
con.out <- file(FILENAME, open="wt")
writeLines(unlist(smooth.regions), con.out, sep="\n")
close(con.out)
warnings()
quit(save="no", status=0)
