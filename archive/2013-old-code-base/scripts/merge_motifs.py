#!/usr/bin/python2.6

import sys as sys
import os as os
import fnmatch as fn
import cStringIO as sio
import gzip as gz
import collections as col

inpath = "/TL/epigenetics2/archive00/pebert/epimodule/module"
outpath = "/TL/epigenetics2/work/pebert/data/crm"

all_files = os.listdir(inpath)
gz_files = fn.filter(all_files, "*.gz")
outbuffers = col.defaultdict(sio.StringIO)

for gzfile in gz_files:
	name_parts = gzfile.split(".")
	file_number = name_parts[4]
	infile = gz.open(os.path.join(inpath, gzfile), "rb")
	region_id = None
	for line in infile:
		if line == "\n": continue
		if line[0] == ">":
			cols = line.strip().split("\t")
			tf_names = "|".join(cols[1].strip("()").split(","))
			tf_ids = "|".join(cols[2].strip("()").split(","))
			region_id = "@".join(["hg18", file_number, tf_names, tf_ids])
			continue
		else:
			cols = line.strip().split("\t")
			positions = []
			chrom = cols[1]
			for coord in cols[2:]:
				positions.append(int(coord.strip("+- ")))
			start = (min(positions)) - 50
			end = (max(positions)) + 50
			outline = "\t".join([chrom, str(start), str(end), region_id])
			outbuffers[chrom].write(outline + "\n")
			continue

for chrom, stream in outbuffers.iteritems():
	outf_name = "_".join(["hg18_crm", chrom, "Cai2010.bed"])
	outfile = open(os.path.join(outpath, outf_name), "w")
	outfile.write(stream.getvalue())
	outfile.close()

sys.exit(0)	
