#!/usr/bin/python3.1

import fnmatch as fnm
import gzip as gz
import collections as col
import sys as sys
import os as os


def count_gzipped_files():
	path = "/TL/epigenetics2/archive00/pebert/epimodule/module"
	all_files = os.listdir(path)
	gz_files = fnm.filter(all_files, "*.gz")
	count_dict = col.defaultdict(int)

	for gzfile in gz_files:
		infile =  gz.open(os.path.join(path, gzfile), "rb")
		for line in infile:
			if line[0] == ">":
				cols = line.strip().split("\t")
				ids = cols[2]
				amount = len(ids.split(","))
				count_dict[amount] += 1
			continue
		infile.close()

	with open("counting.txt" ,"w") as outfile:
		summed = sum(count_dict.values())
		outfile.write("Total number of motifs:\t%s\n" %summed)
		for k,v in count_dict.iteritems():
			outfile.write("Number of %s-combinations:\t%s\n" %(k,v))


sys.exit(0)
