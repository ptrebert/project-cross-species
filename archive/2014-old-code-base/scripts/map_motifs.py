#!/usr/bin/python2.6

import multiprocessing as mp
import os as os
import sys as sys
import re as re
import fnmatch as fn


inpath = "/TL/epigenetics2/work/pebert/data/crm"
outpath = "/TL/epigenetics2/work/pebert/data/crm"
chainfile_susscr2 = "/TL/epigenetics2/work/datasets/chain/hg19/hg19ToSusScr2.over.chain.gz"
liftover = "/TL/epigenetics2/archive00/bin/ucsctools/liftOver -minMatch=0.1 "
fnregexp = re.compile("hg19_crm_chr[0-9]{1,2}_Cai2010.bed")
all_files = os.listdir(inpath)
matched_files = [ os.path.join(inpath, name) for name in all_files if fnregexp.search(name) is not None ]

def perform_lift(filename):
	newname = filename.replace("hg19", "susscr2")
	command = liftover + filename + " " + chainfile_susscr2 + " " + newname + " " + "/dev/null"
	os.system(command)
	return True

workers = mp.Pool(len(matched_files))
workers.map(perform_lift, matched_files)
workers.close()
workers.join()
sys.exit(0)
