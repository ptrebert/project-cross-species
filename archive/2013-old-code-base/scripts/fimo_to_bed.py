#!/usr/bin/python3.1

"""Script to convert txt output of FIMO into BED
format - read from stdin, write to stdout
"""

import sys as sys

while 1:
	line = sys.stdin.readline()
	if not line:
		break
	parts = line.strip().split("\t")
	# no qvalue is computed since output is set to text
	#pattern name   sequence name   start   stop    strand  score   p-value q-value matched sequence
	outline = "\t".join([parts[1], parts[2], parts[3], parts[4], parts[0], parts[8], parts[5], parts[6]]) # parts[7] would be the qvalue
	sys.stdout.write(outline + "\n")
	

#sys.stderr.write("\nRemember to re-format the header of the output file\n")
sys.exit(0)
