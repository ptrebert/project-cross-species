#!/usr/bin/python2.6

import sys as sys
import cStringIO as sio

class TransfacMatrix(object):
	def __init__(self, tfid, altid, alen, values):
		self.motif_id = tfid
		self.motif_alt = altid
		self.motif_alength = alen
		self.motif_length = len(values)
		self.motif_values = values
	def print_base_prob(self, base):
		sys.stdout.write("\nProbabilities for base %s in motif:\n" %base)
		for t in self.motif_values:
			sys.stdout.write(str(self.motif_values.index(t)) + "\t")
			sys.stdout.write(t[base] + "\n")
		return
	def MEME_output(self):
		out = "MOTIF %s %s\nletter-probability matrix: alength= %s w= %s\n" %(self.motif_id, self.motif_alt, self.motif_alength, self.motif_length)
		matrix = map(lambda d: d["A"] + " " + d["C"] + " " + d["G"] + " " + d["T"], self.motif_values)
		out += "\n".join(matrix)
		return out


def parse_transfac_files(path):
	transfac_dict = {}
	with open(path, "r") as infile:
		motif_id = None
		motif_alt = None
		motif_alphlen = None
		motif_values = []
		for line in infile:
			if line == "\n":
				new_tf_motif = TransfacMatrix(motif_id, motif_alt, motif_alphlen, motif_values)
				transfac_dict[motif_id] = new_tf_motif
				continue
			if line[0] == ">":
				cols = line.strip().split("\t")
				motif_id = cols[0].strip(">")
				motif_alt = cols[-1]
				motif_values = []				
				continue
			else:
				cols = line.strip().split(" ")
				assert len(cols) == 4, "Not 4 probs found in line %s" %line
				motif_alphlen = len(cols)
				newline = { "A":cols[0], "C":cols[1], "G":cols[2], "T":cols[3] }
				motif_values.append(newline)
				continue
		new_tf_motif = TransfacMatrix(motif_id, motif_alt, motif_alphlen, motif_values)
		transfac_dict[motif_id] = new_tf_motif
	return transfac_dict
						

infile_loc = "/TL/epigenetics2/archive00/pebert/epimodule/crm_oct.28/data/normalized_transfac_matrices"

tfd = parse_transfac_files(infile_loc)

outbuffer = sio.StringIO()
outbuffer.write("MEME version 4\n\nALPHABET= ACGT\n\nBackground letter frequencies\nA 0.295 C 0.205 G 0.205 T 0.295\n\n")
for mat in tfd.itervalues():
	outbuffer.write(mat.MEME_output())
	outbuffer.write("\n\n")

outfile = open("Transfac_MEME.txt", "w")
outfile.write(outbuffer.getvalue())
outfile.close()
sys.exit(0)
