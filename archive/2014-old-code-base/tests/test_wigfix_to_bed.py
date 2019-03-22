#import path_include

import unittest as unittest
import collections as col
import io as io

# this module is being tested
import scripts.conv_wigfix_to_bed as wfbed


class TestWigFixToBed(unittest.TestCase):
	def setUp(self):
		ConvConfig = col.namedtuple("ConvConfig", "inputfile outputfile resolution aggregate build")
		configuration = ConvConfig(inputfile="INFILE", outputfile="OUTFILE", resolution=5, aggregate="AVG", build="hg19")
		self.test_object = wfbed.WigFixToBed(configuration)
	def tearDown(self):
		del self.test_object
	def test_regular_start(self):
		str_in = "fixedStep chrom=chr1 start=100 step=1\n5\n1\n1\n1\n1\n1\n7\n8\n9\n"
		str_out = "chr1\t100\t105\thg19@chr1@20@1.0\n"
		sio_out = self.test_object._convert(io.StringIO(str_in), io.StringIO(), (lambda lst: sum(lst)/len(lst)))
		self.assertEqual(str_out, sio_out.getvalue())
	def test_negative_start(self):
		str_in = "fixedStep chrom=chr1 start=101 step=1\n5\n1\n1\n1\n1\n7\n8\n9\n6\n"
		str_out = "chr1\t100\t105\thg19@chr1@20@1.8\n"
		sio_out = self.test_object._convert(io.StringIO(str_in), io.StringIO(), (lambda lst: sum(lst)/len(lst)))
		self.assertEqual(str_out, sio_out.getvalue())
	def test_irregular_start(self):
		str_in = "fixedStep chrom=chr1 start=103 step=1\n1\n1\n1\n1\n1\n5\n1\n1\n1\n"
		str_out = "chr1\t105\t110\thg19@chr1@21@1.8\n"
		sio_out = self.test_object._convert(io.StringIO(str_in), io.StringIO(), (lambda lst: sum(lst)/len(lst)))
		self.assertEqual(str_out, sio_out.getvalue())
	def test_continuous_blocks(self):
		str_in = "fixedStep chrom=chr1 start=107 step=1\n107\n108\n109\n110\n1\n1\n1\n1\n1\n2\nfixedStep chrom=chr1 start=117 step=1\n2\n2\n2\n2\n121\n122\n"
		str_out = "chr1\t110\t115\thg19@chr1@22@1.0\nchr1\t115\t120\thg19@chr1@23@2.0\n"
		sio_out = self.test_object._convert(io.StringIO(str_in), io.StringIO(), (lambda lst: sum(lst)/len(lst)))
		self.assertEqual(str_out, sio_out.getvalue())
	def test_discontinuous_chroms(self):
		str_in = "fixedStep chrom=chr1 start=107 step=1\n107\n108\n109\n110\n1\n1\n1\n1\n1\n2\nfixedStep chrom=chr12 start=117 step=1\n117\n118\n119\n120\n2\n2\n2\n2\n2\n"
		str_out = "chr1\t110\t115\thg19@chr1@22@1.0\nchr12\t120\t125\thg19@chr12@24@2.0\n"
		sio_out = self.test_object._convert(io.StringIO(str_in), io.StringIO(), (lambda lst: sum(lst)/len(lst)))
		self.assertEqual(str_out, sio_out.getvalue())
	def test_discontinuous_blocks(self):
		str_in = "fixedStep chrom=chr1 start=107 step=1\n107\n108\n109\n110\n1\n1\n1\n1\n1\n116\nfixedStep chrom=chr1 start=151 step=1\n5\n5\n5\n5\n5\n6\n6\n6\n6\n6\n161\n"
		str_out = "chr1\t110\t115\thg19@chr1@22@1.0\nchr1\t150\t155\thg19@chr1@30@5.0\nchr1\t155\t160\thg19@chr1@31@6.0\n"
		sio_out = self.test_object._convert(io.StringIO(str_in), io.StringIO(), (lambda lst: sum(lst)/len(lst)))
		self.assertEqual(str_out, sio_out.getvalue())
