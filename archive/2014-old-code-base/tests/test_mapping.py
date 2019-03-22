import unittest as unittest
import collections as col
import io as io


from modules.mapping import Mapping
from modules.datastructs import RegionCluster
class TestMappingProcessBlank(unittest.TestCase):
	def setUp(self):
		self.clusters = self.test_build_clusters()
	def tearDown(self):
		del self.clusters
	def test_build_clusters(self):
		region1 = { 'seq_chrom':'chr1', 'seq_start':10 , 'seq_end':23 , 'reg_chrom':'chr18', 'reg_id':101, 'conservation':0.8, 'output_line':'region101' }	
		region2 = { 'seq_chrom':'chr1', 'seq_start':23 , 'seq_end':35 , 'reg_chrom':'chr18', 'reg_id':102, 'conservation':0.8, 'output_line':'region102' }	
		region3 = { 'seq_chrom':'chr1', 'seq_start':37 , 'seq_end':48 , 'reg_chrom':'chr18', 'reg_id':103, 'conservation':0.8, 'output_line':'region103' }	
		region4 = { 'seq_chrom':'chr1', 'seq_start':48 , 'seq_end':58 , 'reg_chrom':'chr18', 'reg_id':104, 'conservation':0.8, 'output_line':'region104' }	
		region5 = { 'seq_chrom':'chr1', 'seq_start':56 , 'seq_end':70 , 'reg_chrom':'chr18', 'reg_id':205, 'conservation':0.9, 'output_line':'region205' }	
		region6 = { 'seq_chrom':'chr1', 'seq_start':70 , 'seq_end':81 , 'reg_chrom':'chr18', 'reg_id':206, 'conservation':0.9, 'output_line':'region206' }	
		region7 = { 'seq_chrom':'chr1', 'seq_start':83 , 'seq_end':95 , 'reg_chrom':'chr18', 'reg_id':207, 'conservation':0.9, 'output_line':'region207' }	
		region8 = { 'seq_chrom':'chr1', 'seq_start':95 , 'seq_end':107 , 'reg_chrom':'chr18', 'reg_id':208, 'conservation':0.9, 'output_line':'region208' }	
		region9 = { 'seq_chrom':'chr1', 'seq_start':17 , 'seq_end':27 , 'reg_chrom':'chr18', 'reg_id':309, 'conservation':0.6, 'output_line':'region309' }	
		region10 = { 'seq_chrom':'chr1', 'seq_start':27 , 'seq_end':40 , 'reg_chrom':'chr18', 'reg_id':310, 'conservation':0.6, 'output_line':'region310' }	
		region11 = { 'seq_chrom':'chr1', 'seq_start':42 , 'seq_end':57 , 'reg_chrom':'chr18', 'reg_id':311, 'conservation':0.6, 'output_line':'region311' }	
		region12 = { 'seq_chrom':'chr1', 'seq_start':57 , 'seq_end':69 , 'reg_chrom':'chr18', 'reg_id':312, 'conservation':0.6, 'output_line':'region312' }
		region13 = { 'seq_chrom':'chr1', 'seq_start':63 , 'seq_end':74 , 'reg_chrom':'chr18', 'reg_id':413, 'conservation':0.5, 'output_line':'region413' }
		region14 = { 'seq_chrom':'chr1', 'seq_start':74 , 'seq_end':90 , 'reg_chrom':'chr18', 'reg_id':414, 'conservation':0.5, 'output_line':'region414' }
		region15 = { 'seq_chrom':'chr1', 'seq_start':97 , 'seq_end':110 , 'reg_chrom':'chr18', 'reg_id':415, 'conservation':0.5, 'output_line':'region415' }
		region16 = { 'seq_chrom':'chr1', 'seq_start':108 , 'seq_end':118 , 'reg_chrom':'chr18', 'reg_id':516, 'conservation':0.3, 'output_line':'region516' }
		region17 = { 'seq_chrom':'chr1', 'seq_start':118 , 'seq_end':130 , 'reg_chrom':'chr18', 'reg_id':517, 'conservation':0.3, 'output_line':'region517' }
		region18 = { 'seq_chrom':'chr1', 'seq_start':132 , 'seq_end':145 , 'reg_chrom':'chr18', 'reg_id':518, 'conservation':0.3, 'output_line':'region518' }
		ovl_regions = [region1, region2, region3, region4, region5, region6, region7, region8, region9, region10, \
						region11, region12, region13, region14, region15, region16, region17, region18]
		clusters = Mapping(None, None)._build_clusters(ovl_regions, 10)
		self.assertEqual(len(clusters), 5)
		new_dict = {}
		for c in clusters.itervalues():
			if c.seq_start == 10:
				self.assertEqual(c.reg_ids, set([101, 102, 103, 104]))
				new_dict[100] = c
			elif c.seq_start == 56:
				self.assertEqual(c.reg_ids, set([205, 206, 207, 208]))
				new_dict[200] = c
			elif c.seq_start == 17:
				self.assertEqual(c.reg_ids, set([309, 310, 311, 312]))
				new_dict[300] = c
			elif c.seq_start == 63:
				self.assertEqual(c.reg_ids, set([413, 414, 415]))
				new_dict[400] = c
			else:
				self.assertEqual(c.reg_ids, set([516, 517, 518]))
				new_dict[500] = c
		return new_dict
	def test_compute_overlap(self):
		conslst, ovl_dict = Mapping(None, None)._compute_overlaps(self.clusters)
		# dont care right now
		#self.assertAlmostEqual(conslst, [ (0.9, 200), (0.8, 100), (0.6, 300), (0.5, 400), (0.3, 500) ])
		#self.assertEqual(conslst, [(1,500), (2,100), (3,200), (3,300), (3,400)])
		comp = sorted([(k,set(v)) for k,v in ovl_dict.iteritems() ])
		self.assertEqual(comp, [ (100,set([200, 300])) , (200, set([100, 300, 400])), (300, set([100, 200, 400])), (400, set([200, 300, 500])), (500, set([400]))])
	def test_select_clusters(self):
		num_ovl, ovl_dict = Mapping(None, None)._compute_overlaps(self.clusters)
		selected = Mapping(None, None)._select_clusters(num_ovl, ovl_dict, self.clusters)
		self.assertEqual(set([500, 200]), selected)

class TestMapFiles(unittest.TestCase):
	def setUp(self):
		self.testObj = Mapping({"degree":3}, None)
	def tearDown(self):
		del self.testObj
	def test_compute_consecutive_regions(self):
		
		result = self.testObj._compute_regions("chr1", [100,101,102,103,104,105,106,107,108,109,110], [0.8, 0.8, 0.8, 0.8, 0.8, 0.5, 0.5, 0.5, 0.5, 0.5], 5)
		self.assertEqual(result,"chr1\t100\t105\t0.8\n")
	def test_extract_region(self):
		pass
	def test_compute_regions_lovl(self):
		pass
	def test_compute_regions_rovl(self):
		pass
