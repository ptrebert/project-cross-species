import unittest as unittest
import collections as col
import io as io


from modules.datastructs import RegionCluster
class TestRegionCluster(unittest.TestCase):
	def setUp(self):
		# region_x = {'seq_chrom': , 'seq_start': , 'seq_end': , 'reg_chrom': , 'reg_id': , 'conservation': , 'output_line': }
		region_a = {'seq_chrom':'chr1' , 'seq_start':100 , 'seq_end':210 , 'reg_chrom':'chr11' , 'reg_id':110 , 'conservation':0.8 , 'output_line':'chr1\t100\t210\thg19@chr11@110' }
		region_b = {'seq_chrom':'chr2' , 'seq_start':200 , 'seq_end':290 , 'reg_chrom':'chr12' , 'reg_id':112 , 'conservation':0.9 , 'output_line':'chr2\t200\t290\thg19@chr12@112' }
		self.test_cluster_a = RegionCluster(100, region_a)
		self.test_cluster_b = RegionCluster(100, region_b)
	def tearDown(self):
		del self.test_cluster_a
		del self.test_cluster_b
	def test_cluster_not_overlapping(self):
		self.assertFalse(self.test_cluster_a.overlaps(self.test_cluster_b))
		self.assertFalse(self.test_cluster_b.overlaps(self.test_cluster_a))
	def test_cluster_self_overlapping(self):
		self.assertTrue(self.test_cluster_a.overlaps(self.test_cluster_a))
	def test_cluster_self_overlap_length(self):
		self.assertEqual(self.test_cluster_a.overlap_length(self.test_cluster_a), 110)
	def test_cluster_overlapping(self):
		region_c = {'seq_chrom':'chr1' , 'seq_start':209 , 'seq_end':309 , 'reg_chrom':'chr11' , 'reg_id':111 , 'conservation':0.7 , 'output_line':'chr1\t209\t309\thg19@chr11@111' }
		test_cluster = RegionCluster(100, region_c)
		self.assertTrue(self.test_cluster_a.overlaps(test_cluster))
		self.assertTrue(test_cluster.overlaps(self.test_cluster_a))
	def test_cluster_overlap_length_right(self):
		region_c = {'seq_chrom':'chr1' , 'seq_start':209 , 'seq_end':309 , 'reg_chrom':'chr11' , 'reg_id':111 , 'conservation':0.7 , 'output_line':'chr1\t209\t309\thg19@chr11@111' }
		test_cluster = RegionCluster(100, region_c)
		self.assertEqual(self.test_cluster_a.overlap_length(test_cluster), 1) # 0-based, half open
		self.assertEqual(test_cluster.overlap_length(self.test_cluster_a), 1) # 0-based, half open
	def test_cluster_overlap_length_left(self):
		region_c = {'seq_chrom':'chr1' , 'seq_start':5 , 'seq_end':100 , 'reg_chrom':'chr11' , 'reg_id':109 , 'conservation':0.7 , 'output_line':'chr1\t5\t100\thg19@chr11@111' }
		test_cluster = RegionCluster(100, region_c)
		self.assertEqual(self.test_cluster_a.overlap_length(test_cluster), 0) # 0-based, half open
		self.assertEqual(test_cluster.overlap_length(self.test_cluster_a), 0) # 0-based, half open
	def test_cluster_overlap_length_both_right(self):
		region_c = {'seq_chrom':'chr1' , 'seq_start':55 , 'seq_end':105 , 'reg_chrom':'chr11' , 'reg_id':109 , 'conservation':0.7 , 'output_line':'chr1\t5\t100\thg19@chr11@111' }
		test_cluster_a = RegionCluster(100, region_c)
		region_d = {'seq_chrom':'chr2' , 'seq_start':280 , 'seq_end':290 , 'reg_chrom':'chr12' , 'reg_id':113 , 'conservation':0.7 , 'output_line':'chr2\t280\t290\thg19@chr12@113' }
		test_cluster_b = RegionCluster(100, region_d)
		self.assertEqual(self.test_cluster_a.overlap_length(test_cluster_a), 5) # 0-based, half open
		self.assertEqual(test_cluster_a.overlap_length(self.test_cluster_a), 5) # 0-based, half open
		self.assertEqual(self.test_cluster_b.overlap_length(test_cluster_b), 10) # 0-based, half open
		self.assertEqual(test_cluster_b.overlap_length(self.test_cluster_b), 10) # 0-based, half open
	def test_cluster_overlap_length_both_left(self):
		region_c = {'seq_chrom':'chr1' , 'seq_start':55 , 'seq_end':210 , 'reg_chrom':'chr11' , 'reg_id':109 , 'conservation':0.7 , 'output_line':'chr1\t55\t210\thg19@chr11@111' }
		test_cluster_a = RegionCluster(100, region_c)
		region_d = {'seq_chrom':'chr2' , 'seq_start':200 , 'seq_end':320 , 'reg_chrom':'chr12' , 'reg_id':113 , 'conservation':0.7 , 'output_line':'chr2\t200\t320\thg19@chr12@113' }
		test_cluster_b = RegionCluster(100, region_d)
		self.assertEqual(self.test_cluster_a.overlap_length(test_cluster_a), 110) # 0-based, half open
		self.assertEqual(test_cluster_a.overlap_length(self.test_cluster_a), 110) # 0-based, half open
		self.assertEqual(self.test_cluster_b.overlap_length(test_cluster_b), 90) # 0-based, half open
		self.assertEqual(test_cluster_b.overlap_length(self.test_cluster_b), 90) # 0-based, half open
	def test_cluster_overlapping_real(self):
		region_a = {'seq_chrom':'chr6' , 'seq_start':118472965 , 'seq_end':118474490 , 'reg_chrom':'chr18' , 'reg_id':148514 , 'conservation':0.7 , 'output_line':"ham" }
		region_b = {'seq_chrom':'chr6' , 'seq_start':118473471 , 'seq_end':118474420 , 'reg_chrom':'chr18' , 'reg_id':142212 , 'conservation':0.5 , 'output_line':"eggs" }
		cluster_a = RegionCluster(100, region_a)
		cluster_b = RegionCluster(100, region_b)
		self.assertTrue(cluster_a.overlaps(cluster_b))
		self.assertTrue(cluster_b.overlaps(cluster_a))
	def test_cluster_boundary_overlap_right(self):
		region_c = {'seq_chrom':'chr1' , 'seq_start':210 , 'seq_end':300 , 'reg_chrom':'chr11' , 'reg_id':111 , 'conservation':0.7 , 'output_line':'chr1\t210\t300\thg19@chr11@111' }
		test_cluster = RegionCluster(100, region_c)
		self.assertFalse(self.test_cluster_a.overlaps(test_cluster))
		self.assertFalse(test_cluster.overlaps(self.test_cluster_a))
	def test_cluster_boundary_overlap_left(self):
		region_c = {'seq_chrom':'chr1' , 'seq_start':10 , 'seq_end':100 , 'reg_chrom':'chr11' , 'reg_id':109 , 'conservation':0.9 , 'output_line':'chr1\t10\t100\thg19@chr11@109' }
		test_cluster = RegionCluster(100, region_c)
		self.assertFalse(self.test_cluster_a.overlaps(test_cluster))
		self.assertFalse(test_cluster.overlaps(self.test_cluster_a))
	def test_gapless_adjacent_right(self):
		region_c = {'seq_chrom':'chr1' , 'seq_start':210 , 'seq_end':300 , 'reg_chrom':'chr11' , 'reg_id':111 , 'conservation':0.8 , 'output_line':'chr1\t210\t300\thg19@chr11@111' }
		self.assertTrue(self.test_cluster_a.is_adjacent(region_c))
	def test_gapless_adjacent_left(self):
		region_c = {'seq_chrom':'chr1' , 'seq_start':0 , 'seq_end':100 , 'reg_chrom':'chr11' , 'reg_id':109 , 'conservation':0.8 , 'output_line':'chr1\t0\t100\thg19@chr11@109' }
		self.assertTrue(self.test_cluster_a.is_adjacent(region_c))
	def test_gap_adjacent_left(self):
		region_c = {'seq_chrom':'chr1' , 'seq_start':3 , 'seq_end':95 , 'reg_chrom':'chr11' , 'reg_id':109 , 'conservation':0.8 , 'output_line':'chr1\t3\t95\thg19@chr11@109' }
		self.assertTrue(self.test_cluster_a.is_adjacent(region_c))
	def test_gap_adjacent_right(self):
		region_c = {'seq_chrom':'chr1' , 'seq_start':243 , 'seq_end':367 , 'reg_chrom':'chr11' , 'reg_id':111 , 'conservation':0.8 , 'output_line':'chr1\t243\t367\thg19@chr11@111' }
		self.assertTrue(self.test_cluster_a.is_adjacent(region_c))
	def test_only_seq_adjacent(self):
		region_c = {'seq_chrom':'chr1' , 'seq_start':243 , 'seq_end':367 , 'reg_chrom':'chr8' , 'reg_id':234 , 'conservation':0.8 , 'output_line':'chr1\t243\t367\thg19@chr8@234' }
		self.assertFalse(self.test_cluster_a.is_adjacent(region_c))
	def test_only_reg_adjacent(self):
		region_c = {'seq_chrom':'chr3' , 'seq_start':56 , 'seq_end':143 , 'reg_chrom':'chr11' , 'reg_id':111 , 'conservation':0.8 , 'output_line':'chr3\t56\t143\thg19@chr11@111' }
		self.assertFalse(self.test_cluster_a.is_adjacent(region_c))
	def test_add_region_left(self):
		region_c = {'seq_chrom':'chr1' , 'seq_start':5 , 'seq_end':95 , 'reg_chrom':'chr11' , 'reg_id':109 , 'conservation':0.8 , 'output_line':'chr1\t5\t95\thg19@chr11@109' }
		self.test_cluster_a.add_region(region_c)
		self.assertEqual(self.test_cluster_a.seq_chrom, 'chr1')
		self.assertEqual(self.test_cluster_a.seq_start, 5)
		self.assertEqual(self.test_cluster_a.seq_end, 210)
		self.assertEqual(self.test_cluster_a.reg_ids, set([109, 110]))
		self.assertEqual(self.test_cluster_a.cons_values, [0.8, 0.8])
		self.assertEqual(self.test_cluster_a.output, ['chr1\t100\t210\thg19@chr11@110', 'chr1\t5\t95\thg19@chr11@109'])
		self.assertEqual(self.test_cluster_a.get_cluster_regions(), "chr1\t100\t210\thg19@chr11@110\nchr1\t5\t95\thg19@chr11@109")
	def test_add_region_right(self):
		region_c = {'seq_chrom':'chr1' , 'seq_start':210 , 'seq_end':320 , 'reg_chrom':'chr11' , 'reg_id':111 , 'conservation':0.9 , 'output_line':'chr1\t210\t320\thg19@chr11@111' }
		self.test_cluster_a.add_region(region_c)
		self.assertEqual(self.test_cluster_a.seq_chrom, 'chr1')
		self.assertEqual(self.test_cluster_a.seq_start, 100)
		self.assertEqual(self.test_cluster_a.seq_end, 320)
		self.assertEqual(self.test_cluster_a.reg_ids, set([111, 110]))
		self.assertEqual(self.test_cluster_a.cons_values, [0.8, 0.9])
		self.assertEqual(self.test_cluster_a.output, ['chr1\t100\t210\thg19@chr11@110', 'chr1\t210\t320\thg19@chr11@111'])
		self.assertEqual(self.test_cluster_a.get_cluster_regions(), "chr1\t100\t210\thg19@chr11@110\nchr1\t210\t320\thg19@chr11@111")
