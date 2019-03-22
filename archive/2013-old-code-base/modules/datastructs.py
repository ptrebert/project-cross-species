"""Module to hold special purpose
data structures or classes
"""

import random as rand

class RegionCluster:
    """ To resolve conflicting mappings with overlapping regions
    An entry has the form <sequence position> in target species
    followed by <region id> from source species
    Need to make sure that both positional information pertain to the
    intuition of a neighbouring region (w/o gaps)
    """
    def __init__(self, res, region):
        self.resolution = int(res)
        self.seq_chrom = region['seq_chrom']
        self.seq_start = int(region['seq_start'])
        self.seq_end = int(region['seq_end'])
        self.reg_chrom = region['reg_chrom']
        self.reg_ids = set([region['reg_id']])
        self.cons_values = [region['conservation']]
        self.output = [region['output_line']]
        self.cluster_cons = 0
    def overlaps(self, cluster):
        """Checks for possible overlaps of two clusters based on
        the genomic position in the target species (overlapping mappings)
        """
        if self.seq_chrom != cluster.seq_chrom: return False
        if self.seq_start < cluster.seq_end and self.seq_end > cluster.seq_start: return True
        return False
    def overlap_length(self, cluster):
        assert self.seq_chrom == cluster.seq_chrom, "Called RegionCluster.overlap_length with clusters from different chromosomes"
        # difficult question: does liftOver produce 0-based, half open intervals?
        # Assume 'yes' for calculation of region length
        if self.seq_start == cluster.seq_start and self.seq_end == cluster.seq_end: return (self.seq_end - self.seq_start)
        elif cluster.seq_end >= self.seq_start > cluster.seq_start: return (cluster.seq_end - self.seq_start)
        elif self.seq_start <= cluster.seq_end and self.seq_end <= cluster.seq_end: return (self.seq_end - cluster.seq_start)
        elif self.seq_start <= cluster.seq_start and self.seq_end > cluster.seq_end: return (cluster.seq_end - cluster.seq_start)
        elif self.seq_start > cluster.seq_start and self.seq_end < cluster.seq_end: return (self._seq_end - self.seq_start)
        return 0
    def is_adjacent(self, region):
        """Check if region is adjacent to the cluster
        and should be added; based on both positional
        information (source and target species)
        """
        if self.reg_chrom != region['reg_chrom'] or self.seq_chrom != region['seq_chrom']: return False
        if ((region['reg_id'] - 1) in self.reg_ids) or ((region['reg_id'] + 1) in self.reg_ids):
            # allow minor gaps in sequence continuity since the mapping results in target regions of various size
            # However, make sure no overlaps occur
            # Check right end - reminder: 0-base, half open intervals
            if self.seq_end <= int(region['seq_start']) and abs(self.seq_end - int(region['seq_start'])) <= self.resolution: return True
            # Check left end of cluster
            if self.seq_start >= int(region['seq_end']) and abs(self.seq_start - int(region['seq_end'])) <= self.resolution: return True
        return False
    def add_region(self, region):
        assert self.reg_chrom == region['reg_chrom'], "RegionCluster: adding source region from different chromosome: %s - region is: %s" %(self.reg_chrom, region['output_line'])
        assert self.seq_chrom == region['seq_chrom'], "RegionCluster: adding target region from different chromosome: %s - region is: %s" %(self.seq_chrom, region['output_line'])
        self.reg_ids.add(region['reg_id'])
        self.cons_values.append(region['conservation'])
        self.output.append(region['output_line'])
        # is_adjacent guarantees that region is not overlapping
        if int(region['seq_start']) >= self.seq_end:
            self.seq_end = int(region['seq_end'])
        else:
            self.seq_start = int(region['seq_start'])
        return      
    def compute_conservation(self):
        """compute_conservation is called when all regions
        are added to a cluster
        """
        self.cluster_cons = sum(self.cons_values) / len(self.cons_values)
        return
    def get_cluster_regions(self):
        out = "\n".join(self.output)
        return out
    def summary(self):
        output = "\nThis is cluster %s\nChrom: %s\nStart: %s\nEnd: %s\nWill print: %s\n" %(id(self), self.seq_chrom, self.seq_start, self.seq_end, "\t".join(self.output))
        return output

class RepConTree:
    def __init__(self, region, relax):
        self._left_child = None
        self._right_child = None
        self._bound_low = region['perc_Rep'] - relax
        self._bound_high = region['perc_Rep'] + relax
        self._relax = relax
        self._regions = [region]
    def insert_region(self, region):
        if region['perc_Rep'] < self._bound_low:
            if self._left_child:
                self._left_child.insert_region(region)
            else:
                self._left_child = RepConTree(region, self._relax)
        elif region['perc_Rep'] > self._bound_high:
            if self._right_child:
                self._right_child.insert_region(region)
            else:
                self._right_child = RepConTree(region, self._relax)
        else:
            self._regions.append(region)
        return
    def find_match(self, query_region):
        if query_region['perc_Rep'] < self._bound_low:
            if self._left_child:
                return self._left_child.find_match(query_region)
            return None
        elif query_region['perc_Rep'] > self._bound_high:
            if self._right_child:
                return self._right_child.find_match(query_region)
            return None
        else:
            return (rand.sample(self._regions, 1))[0]

class GCTree:
    def __init__(self, region, relax):
        self._left_child = None
        self._right_child = None
        self._bound_low = region['perc_GC'] - relax
        self._bound_high = region['perc_GC'] + relax
        self._relax = relax
        self._repcontree = RepConTree(region, relax)
    def insert_region(self, region):
        if region['perc_GC'] < self._bound_low:
            if self._left_child:
                self._left_child.insert_region(region)
            else:
                self._left_child = RepConTree(region, self._relax)
        elif region['perc_GC'] > self._bound_high:
            if self._right_child:
                self._right_child.insert_region(region)
            else:
                self._right_child = RepConTree(region, self._relax)
        else:
            self._repcontree.insert_region(region)
        return
    def find_match(self, query_region):
        if query_region['perc_GC'] < self._bound_low:
            if self._left_child:
                return self._left_child.find_match(query_region)
            return None
        elif query_region['perc_GC'] > self._bound_high:
            if self._right_child:
                return self._right_child.find_match(query_region)
            return None
        else:
            return self._repcontree.find_match(query_region)

class SizeTree:
    def __init__(self, region, relax):
        self._left_child = None
        self._right_child = None
        self._bound_low = region['length'] * (1. - relax / 100)
        self._bound_high = region['length'] * (1. + relax / 100)
        self.relax = relax
        self.num_regions = 1
        self._gctree = GCTree(region, relax)
    def insert_region(self, region):
        self.num_regions += 1
        if region['length'] < self._bound_low:
            if self._left_child:
                self._left_child.insert_region(region)
            else:
                self._left_child = SizeTree(region, self.relax)
        elif region['length'] > self._bound_high:
            if self._right_child:
                self._right_child.insert_region(region)
            else:
                self._right_child = SizeTree(region, self.relax)
        else:
            self._gctree.insert_region(region)
        return
    def find_match(self, query_region):
        if query_region['length'] < self._bound_low:
            if self._left_child:
                return self._left_child.find_match(query_region)
            return None
        elif query_region['length'] > self._bound_high:
            if self._right_child:
                return self._right_child.find_match(query_region)
            return None
        else:
            return self._gctree.find_match(query_region)

class OverlapTree:
    def __init__(self, level=None, region=None):
        self._left_child = None
        self._right_child = None
        self._level = level
        self._region = region
        self._subtree = None
    def is_overlapping(self, region):
        if self._level == "chromosome":
            if region['chrom'] < self._region['chrom']:
                if self._left_child: return self._left_child.is_overlapping(region)
                self._left_child = OverlapTree(level="chromosome", region=region)
                self._left_child._subtree = OverlapTree(level="nucleotide", region=region)
                return False
            if region['chrom'] > self._region['chrom']:
                if self._right_child: return self._right_child.is_overlapping(region)
                self._right_child = OverlapTree("chromosome", region=region)
                self._right_child._subtree = OverlapTree(level="nucleotide", region=region)
                return False
            return self._subtree.is_overlapping(region)         
        if self._level == "nucleotide":
            if region['end'] < self._region['start']:
                if self._left_child: return self._left_child.is_overlapping(region)
                self._left_child = OverlapTree(level="nucleotide", region=region)
                return False
            if region['start'] > self._region['end']:
                if self._right_child: return self._right_child.is_overlapping(region)
                self._right_child = OverlapTree(level="nucleotide", region=region)
                return False
            return True
        self._level = "chromosome"
        self._region = region
        self._subtree = OverlapTree(level="nucleotide", region=region)
        return False
    def collect_matches(self):
        if self._level == "chromosome":
            ret = self._subtree.collect_matches()
            if self._left_child:
                ret.extend(self._left_child.collect_matches())
            if self._right_child:
                ret.extend(self._right_child.collect_matches())
            return ret
        if self._level == "nucleotide":
            ret = [self._region]
            if self._left_child:
                ret.extend(self._left_child.collect_matches())
            if self._right_child:
                ret.extend(self._right_child.collect_matches())
            return ret
        raise ValueError("OverlapTree level was found to be {}".format(self._level))
