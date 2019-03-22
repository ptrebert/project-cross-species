"""Module to handle preprocessing or conversion
that is only optional (data conversion, annotation of peak regions...)
Should remain useable as imported module as well as usecase for the pipeline
(e.g. for batch conversion)
"""

import collections as col
import logging as log
import gzip as gz
import subprocess as sp
import os as os
import tempfile as tf

import modules.osinterface as osi

class Toolbox(object):
    def __init__(self, configuration):
        self.config = configuration
        self.tmp_delete = []
        self.logger = log.getLogger(__name__)
    def stand_alone(self):
        """Calling this function means the toolbox
        is used like a pipeline usecase - configuration
        needs to be sufficient in that case
        """
        return True
    def clean_up(self):
        """Try to delete all tempfiles created
        Errors are reported to logfile
        """
        status = True
        self.logger.debug("Attempting to remove all {} temp files".format(len(self.tmp_delete)))
        for entry in self.tmp_delete:
            try:
                os.remove(entry)
            except OSError:
                self.logger.warning("Given path {} might be a directory or is no longer existing".format(entry))
                status = False
            except Exception as e:
                self.logger.warning("Caught exception {} while removing temp files".format(e))
                status = False
        return status
    def converter(self, which):
        pass
    def ant_sequence(self, region_lst, genome_location, filterN):
        """Given a list of region dictionaries, annotate each region with the
        correct sequence (FASTA files are assumed to be in path 'genome_location')      
        """
        self.logger.debug("Annotating genomic regions with sequences")
        allseqfiles = osi.collect_all_files(genome_location)
        chromseqfiles = {}
        for f in allseqfiles:
            chrom = os.path.basename(f).strip(".fagz")
            chromseqfiles[chrom] = f
        split_to_chrom = col.defaultdict(list)
        for rd in region_lst:
            split_to_chrom[rd['chrom']].append(rd)
        assert set(split_to_chrom.keys()).issubset(set(chromseqfiles.keys())), "Cannot find all appropriate sequence files for chromosomes {}\nFiles: {}".format(split_to_chrom.keys(), chromseqfiles.keys())
        annotated = []
        for c, rlst in split_to_chrom.items():
            chromfile = chromseqfiles[c]
            seq = self.prs_fasta_file(chromfile)
            for regiond in rlst:
                # this is supposed to annotate BED regions, so the correct sequence is 1-based and includes the last position
                regseq = seq[ int(regiond['start']) + 1 : int(regiond['end']) + 1 ]
                if filterN: # this is a float; None or 0.0 will be False
                    seqlen = len(regseq)
                    countN = float(regseq.count("N"))
                    percN = (countN / seqlen) * 100.
                    if percN > filterN: continue
                regiond['seq'] = regseq
                annotated.append(regiond)
        self.logger.debug("Annotation complete, returning from normal operation")
        return annotated
    def prs_chrom_sizes(self, dbchromsizes, genome_assembly, chrom_matcher):
        """Parses a file containing chromosome sizes for a specific assembly
        Returns a tuple list with tuples (chrom, size)
        chrom_matcher is used to restrict accepted chromosomes
        """
        try:
            dbentry = dbchromsizes[genome_assembly]
            chromsizes = []
            with open(os.path.join(dbentry['dir'], dbentry['filename']), 'r') as f:
                for line in f:
                    cols = line.strip().split('\t')
                    if chrom_matcher(cols[0]):
                        try:
                            chromsizes.append((cols[0], int(cols[1])))
                        except ValueError:
                            self.logger.warning("Found malformed line > {} < in chromsize file {}".format(line, f.name))
                            continue
            return chromsizes
        except KeyError:
            raise Exception("Could not find chromosome sizes file for assembly {}".format(genome_assembly))
        return None
    def prs_fasta_file(self, seqfile):
        """Read a FASTA sequence file (.fa or .fa.gz) and return a single
        sequence string
        """
        self.logger.debug("Parsing FASTA file: {}".format(seqfile))
        chromfile = None
        if ".fa.gz" in seqfile:
            chromfile = gz.open(seqfile, "rb")
        else:
            chromfile = open(seqfile, "r")
        # For python3x, str with utf-8 encoding is required since read() returns a byte object
        seq = chromfile.read()
        if isinstance(seq, bytes):
            seq = str(seq, "utf-8")
        chromfile.close()
        if seq[0] == ">":
            # regular header in fasta file
            idx = seq.find("\n")
            seq = seq[idx: ].replace("\n", "")
        seq = seq.replace("\n", "")
        self.logger.debug("FASTA file parsed, returning a sequence of length {}".format(len(seq)))
        return seq
    def prs_genome_region_file(self, infile, ncol, lkey, ltype, rtype):
        """Read a file containing genomic regions (line sep is '\n', col sep is '\t')
        Default is 3 columns (chrom, start, end) and list of tuples as return value
        rtype can be set to dict [return list of dicts], in case of 4 columns, need to specify key for last column
        and possibe coercion function ltype (like str, float)
        """
        ret = []
        with open(infile, "r") as inf:
            for line in inf:
                cols = line.strip().split("\t")
                cols = tuple(cols[:ncol])
                ret.append(cols)
        keys = []
        if ncol == 3:
            keys = ["chrom", "start", "end"]            
            ret = [ (t[0],int(t[1]),int(t[2])) for t in ret ]
        else:
            keys = ["chrom", "start", "end", lkey]
            ret = [ (t[0],int(t[1]),int(t[2]),ltype(t[3])) for t in ret ]
        if rtype == "dict":
            ret = [ dict(zip(keys, t)) for t in ret ]
        return ret
    def prs_bed_wheader_wid(self, filename):
        """Parse the given BED file with header
        Set ID for each entry, return a list of dicts
        """
        retlst = []
        idcount = 0
        with open(filename, 'r') as infile:
            header = infile.readline().strip('\n# ').split('\t')
            for line in infile:
                idcount += 1
                entry = dict(zip(header, line.strip().split('\t')))
                entry['id'] = str(idcount)
                entry['length'] = str(int(entry['end']) - int(entry['start']))
                retlst.append(entry)
        return retlst, idcount
    def wrt_str_tempfile(self, outstr):
        """Write a string to a tempfile and return the file object
        """
        tmpf = tf.NamedTemporaryFile(mode="r+", buffering=-1, suffix=".tmp", prefix="Py_", dir="/tmp", delete=False)
        tmpf.write(outstr)
        tmpf.flush()
        os.fsync(tmpf.fileno())
        tmpf.close()
        self.tmp_delete.append(tmpf.name)
        return tmpf
    def scl_bed_complement(self, bedfile, genomefile, outfile):
        if outfile == "TEMPFILE":
            tmpf = tf.NamedTemporaryFile(mode="r+", buffering=-1, suffix=".tmp", prefix="Py_", dir="/tmp", delete=False)
            outfile = tmpf
            self.tmp_delete.append(tmpf.name)
        exe = self.config['bedcomplement']
        p = sp.Popen([exe, "-i", bedfile.name, "-g", genomefile.name], stdout=outfile)
        ret = p.communicate()
        return outfile
    def scl_bed_intersect_overlap(self, afile, bfile, outfile, issorted):
        if outfile == "TEMPFILE":
            tmpf = tf.NamedTemporaryFile(mode="r+", buffering=-1, suffix=".tmp", prefix="Py_", dir="/tmp", delete=False)
            outfile = tmpf
            self.tmp_delete.append(tmpf.name)
        exe = self.config['intersectbed']
        param_lst = [exe, "-a", afile, "-b", bfile, "-wo", "-sorted"] if issorted else [exe, "-a", afile, "-b", bfile, "-wo"]
        p = sp.Popen(param_lst, stdout=outfile)
        ret = p.communicate()
        return outfile          
    def scl_bed_intersect_count(self, afile, bfile, outfile, issorted):
        if outfile == "TEMPFILE":
            tmpf = tf.NamedTemporaryFile(mode="r+", buffering=-1, suffix=".tmp", prefix="Py_", dir="/tmp", delete=False)
            outfile = tmpf
            self.tmp_delete.append(tmpf.name)
        exe = self.config['intersectbed']
        param_lst = [exe, "-c", "-a", afile, "-b", bfile]
        p = sp.Popen(param_lst, stdout=outfile)
        ret = p.communicate()
        return outfile
    def scl_bed_sort(self, infile, outfile):
        if outfile == "TEMPFILE":
            tmpf = tf.NamedTemporaryFile(mode="r+", buffering=-1, suffix=".tmp", prefix="Py_", dir="/tmp", delete=False)
            outfile = tmpf
            self.tmp_delete.append(tmpf.name)
        exe = self.config['bedsort']
        param_lst = [exe, infile.name, outfile.name]
        p = sp.Popen(param_lst, stdout=outfile)
        ret = p.communicate()
        return outfile
