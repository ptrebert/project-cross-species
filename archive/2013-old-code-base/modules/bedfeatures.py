"""This module is designed to compute features
for genomic regions in BED format
So, all data supplied to this class has the [minimal]
format chr - start - end - sequence
The module uses multicore processors since all computations
are assumed to be independent
"""
import os as os
import logging as log
import sys as sys
import time as ti
import multiprocessing as mp
import logging as log
import functools as fnt
import itertools as itr
import collections as col
import re as re
from queue import Empty as EmptyException
from queue import Full as FullException

class BEDFeatures(object):
    def __init__(self, features, sigperc, kmers):
        self.features = features
        self.kmers = kmers
        self.signal_perc = sigperc
        self.avail_features = {"gccpg":self._GC_CpG_content, "cgioverlap":self._CGI_overlap,
                                "kmerfrequency":self._kmer_frequency , "repcontent":self._repetitive_content,
                                "signaloverlap":self._signal_overlap , "crmmotifs":self._CRM_motifs ,
                                "tfbsmotifs":self._TFBS_motifs, "peakoverlap":self._peak_overlap }
        self.logger = log.getLogger(__name__)
        # the two following members are hardcoded since they are identical for all species
        # they just point to the files where the info about all possible CRM/TFBS features is stored
        self._crm_all_loc = "/TL/epigenetics2/work/pebert/data/crm/crm_2-12comb.txt"
        self._tfbs_all_loc = "/TL/epigenetics2/work/pebert/data/tfbs_meme/Transfac_Motifs.txt"
        self.crm = []
        self.tfbs = []
        if "tfbsmotifs" in features or "crmmotifs" in features:
            try:
                self._generate_all_info()
            except IOError:
                self.logger.error("Error parsing CRM/TFBS info, check file locations:\n{}\n{}".format(self._crm_all_loc, self._tfbs_all_loc))
                raise IOError("BEDFeatures object: cannot parse files")
        self.num_regions = 0
    def get_online_version(self):
        """Special version for serial feature computation
        Returns a closure
        """
        exec_functions = self._feature_functions()
        def compute_features(d):
            for f in exec_functions:
                d = f(d)
            return d
        return compute_features     
    def execute(self, bedregions, num_proc):
        self.num_regions = len(bedregions)
        self.logger.debug("Received {} regions for feature computation".format(len(bedregions)))
        assert self.num_regions > 0, "No genomic regions found for feature computation"
        mpman = mp.Manager()
        indata = mpman.Queue()
        active = mpman.Value("i", 0)
        active.value = num_proc
        lock = mpman.Lock()
        processed = mpman.Queue()
        to_exec = self._feature_functions()
        proc = [ mp.Process(target=self._compute_features, name=str(i), args=(to_exec, indata, processed, active, lock)) for i in range(num_proc) ]
        self.logger.info("Created {} child processes - starting computation now...".format(num_proc))
        for p in proc:
            p.daemon = True
            p.start()
            p.join(0)
        try: # accept dict or list
            for d in bedregions.values():
                indata.put(d)
        except AttributeError:
            for d in bedregions:
                indata.put(d)
        for i in range(num_proc):
            indata.put("SENTINEL") 
        results = []
        self.logger.debug("Start collecting results")
        for item in iter(processed.get, "DONE"):
            results.append(item)
        self.logger.info("Computation finished, checking for survivors...")
        for p in proc:
            if p.is_alive():
                p.join(0.1)
                p.terminate()
        self.logger.debug("Aggregated {}".format(len(results)))
        assert len(results) == self.num_regions, "Lost some genomic regions during feature calculation, multiprocessing.Queue corrupted?"
        return results
    def _generate_all_info(self):
        """Parse the two files with all possible inputs of CRM/TFBS
        """
        with open(self._crm_all_loc, "r") as crmfile:
            for line in crmfile:
                cols = line.split("\t")
                self.crm.append(cols[0])
        with open(self._tfbs_all_loc, "r") as tfbsfile:
            for line in tfbsfile:
                cols = line.split("\t")
                self.tfbs.append(cols[0])
        return
    def _compute_features(self, execFun, indata, processed, active, lock):
        for d in iter(indata.get, "SENTINEL"):
            #tmpd = d
            for f in execFun:
                #tmpd = f(tmpd)
                d = f(d)
            #processed.put(tmpd)
            processed.put(d)
        me = mp.current_process()
        nm = me.name
        pid = me.pid
        with lock:
            self.logger.debug("Process {} is checking out".format(pid))
            active.value -= 1
            self.logger.debug("Number of still active processes: {}".format(active.value))
            if active.value <= 0:
                self.logger.debug("Putting DONE signal to PROCESSED queue")
                processed.put("DONE")
        return
    def _feature_functions(self):
        comp_feat = self.features.lower().split("-")
        kmers = map(int, self.kmers.split("-"))
        to_exec = []
        for fn in comp_feat:
            try:
                if fn == "kmerfrequency":
                    for k in kmers:
                        pfn = fnt.partial(self.avail_features[fn], k)
                        to_exec.append(pfn)
                else:
                    pfn = fnt.partial(self.avail_features[fn])
                    to_exec.append(pfn)
            except KeyError:
                self.logger.warning("Feature {} is requested but not yet implemented".format(f))
                continue
        return to_exec
    def _GC_CpG_content(self, d):
        try:
            ham = d['perc_rawGC']
            eggs = d['perc_GC']
            bread = d['perc_CpG']
            return d
        except KeyError:
            seq = (d['seq']).lower()
            seqlen = float(len(seq))
            raw_GC = seq.count("c") + seq.count("g")
            total_CpG = seq.count("cg")
            total_GC = raw_GC - total_CpG * 2
            d['perc_rawGC'] = round(raw_GC/seqlen * 100., 3) 
            d['perc_GC'] = round(total_GC/seqlen * 100., 3)
            d['perc_CpG'] = round(total_CpG/seqlen * 100., 3)
            return d
    def _repetitive_content(self, d):
        try:
            ham = d['perc_Rep']
            return d
        except KeyError:
            seq = d['seq']
            seqlen = float(len(seq))
            repetitive = seq.count('a') + seq.count('c') + seq.count('g') + seq.count('t')
            d['perc_Rep'] = round(repetitive/seqlen * 100., 3)
            return d
    def _CGI_overlap(self, d):
        try:
            ham = d['perc_CGIovl']
            eggs = d['num_CGIovl']
            return d
        except KeyError:
            overlaps = d['cgioverlap']
            del d['cgioverlap']
            ovl_num = round((len(overlaps) / float(d['length'])) * 100., 3)
            ovl_perc = round((sum(map(int, overlaps)) / float(d['length'])) * 100., 3)
            d['num_CGIovl'] = ovl_num
            d['perc_CGIovl'] = ovl_perc
            return d
    def _peak_overlap(self, d):
        """Since this is a binary feature, there is nothing to do
        """
        return d
    def _CRM_motifs(self, d):
        try:
            # just check existence of one CRM type
            ham = d[self.crm[0]]
            return d
        except KeyError:
            motifs = d['crmmotifs']
            del d['crmmotifs']
            regflen = float(d['length'])
            countdict = dict(zip(self.crm, [0] * len(self.crm)))
            for m in motifs:
                countdict[m] += 1
            for k,v in countdict.items():
                d[k] = round(v / regflen * 100., 3)
            return d
    def _TFBS_motifs(self, d):
        try:
            ham = d[self.tfbs[0]]
            return d
        except KeyError:
            motifs = d['tfbsmotifs']
            del d['tfbsmotifs']
            regflen = float(d['length'])
            countdict = dict(zip(self.tfbs, [0] * len(self.tfbs)))
            pos0 = set(map(lambda s: s.strip("_cont"), motifs[0].split()))
            for m in pos0:
                countdict[m] += 1
            for pos in motifs[1:]:
                mlst = [ m for m in pos.split() if "_cont" not in m ] # do not double count continued motifs
                for m in mlst:
                    countdict[m] += 1
            for k,v in countdict.items():
                d[k] = round(v / regflen * 100., 3)
            return d
    def _signal_overlap(self, d):
        # since we don't know a priori which signal tracks were overlapped
        # we cannot check for already calculated overlaps
        all_sigfiles = [ k for k in d.keys() if "sigovl" in k ]
        regflen = float(d['length'])
        if all_sigfiles:
            for sigovl in all_sigfiles:
                keyname = sigovl.strip("_sigovl") # rest is just numeric file key
                values = d[sigovl]
                values = map(float, values)
                del d[sigovl]
                try:
                    perc = self.signal_perc[keyname]
                    keyname = keyname + "_"
                    p20, p40, p60, p80, p90, p95 = map(float, [perc['p20'], perc['p40'], perc['p60'], perc['p80'], perc['p90'], perc['p95']])
                    # there is probably a more pythonic way for the following...
                    d[(keyname + "p20")] = round((sum(1 for v in values if v <= p20)) / regflen * 100., 3)
                    d[(keyname + "p40")] = round((sum(1 for v in values if p20 < v <= p40)) / regflen * 100., 3)
                    d[(keyname + "p60")] = round((sum(1 for v in values if p60 < v <= p80)) / regflen * 100., 3)
                    d[(keyname + "p80")] = round((sum(1 for v in values if p80 < v <= p90)) / regflen * 100., 3)
                    d[(keyname + "p90")] = round((sum(1 for v in values if p90 < v <= p95)) / regflen * 100., 3)
                    d[(keyname + "p95")] = round((sum(1 for v in values if v > p95)) / regflen * 100., 3)
                except KeyError:
                    # Technically, the logging module is supposed to be threadsafe...
                    # However, this could result in thousands of error messages, so it is not
                    # the best solution...
                    self.logger.error("No percentiles available for signal track {}".format(sigovl))
                    continue
        return d
    def _kmer_frequency(self, k, d):
        kmerit = itr.product("ACGTN", repeat=k)
        kmerdict = {}
        while 1:
            try:
                kmerdict[ ("".join(next(kmerit))) ] = 0
            except StopIteration: break
        regflen = float(d['length'])
        seq = d['seq'].upper()
        wordfreqs = col.defaultdict(int)
        for i in range(0,k):
            words = re.findall('.{%i}' %k, seq[i:])
            for w in words:
                wordfreqs[w] += 1
        for k,v in wordfreqs.items():
            v = round(v / regflen * 100., 3)
            kmerdict[k] = v
        d.update(kmerdict)
        return d
            
        


"""
import multiprocessing as mp
import sys as sys
import logging as log
import functools as fnt

class CompFeat(object):
    def __init__(self, chunks):
        self.chunks = mp.Queue()
        self.logger = log.getLogger(__name__)
        log.basicConfig(stream=sys.stdout, level=log.DEBUG)
        self.out_buffer = mp.Queue()
        for d in chunks:
            self.chunks.put(d)
        self.chunks.put("STOP")
        self.chunks.put("STOP")
        self.workers = [mp.Process(target=self.run, kwargs={"k":[1,2,3], "j":"foo"}) for i in range(2) ]
        for p in self.workers:
            p.start()
        for p in self.workers:
            p.join()
        return
    def run(self, k, j):
        self.logger.debug("Process logging")
        print "started"
        print k
        print j
        p = mp.current_process()
        print p.pid
        for d in iter(self.chunks.get, "STOP"):
            clos = fnt.partial(self.triple, 3)  
            self.out_buffer.put(clos(d))
            clos = fnt.partial(self.triple, 5)
            self.out_buffer.put(clos(d))
            self.out_buffer.put(self.square(d))
            self.out_buffer.put(self.plusone(d))
        print "done"
        return
    def triple(self, l, d):
        for k,v in d.iteritems():
            d[k] = v * l
        return d
    def square(self, d):
        print d
        for k,v in d.iteritems():
            d[k] = v ** 2
        return d
    def plusone(self, d):
        for k,v in d.iteritems():
            d[k] = v + 1
        return d
    def get_result(self):
        self.out_buffer.put("STOP")
        outdict = {}
        for d in iter(self.out_buffer.get, "STOP"):
            outdict.update(d)
        return outdict

chunks = [ {'a':2, 'b':4}, {'c':3, 'd':5} ]
p = mp.current_process()
print p.pid
obj = CompFeat(chunks)
result = obj.get_result()
print result
print "done"
log.shutdown()
sys.exit(0)
"""
