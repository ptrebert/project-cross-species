""" seqmatchfinder module handles the computational burden
to identify matching [similar] sequences for a given set of
sequences.
Input is expected to be list of region dictionaries at least
holding the entries 'chrom', 'seq', 'perc_GC', 'perc_Rep' and 'length'
(as long as selection is based on these features)
"""

import os as os
import time as ti
import random as rand
import logging as log
import multiprocessing as mp
from queue import Empty as EmptyException
from queue import Full as FullException

from modules.datastructs import SizeTree
from modules.datastructs import OverlapTree

class TreeSeqMatchFinder:
    def __init__(self, config, featfunc):
        self.logger = log.getLogger(__name__)
        self._feature_functions = featfunc
        self.config = config
    def execute(self, queryregions, complregions):
        NUM_MATCH = len(queryregions)
        suglim = int(2*NUM_MATCH if self.config['suggestlimit'] == "auto" else self.config['suggestlimit'])
        mpman = mp.Manager()
        finished = mpman.Event()
        fragmentq = mpman.Queue()
        matchq = mpman.Queue()
        testq = mpman.Queue()
        matchedids = mpman.dict()
        curr_relax = float(self.config['initrelax'])
        querytree, len_low, len_high = self._build_querytree(queryregions, matchedids, curr_relax)
        overlaplst = mpman.list()
        overlaplst.append(OverlapTree())
        self.logger.debug("Start looping...")
        start_time = ti.time()
        while not finished.is_set():
            generator = mp.Process(target=self.generate_fragments,
                    args=(complregions, fragmentq, finished, curr_relax, suglim, len_low, len_high))
            generator.daemon = True
            generator.start()
            generator.join(0)
            matcher = mp.Process(target=self.match_query_fragment,
                            args=(querytree, fragmentq, matchq, finished))
            matcher.daemon = True
            matcher.start()
            matcher.join(0)
            ovltester = mp.Process(target=self.test_fragment_overlap,
                            args=(matchq, overlaplst, matchedids, finished, NUM_MATCH))
            ovltester.daemon = True
            self.logger.debug("Starting and joining overlap checker - waiting "
                                "for process to return...")
            ovltester.start()
            ovltester.join()
            self.logger.info("Returned - found {} matches so far".format(len(matchedids)))
            curr_relax += float(self.config['stepsize'])
            if curr_relax > float(self.config['relaxlimit']):
                curr_relax = float(self.config['initrelax'])
            querytree, len_low, len_high = self._build_querytree(queryregions, matchedids, curr_relax)
            time_now = ti.time()
            if (time_now - start_time) / 3600 > float(self.config['timelimit'])\
                and (len(matchedids) / NUM_MATCH) > (float(self.config['minmatch']) / 100.):
                finished.set()
                self.logger.info("Reached limit for runtime!")
            continue
        # finished was set, collect results
        matchlst = overlaplst[0].collect_matches()
        self.logger.debug("Collected {} matches from OverlapTree".format(len(matchlst)))
        mergedlst = self._merge_regions(queryregions, matchlst)
        return mergedlst
    def _merge_regions(self, queryregions, matchlst):
        assembly = queryregions[0]['assembly']
        last_id = max([qr['id'] for qr in queryregions])
        posdict = { qr['id']:qr for qr in queryregions }
        reslst = []
        for m in matchlst:
            last_id += 1
            pos = posdict[m['query_id']]
            pos['relax'] = -1
            pos['matching_id'] = last_id
            del m['complement_id']
            del m['query_id']
            m['assembly'] = assembly
            m['id'] = last_id
            m['matching_id'] = pos['id']
            reslst.extend([pos, m])
            assert pos['id'] == m['matching_id'] and m['id'] == pos['matching_id'],\
                "Mixing up ids/matching ids"
            assert pos.keys() == m.keys(), ("Different sets of keys: \
                                        {} vs {}".format(pos.keys(), m.keys()))
        assert len(reslst) <= 2*len(queryregions), ("Merging resulted in "
                                        "several matches per query")
        return reslst
    def generate_fragments(self, complreglst, outqueue, finished, relax,
                             suglim, len_low, len_high):
        '''Generate fragment suggestions for matches out of the list of
        complementary fragments; do so until LIMIT is reached or
        FINISHED flag is set

        complreglst: list of complementary regions
        finished: finished flag, shared event
        outqueue: shared queue to put suggestions
        relax: current relaxation in perc. points
        suglim: number of suggestions to create before returning
        '''
        num_sug = 0
        posidxseq = list(range(len(complreglst)))
        add_features = self._feature_functions
        LIMIT = int(suglim)
        FILTERN = float(self.config['filtern'])
        while not finished.is_set() and num_sug <= LIMIT:
            rand.shuffle(posidxseq)
            for idx in posidxseq:
                complreg = complreglst[idx]
                complreg_len = int(complreg['end']) - int(complreg['start'])
                fraglen = rand.randint(len_low, len_high)
                try:
                    fragstart = rand.randint(0, complreg_len - fraglen)
                except ValueError:
                    continue
                fragseq = complreg['seq'][fragstart:fragstart + fraglen]
                undef = fragseq.count('N') + fragseq.count('n')
                if (undef / fraglen) * 100. >= FILTERN: continue
                fragstart = complreg['start'] + fragstart
                fragend = fragstart + fraglen
                assert fraglen == (fragend - fragstart) == len(fragseq), \
                    ("Possible OBO error: fragment l/e-s/seql {}/{}/{}"\
                    .format(fraglen, fragend-fragstart, len(fragseq)))
                suggest = {'chrom':complreg['chrom'], 'start':fragstart,
                            'end':fragend, 'length':fraglen, 'seq':fragseq,
                            'complement_id':complreg['id'], 'relax':relax,
                            'class':"NEG", 'label':-1}
                suggest = add_features(suggest)
                outqueue.put(suggest)
                num_sug += 1
        if finished.is_set(): return
        outqueue.put("LIMIT")
        return
    def match_query_fragment(self, querytree, suggestions, matches, finished):
        '''Check if a suggestion created by generator is indeed a match for
        any query region

        querytree: tree of query regions, shrinks over time
        suggestions: shared queue with suggestions created by generator
        matches: shared queue to put fragments that match a query
        finished: finished flag, shared event
        '''
        while not finished.is_set():
            suggest = suggestions.get()
            if suggest == "LIMIT": break
            match = querytree.find_match(suggest) # can be NONE
            if not match: continue
            suggest['query_id'] = match['id']
            matches.put(suggest)
        if finished.is_set(): return
        matches.put("LIMIT")
        return
    def test_fragment_overlap(self, matches, ovltreelst, matchedids, finished, num_match):
        '''Build a tree with the matching fragments to ensure they do not
        overlap; discard fragment if it overlaps with previous one

        matches: shared queue where matching fragments are put
        chromtreelst: shared list with OverlapTree
        matchedids: shared dict to keep track of what query has been matched by
                    which fragment
        finished: finished flag, shared event
        num_match: number of queryregions
        '''
        ovltree = ovltreelst[0]
        while not finished.is_set():
            match = matches.get()
            if match == "LIMIT": break
            if match['query_id'] in matchedids: continue
            ovl_chk = ovltree.is_overlapping(match)
            if ovl_chk: continue
            matchedids[match['query_id']] = match['complement_id']
            if len(matchedids) >= num_match: # make sure it stops
                finished.set()
                assert len(matchedids) == num_match, ("Found more matches than \
                    queries: {} vs {}").format(len(matchedids), num_match)
        ovltreelst[0] = ovltree
        return
    def _build_querytree(self, queryregions, matchedids, relaxation):
        '''Build the search tree for the query regions
        
        queryregions: list with query regions
        matchedids: shared list with already matched ids [matched queries]
        relaxation: relaxation in percentage points
        '''
        if len(queryregions) == len(matchedids): return None, 0, 0
        to_be_matched = [qr for qr in queryregions\
                            if qr['id'] not in matchedids.keys()]
        item = to_be_matched[0]
        llst = [item['length']]
        stree = SizeTree(item, relaxation)
        for qr in to_be_matched[1:]:
            llst.append(qr['length'])
            stree.insert_region(qr)
        return stree, min(llst), max(llst)
    

class IndexSeqMatchFinder:
	def __init__(self, initrelax, stepsize, maxrelax, filterN, featurefunc):
		self.logger = log.getLogger(__name__)
		self._relax_values = list(self._frange(initrelax, maxrelax, stepsize))
		self._filterN = float(filterN)
		self._featurefunc = featurefunc
	def _frange(self, init, maxval, stepsize):
		while init <= maxval:
			yield init
			init += stepsize
	def execute(self, genomic_regions, complements):
		complements = self._enum_complements(complements)
		man = mp.Manager()
		finished = man.Event()
		unfinished = man.Event()
		matched_regions = man.dict()
		relax_gen = map(lambda x: x, self._relax_values)
		while not finished.is_set():
			try:
				current_relax = next(relax_gen)
			except StopIteration:
				relax_gen = map(lambda x: x, self._relax_values) # start over again
				current_relax = next(relax_gen)
			self.logger.debug("Building search index")
			search_index, minlen, maxlen = self._build_search_index(genomic_regions, current_relax)
			self.logger.debug("Search index built with {} entries".format(search_index.num_regions))
			suggested = man.Queue()
			suggest = mp.Process(target=self._generate_suggestions, args=(complements, minlen, maxlen, suggested, 2*len(genomic_regions), finished) , daemon=True)
			self.logger.debug("Starting suggestion generator")
			suggest.start()
			processed = man.Queue()
			sugprocess = mp.Process(target=self._process_suggestions, args=(suggested, processed, self._filterN, self._featurefunc) , daemon=True)
			self.logger.debug("Starting suggestion processor")
			sugprocess.start()
			matchfinder = mp.Process(target=self._search_matches, args=(search_index, processed, matched_regions, len(genomic_regions), finished, unfinished) , daemon=True)
			self.logger.debug("Starting match finder")
			matchfinder.start()
			# can two event.wait() calls be combined???
			while (not finished.is_set()) or (not unfinished.is_set()):
				ti.sleep(10)
			self.logger.debug("Joining all processes")
			suggest.join(1)
			sugprocess.join(1)
			matchfinder.join(1)
			[p.terminate() for p in [suggest, sugprocess, matchfinder] if p.is_alive()]
			self.logger.debug("Terminated all processes (if necessary)")
			self.logger.info("Found {} matches in this iteration".format(len(matched_regions)))
			# DEBUG
			self.logger.debug("Exiting here...")
			break
		return False
	def _enum_complements(self, complements):
		count = 1
		res = []
		for d in complements:
			d["complement_id"] = count
			res.append(d)
			count += 1
		return res
	def _build_search_index(self, genomic_regions, relax):
		genomic_regions.sort(key=lambda x: x['length'])
		median_element = genomic_regions[len(genomic_regions) // 2]
		min_length = genomic_regions[0]['length']
		max_length = genomic_regions[-1]['length']
		assert min_length <= max_length, "Region lengths are inconsistent - min is {} and max is {}".format(min_length, max_length)
		genomic_regions.remove(median_element)
		search_tree = SizeTree(median_element, relax)
		insert = search_tree.insert_region
		for gr in genomic_regions:
			insert(gr)
		return search_tree, min_length, max_length
	def _generate_suggestions(self, complements, minlen, maxlen, outqueue, num_suggest, finished):
		sugg_submit = 0
		while not finished.is_set() or sugg_submit < num_suggest:
			rand.shuffle(complements)
			for c in complements:
				frag_len = rand.randint(minlen, maxlen)
				if c['length'] - (frag_len + 1) <= frag_len: continue
				start_pos = rand.randint(c['start'], c['end'] - (frag_len + 1))
				suggestion = { "chrom":c["chrom"], "start":start_pos, "end":start_pos + frag_len, "seq":c["seq"][start_pos - c["start"]:start_pos - c["start"] + frag_len + 1], "complement_id":c["complement_id"] }
				outqueue.put(suggestion)
				sugg_submit += 1
		outqueue.put("Done")
		return
	def _process_suggestions(self, inqueue, outqueue, filterN, featurefunc):
		for s in iter(inqueue.get, "Done"):
			if (s["seq"].lower().count("n") / s["length"]) * 100. > filterN: continue
			s = featurefunc(s)
			outqueue.put(s)
		outqueue.put("Done")
		return
	def _search_matches(self, search_tree, inqueue, outdict, max_match, finished, unfinished):
		found_matches = 0
		for s in iter(inqueue.get, "Done"):
			match = search_tree.find_match(s)
			if match:
				outdict[match['id']] = s
				found_matches += 1
		if found_matches >= max_match:
			finished.set()
		else:
			unfinished.set()
		return

class SeqMatchFinder(object):
	def __init__(self, genome_regions, complements, initrelax, stepsize, relaxlimit, nindex, cyclelimit, min_limit_matches):
		self.logger = log.getLogger(__name__)
		self.regions = mp.Queue()
		self.complements = mp.Queue()
		self.matches = mp.Queue()
		# count number of matches backwards
		self.match_counter = mp.Value("i", len(genome_regions))
		self.last_match = mp.Value("f", 0.)
		# cycle until at least this many matches are found
		self.min_limit_matches = int(len(genome_regions) * (1 - (float(min_limit_matches) / 100.)))
		self.logger.info("Processes are allowed to abort computation if there are less than {} regions unmatched".format(self.min_limit_matches))
		assert self.min_limit_matches <= len(genome_regions), "Cannot find more matches than there are genomic regions"
		# how many times was the relaxation reset?
		self.reset_counter = mp.Value("i", 0)
		self.relax = mp.Value("f", float(initrelax))
		self.checkout = mp.Value("i", 0)
		self.finished = mp.Event()
		self.stepsize = float(stepsize)
		self.relaxlimit = float(relaxlimit)
		self.nindex = nindex
		# how many times can the relaxation be reset before accepting mediocre solutions?
		self.cyclelimit = cyclelimit
		self.lock = mp.Lock()
		for gr in genome_regions:
			self.regions.put(gr)
		self.regions.put("SENTINEL")
		# decrease possible "bias" towards complementary regions
		# at the beginning of the queue
		rand.shuffle(complements)
		for c in complements:
			self.complements.put(c)
		self.logger.info("Counter starts at: {}".format(self.match_counter.value))
	def execute(self, ncpu, assembly, last_id, regclass, label):
		"""Control over the child processes
		Ensure proper termination of processes
		in case of exceptions
		ncpu: number of free CPUs / processes to create
		assembly: to which genomic assembly are the coordinates compliant
		last_id: last_id + 1 gives first id for sequence matches (e.g. POS 1-34721 NEG 34722 - ...)
		regclass: class of the identified regions (usually "NEG")
		label: class label for these regions (usually -1)
		"""
		processes = []
		try:
			self.checkout.value = ncpu
			for n in range(ncpu):
				p = mp.Process(target=self._find_match, name=str(n))
				p.daemon = True
				processes.append(p)
			self.logger.debug("Created {} child processes, starting computation now".format(ncpu))
			for p in processes:
				p.start()
			for p in processes:
				p.join(0)
			# collect results
			results = []
			current_id = last_id + 1
			while True:
				try:
					while True:
						item = self.matches.get_nowait()
						item = self._meta_annotate(item, current_id, assembly, regclass, label)
						results.append(item)
						current_id += 1
				except EmptyException:
					if self.finished.is_set():
						ti.sleep(1)
						if self.matches.qsize() <= 0 and self.match_counter.value <= 0: break
						elif self.checkout.value <= 0:
							self.logger.info("All workers checked out, assume enough matches have been found")
							break
						continue
					else:
						ti.sleep(5)
						continue
			self.logger.info("Finished flag was set, checking for survivors...")
			for p in processes:
				if p.is_alive():
					self.logger.debug("Process still alive: {}".format(p.name))
					p.join(1)
					p.terminate()
			return results
		except Exception as e:
			self.logger.error("Caught exception {}, terminating all child processes".format(e))
			for p in processes:
				if p.is_alive():
					p.join(1)
					p.terminate()
			self.logger.warning("All processes terminated - possible data corruption, do not use the data of this run")
			return []
	def _meta_annotate(self, item, current_id, assembly, regclass, label):
		item['assembly'] = assembly
		item['class'] = regclass
		item['label'] = label
		item['id'] = current_id
		return item
	def _find_match(self):
		"""Iterate over both queues (regions and complements)
		If SENTINEL in regions is found, increase relaxation parameter
		"""
		initialrelax = self.relax.value
		myname = mp.current_process().pid
		last_sentinel_catch = 0.0
		sentinel_catch_counter = 0
		while True:
			try:
				if self.match_counter.value <= 0:
					with self.lock:
						self.logger.debug("{} exits - match counter reached zero".format(myname))
						self.checkout.value -= 1
						if self.checkout.value <= 0:
							self.logger.debug("{} is setting the finished flag".format(myname))
							self.finished.set()
					break
				if sentinel_catch_counter > 100 and self.checkout.value > self.match_counter.value:
					# there are more processes active than regions that need to be processed
					with self.lock:
						self.logger.debug("{} exits, too many processes around - {} regions still need to be matched".format(myname, self.match_counter.value))
						self.checkout.value -= 1
						assert self.checkout.value > 0, "Last process terminated, but work is not done"
					break
				# need a second break condition to avoid cycling for too long if some regions are very hard
				# to be matched - so if no match was found within the last 15 minutes, leave the stage...
				elif self.match_counter.value < self.min_limit_matches and (ti.time() - self.last_match.value) > 900.:
					# found enough matches
					with self.lock:
						self.logger.debug("{} exits - found already minimum number of required matches".format(myname))
						self.checkout.value -= 1
						if self.checkout.value <= 0:
							self.logger.debug("%s is setting the finished flag" %myname)
							self.finished.set()
					break						
				region = None
				try:
					region = self.regions.get_nowait()
				except EmptyException:
					continue
				if region == "SENTINEL":
					if (ti.time()) - last_sentinel_catch > 2.: # was "long" ago, assume not the same process
						with self.lock:
							self.relax.value += self.stepsize
							if self.relax.value > self.relaxlimit:
								self.logger.debug("{} is resetting relaxation parameter".format(myname))
								self.relax.value = initialrelax
								self.reset_counter.value += 1
					while True:
						try:
							self.regions.put("SENTINEL", True, 1)
							break
						except FullException:
							continue
					last_sentinel_catch = ti.time()
					sentinel_catch_counter += 1
					continue # sentinel set, continue with outer loop
				# have a region, need a complement
				complement = None
				while True:
					try:
						complement = self.complements.get_nowait()
						break
					except EmptyException:
						continue
				match = self._check_match(region, complement)
				if match:
					match, remains = self._extract_match(match, complement)
					with self.lock:
						while True:
							try:
								self.matches.put(match, True, 1)
								break
							except FullException:
								continue
						while True:
							for i in range(len(remains)): # this should be at most two
								try:
									self.complements.put(remains[i], True, 1)
								except FullException:
									remains = remains[i:]
									continue
							break
						self.match_counter.value -= 1
						self.logger.debug("{} found match - counter at: {}".format(myname,self.match_counter.value))
						self.last_match.value = ti.time()
						continue # with outer loop
				while True:
					try:
						self.regions.put(region, True, 1)
						break
					except FullException:
						continue
				while True:
					try:
						self.complements.put(complement, True, 1)
						break
					except FullException:
						continue
				continue
			except Exception as e: # outer try block 
				with self.lock:
					self.logger.error("{} caught exception {}".format(myname, e))
				raise e
		return None				
	def _check_match(self, region, complement):
		"""Generate list of random start indices and check if selected
		sequence is a possible match
		"""
		compl_len = complement['end'] - complement['start']
		region_len = region['length']
		# allow for 10% variation wrt region length
		if (region_len * 1.1) >= compl_len: return None
		remaining_length = int(compl_len - region_len * 1.1)
		indices = rand.sample(xrange(0, remaining_length, 1), min(self.nindex, remaining_length)) # in case complementary is only a bit longer than region
		compl_seq = (complement['seq']).lower()
		region_gc = region['perc_rawGC']
		limit_low = int(region_len * 0.9)
		limit_high = int(region_len * 1.1)
		region_rep = region['perc_Rep']
		# self.relax gives percent of allowed variation as percentage points
		relax = self.relax.value
		reset_counter = self.reset_counter.value
		for i in indices:
			match_len = rand.randint(limit_low, limit_high)
			test_seq = compl_seq[i:i + match_len]
			test_gc = round(((test_seq.count('c') + test_seq.count('g')) / float(match_len)) * 100., 3)
			test_rep = round(((test_seq.count('a') + test_seq.count('c') + test_seq.count('g') + test_seq.count('t')) / float(match_len)) * 100., 3)
			if region_gc - relax <= test_gc <= region_gc + relax:
				# GC content is considered more important, so if no matches are found for too long
				# simply return a half-way appropriate match
				if reset_counter > self.cyclelimit:
					# turn the id into a negative one to indicate "weak" matches
					return (i, i + match_len, (region['id']) * -1, reset_counter, relax)
				if test_rep >= region_rep - relax and test_rep <= region_rep + relax:
					return (i, i + match_len, region['id'], reset_counter, relax)
			continue
		return None
	def _extract_match(self, match, complement):
		"""Main task is to convert between string coordinates
		and the actual genomic coordinates
		"""
		str_start = match[0]
		str_end = match[1]
		matching_id = match[2]
		reset_counter = match[3]
		relax = match[4]
		match = {}
		match['seq'] = (complement['seq'])[str_start : str_end]
		match['start'] = complement['start'] + str_start # this needs to be checked for 0-1 based problems
		match['end'] = complement['start'] + str_end # this needs to be checked for 0-1 based problems
		match['chrom'] = complement['chrom']
		match['length'] = match['end'] - match['start']
		match['matching_id'] = matching_id
		match['reset_counter'] = reset_counter
		match['relax'] = relax
		remainder = []
		c1 = {}
		c1['start'] = complement['start']
		c1['end'] = match['start']
		c1['length'] = c1['end'] - c1['start']
		if c1['length'] >= match['length'] * 0.8: # avoid that many very small fragments clog the queue
			c1['seq'] = (complement['seq'])[ : c1['end']]
			c1['chrom'] = complement['chrom']
			remainder.append(c1)
		c2 = {}
		c2['start'] = match['end']
		c2['end'] = complement['end']
		c2['length'] = c2['end'] - c2['start']
		if c2['length'] >= match['length'] * 0.8:
			c2['seq'] = (complement['seq'])[c2['start'] : ]
			c2['chrom'] = complement['chrom']
			remainder.append(c2)
		return match, remainder	
