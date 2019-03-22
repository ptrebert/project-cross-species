"""Prediction module
Binary classification problem: presence/absence of peak in genomic region
Two use cases: prepare_dataset and predict_peaks
Expects as input set of genomic regions (positives) as BED format file
Computes matching set of negative regions (roughly same feature composition)
Or
Expects as input complete dataset with pos/neg samples and sends this to R
for machine learning part of pipeline
"""

import os as os
import logging as log
import collections as col
import io as io
import subprocess as sp
import tempfile as tf
import datetime as dat
import numpy as np
import gc as gc
import pickle as pckl
import gzip as gz
import resource as rsrc

import modules.dbmanager as dbm
import modules.toolbox as tlbx
import modules.regexpmanager as regexpman
import modules.osinterface as osi
import modules.bedfeatures as bedf
import modules.seqmatchfinder as seqmf
import modules.debugging as dbg

class Prediction(object):
    def __init__(self, configuration):
        self.config = configuration
        self.logger = log.getLogger(__name__)
    # ==============================================
    # first use case: prepare datasets: read peak file, find negatives
    def prepare_datasets(self, dbmanager, toolbox):
        """prepare_datasets computes a complete dataset (finding negatives)
        for given input datafiles in BED format
        """
        status = True
        self.logger.debug("prepare_datasets started")
        peakfiles_dict = dbmanager.create_DBdict("genomicregions")
        datasets_dict = dbmanager.create_DBdict("datasets")
        peakfiles_lst = self._check_existing(peakfiles_dict, datasets_dict)
        if not peakfiles_lst:
            # all peak files have been processed already - or error
            self.logger.info("Found a complete dataset for each peak file or none marked as usable - returning...")
            return status
        # we are only interested in regions on specific chromosomes
        match_chrom = regexpman.RegExpManager().get_chrom_matcher(self.config['predchrom'])
        for pkfd in peakfiles_lst:
            pkreg_lst, region_counter = self._parse_datafile(pkfd, match_chrom)
            pkreg_lst = self._annotate_sequence(pkreg_lst, pkfd['assembly'], dbmanager, toolbox)
            pkreg_lst = self._compute_select_features(pkreg_lst, self.config['selectfeat'], self.config['kmers'])
            complreg_lst = self._build_sequence_complement(pkreg_lst, dbmanager, toolbox, match_chrom)
            complreg_lst = self._annotate_sequence(complreg_lst, pkfd['assembly'], dbmanager, toolbox)
            merged_lst = self._find_similar_regions(pkreg_lst, complreg_lst)
            # Merge positives and negatives and prepare dumping to file on disk
            if (self.config['compfeat']).lower() in ['y', '1', 'yes', 'oh yeah', 'please', 'sure', 'of course', 'positive']:
                merged_dict = { d['id']:d for d in merged_lst }
                tuplist = self._prepare_feature_computation(merged_dict, pkfd, self.config['features'], dbmanager, toolbox)
                newdbentry_noid = self._create_prelim_metadata(pkfd, self.config['features'], "genomicregions", len(merged_lst))
                status = self._perform_feature_computation(tuplist, self.config['features'], self.config['kmers'], dbmanager, newdbentry_noid)
            else: # write a genomic regions file for later computation of features; is yet a valid entry of completedatasets
                status = self._write_clean_completedataset(merged_lst, pkfd, "genomicregions", dbmanager)  
            toolbox.clean_up()
            self.logger.info("Preparation of dataset for peakfile {} completed".format(pkfd['filename']))
        return status
    def _check_existing(self, peakfiles_dict, datasets_dict):
        """Check if for a given genomic regions files, there is already a complete dataset
        Avoid looking again for matching set of negative regions
        """
        # TODO: this function is no longer doing the right thing
        #existing = set([ int(d['link']) for d in datasets_dict.values() ])
        #peakfile_lst = [ pkfd for pkfd in peakfiles_dict.values() if pkfd['status'] == "use" and pkfd['key'] not in existing ]
        #self.logger.debug("Found {} region files to process".format(len(peakfile_lst)))
        peakfile_lst = [ pkfd for pkfd in peakfiles_dict.values() if pkfd['status'] == "use"]
        return peakfile_lst
    def _parse_datafile(self, filedict, matcher):
        """Parse a BED-like file, assign an ID to each row
        and returns a list of dictionaries describing each genomic region
        """
        all_regions = []
        region_counter = 0
        with open(os.path.join(filedict['dir'], filedict['filename']), "r") as infile:
            assembly = filedict['assembly']
            for line in infile:
                try:
                    cols = line.strip().split("\t")
                    if not matcher(cols[0]): continue # this probably takes care of any comment already...
                    region_counter += 1
                    region_dict = {}
                    region_dict['chrom'] = cols[0]
                    region_dict['start'] = int(cols[1])
                    region_dict['end'] = int(cols[2])
                    region_dict['length'] = region_dict['end'] - region_dict['start']
                    region_dict['id'] = region_counter
                    region_dict['assembly'] = assembly
                    region_dict['class'] = "POS"
                    region_dict['label'] = 1
                    all_regions.append(region_dict)
                except ValueError:
                    # in case there are comments/weird entries in the file
                    continue
        self.logger.debug("Parsed file {} - {} regions matching criterion".format(filedict['filename'], region_counter))
        return all_regions, region_counter
    def _annotate_sequence(self, region_lst, assembly, dbmanager, toolbox):
        """Receives a list of region dictionaries and annotates these with the 
        corresponding sequence (adding 'seq' entry to each dictionary)
        """
        self.logger.debug('Annotating regions with sequence information')
        dbgenomes = dbmanager.create_DBdict("genomes")
        genome_loc = (dbgenomes[assembly])['dir']
        annotated_regions = toolbox.ant_sequence(region_lst, genome_loc, float(self.config['filtern']))
        return annotated_regions
    def _compute_select_features(self, region_lst, features, kmers):
        """Create a BEDFeatures object and execute feature calculation
        Assumes that the dictionaries in region_lst have at least the entries
        'chrom' 'start' 'end' 'seq'
        """
        bedfeat = bedf.BEDFeatures(features, None, kmers)
        cpu_est = osi.estimate_free_cpus()
        self.logger.debug("Control handed over to BEDFeatures object")
        region_lst = bedfeat.execute(region_lst, min(cpu_est, len(region_lst)))
        return region_lst
    def _build_sequence_complement(self, region_lst, dbmanager, toolbox, chrom_matcher):
        """Builds the sequence complement for a given set of BED regions
        """
        self.logger.debug("Building sequence complement...")
        bedout = "\n".join(map(lambda reg: "\t".join([reg['chrom'], str(reg['start']), str(reg['end'])]), region_lst))
        bedfile = toolbox.wrt_str_tempfile(bedout)
        assembly = region_lst[0]['assembly']
        dbchromsize = dbmanager.create_DBdict("chromsizes")
        chromsz_lst = toolbox.prs_chrom_sizes(dbchromsize, assembly, chrom_matcher)
        # chrom sizes: tuples (chrom, size)
        chromout = "\n".join(map(lambda t: "\t".join([t[0], str(t[1])]), chromsz_lst))
        chszfile = toolbox.wrt_str_tempfile(chromout)
        complreg_file = toolbox.scl_bed_complement(bedfile, chszfile, "TEMPFILE")
        complregs = toolbox.prs_genome_region_file(complreg_file.name, 3, None, None, "dict")
        tmp = []
        for i in range(1, len(complregs) + 1):
            td = complregs[i-1]
            td['id'] = i
            tmp.append(td)
        complregs = tmp 
        self.logger.debug("Generated complementary regions for given BED file")
        return complregs
    def _find_similar_regions(self, pkreg_lst, complreg_lst):
        """Given a list of positives, find a matching set of negative regions
        Delegate process to SeqMatchFinder
        """
        self.logger.debug("Initializing TreeSeqMatchFinder...")
        compfeat = bedf.BEDFeatures(self.config['selectfeat'], None,
                                    self.config['kmers']).get_online_version()
        tsmf = seqmf.TreeSeqMatchFinder(self.config, compfeat)
        self.logger.debug("Executing")
        reslst = tsmf.execute(pkreg_lst, complreg_lst)
        self.logger.debug("Similar regions found, returning from normal operation")
        return reslst
    def _write_clean_completedataset(self, regions_lst, orig_dbentry, link_table, dbmanager):
        """Write dataset without features on disk
        Makes use of already existing methods
        """
        self.logger.debug("Writing a complete dataset to disk...")
        newdbentry_noid = self._create_prelim_metadata(orig_dbentry, "notcomputed", "genomicregions", len(regions_lst))
        compldat_id = str(dbmanager.get_next_id("completedatasets"))
        new_dbentry = newdbentry_noid
        new_dbentry['key'] = compldat_id
        new_dbentry['status'] = "ignore"
        new_dbentry['filename'] = new_dbentry['filename'] + compldat_id + ".bed"
        self.logger.debug("Updating DB with new entry")
        dbmanager.update_DBfile_no_id("completedatasets", [new_dbentry])
        # write complete dataset to file and continue
        good = self._write_complete_dataset(new_dbentry, regions_lst)
        self.logger.info("Wrote complete dataset {} to disk".format(os.path.join(new_dbentry['dir'], new_dbentry['filename'])))
        return good
    # ==============================================
    # use case: prepare expression data
    # create 2 complete datasets for each file [high/high_strict], assign ID to each entry
    def prepare_expression_data(self, dbmanager, toolbox):
        status = True
        self.logger.debug('prepare_expression_data started')
        dbexpr = dbmanager.create_DBdict('expression')
        processable_files = [ dbentry for dbentry in dbexpr.values() if dbentry['status'] == 'use' ]
        if not processable_files:
            self.logger.debug('No expression datasets marked as active, returning')
            return status
        for dbentry in processable_files:
            exprdata_lst, size = toolbox.prs_bed_wheader_wid(os.path.join(dbentry['dir'], dbentry['filename']))
            self.logger.debug('Read file {} with {} entries'.format(dbentry['filename'], len(exprdata_lst)))
            exprdata_lst = self._annotate_sequence(exprdata_lst, dbentry['assembly'], dbmanager, toolbox)
            self.logger.debug('Annotation with sequence finished')
            for exprval in ['high_expr', 'high_expr_strict']:
                fnid = 'normal' if exprval == 'high_expr' else 'strict'
                exprdata_lst = self._create_expression_label(exprdata_lst, float(dbentry[exprval]), dbentry['assembly'])
                compldat_id = str(dbmanager.get_next_id('completedatasets'))
                new_dbentry = self._create_expression_metadata(dbentry, 'notcomputed', 'expression', size)
                new_dbentry['key'] = compldat_id
                new_dbentry['filename'] = new_dbentry['filename'] + fnid + ".bed"
                self.logger.debug('Updating DB entry')
                dbmanager.update_DBfile_no_id("completedatasets", [new_dbentry])
                self.logger.debug("DB updated...")
                # write complete dataset to file and continue
                self.logger.debug('Writing complete dataset')
                status = status and self._write_complete_dataset(new_dbentry, exprdata_lst)
        toolbox.clean_up()
        self.logger.info('All expression datasets prepared')
        return status
    def _create_expression_label(self, exprdata_lst, threshold, assembly):
        """If RPKM of entry is greater than threshold -> high/1
        else -> low/-1
        """
        for entry in exprdata_lst:
            entry['assembly'] = assembly
            if float(entry['RPKM']) > threshold:
                entry['class'] = 'HIGH'
                entry['label'] = str(1)
            else:
                entry['class'] = 'LOW'
                entry['label'] = str(-1)
        assert 'class' in exprdata_lst[0], 'Setting labels for expression data failed: {}'.format(exprdata_lst[0])
        return exprdata_lst
    def _create_expression_metadata(self, dbentry, features, table, size):
        """Need key timestamp   dir filename    link    table   features    status  assembly    type    cell    size
        Key is not yet determined
        """
        metadict = {}
        metadict['key'] = "-1"
        metadict['timestamp'] = dat.date.today().isoformat()
        metadict['dir'] = self.config['outdir']
        metadict['filename'] = (dbentry['filename']).strip('.bed') + '_'
        metadict['link'] = dbentry['key']
        metadict['table'] = table
        metadict['size'] = str(size)
        metadict['assembly'] = dbentry['assembly']
        metadict['type'] = dbentry['element']
        metadict['cell'] = dbentry['cell']
        metadict['features'] = features
        metadict['status'] = "ignore"
        return metadict
    # ==============================================
    # use case: compute all features for a given set of complete
    # datasets (existence of pos/neg regions is implicitly assumed but not required
    def compute_features(self, dbmanager, toolbox):
        self.logger.debug("compute_features started")
        status = True
        dbcompdat = dbmanager.create_DBdict("completedatasets")
        processable_files = [ dbentry for dbentry in dbcompdat.values() if dbentry['status'] == 'use' and dbentry['features'] == 'notcomputed' ] # for adding more features to an already existing datset, this would need to be changed
        if not processable_files:
            self.logger.info("No complete datasets marked as active, returning")
            return status
        last_assembly = "godzilla"
        tfbs_chrom = None
        for dbentry in processable_files:
            region_dict = {}
            with open(os.path.join(dbentry['dir'], dbentry['filename']), "r") as infile:
                header = infile.readline().strip("\n#").split("\t") # strip also comment char; a header is implicitly assumed (danger?)
                for line in infile:
                    region = dict(zip(header, line.strip().split("\t")))
                    region_dict[region['id']] = region
            self.logger.debug("Read {} regions from file {}".format(len(region_dict), dbentry['filename']))
            region_dict = self._prepare_feature_computation(region_dict, dbentry, self.config['features'], dbmanager, toolbox)
            self.logger.debug("Preparation for feature computation done")
            newdbentry_noid = self._create_prelim_metadata(dbentry, self.config['features'], "completedatasets", len(region_dict))
            status = status and self._perform_feature_computation(region_dict, self.config["features"], self.config["kmers"], dbmanager, newdbentry_noid)
            self.logger.debug(dbg.print_mem_io_load())
        self.logger.debug("All features computed, returning from normal operation")
        toolbox.clean_up()
        return status
    def _prepare_tfbs_chromlists(self, assembly, motifdbentry, dbmanager, toolbox):
        """Get chromosome sizes for assembly and create lists of empty strings
        """
        self.logger.debug("Preparing lists with TF binding sites...")
        dbchromsizes = dbmanager.create_DBdict("chromsizes")
        chromsizes = toolbox.prs_chrom_sizes(dbchromsizes, assembly, regexpman.RegExpManager().get_chrom_matcher("autosomal"))
        resdict = {}
        for t in chromsizes:
            resdict[t[0]] = [""] * t[1]
        self.logger.debug("Preparing empty lists done")
        lines = []
        with open(os.path.join(motifdbentry['dir'], motifdbentry['filename']), "r") as infile:
            # emphasize speed over memory consumption here
            lines = infile.read().strip().split("\n")
        self.logger.debug("Reading content of TFBS file done - {} lines read".format(len(lines)))
        tfbs_chrom = self._split_tfbs_records(lines)
        assert tfbs_chrom, "Dictionary with TF binding sites sorted by chromosome is empty"
        del lines
        self.logger.debug("Generating lookup data structure")
        for chrom in resdict.keys():
            emptylist = resdict[chrom]
            del resdict[chrom]
            tfbslist = []
            try:
                tfbslist = tfbs_chrom[chrom]
                del tfbs_chrom[chrom]
            except KeyError:
                continue
            for t in tfbslist:
                emptylist[t[1]] += t[3] # Setting an item is O(1), setting a k-slice is O(k+n)
                motif = t[3] + "_cont"
                for i in range(t[1] + 1, t[2] + 1, 1):
                    emptylist[i] += motif
            resdict[chrom] = emptylist
        self.logger.debug("Annotation process done, freeing memory")
        del tfbs_chrom
        assert resdict, "Lists with TF binding sites are empty"
        self.logger.debug("Returning from normal operation")
        return resdict
    def _split_tfbs_records(self, bedlines):
        bed_chrom = col.defaultdict(list)
        for bedline in bedlines:
            cols = bedline.split("\t")
            bed_chrom[cols[0]].append( (cols[0], int(cols[1]), int(cols[2]), " " + cols[4]) ) # last part is fastest for simple concat
        return bed_chrom
    def _prepare_feature_computation(self, genreg_dict, dbentry, features, dbmanager, toolbox):
        """Prepare to add requested features to genomic regions
        It is assumed that the regions (dictionaries)
        have all necessary keys, i.e. as minimum chr, start,
        end, seq and unique identifier
        Checks are required for: signal/CRM/CGI overlap, TFBS motifs
        """
        self.logger.debug("Begin preparation for feature computation...")
        # dump the regions to tempfile for intersecting with data files
        outstr = [ "\t".join([d['chrom'], str(d['start']), str(d['end']), str(d['id'])]) for d in genreg_dict.values() ]
        outstr = "\n".join(outstr)
        regions_tmpf = toolbox.wrt_str_tempfile(outstr)
        self.logger.debug(dbg.print_mem_io_load())
        del outstr
        if "cgioverlap" in features:
            dbcpgislands = dbmanager.create_DBdict("cpgislands")
            dbentry_cgi = dbcpgislands[dbentry['assembly']]
            genreg_dict = self._add_overlap(genreg_dict, regions_tmpf, dbentry_cgi, toolbox, 3, 8, "cgioverlap", False)
        self.logger.debug(dbg.print_mem_io_load())
        if "crmmotifs" in features:
            dbcrms = dbmanager.create_DBdict("crm")
            dbentry_crm = dbcrms[dbentry['assembly']]
            genreg_dict = self._add_overlap(genreg_dict, regions_tmpf, dbentry_crm, toolbox, 3, 7, "crmmotifs", False)
        if "tfbsmotifs" in features:
            genreg_dict = self._add_tfbs_overlap(genreg_dict, dbmanager, toolbox)
        if "peakoverlap" in features:
            isect_files = self._get_peak_intersect(dbentry['key'], dbmanager, toolbox)
            for f in isect_files:
                feattype = f['key'] + "_peak"
                genreg_dict = self._add_overlap_count(genreg_dict, regions_tmpf, f, toolbox, 3, feattype,  False)
        self.logger.debug(dbg.print_mem_io_load())
        result = []
        if "signaloverlap" in features:
            tuplst = self._get_signal_intersect(dbentry['key'], dbmanager, toolbox)
            for t in tuplst:
                isect_files = t[0]
                percentiles = t[1]
                for f in isect_files:
                    feattype = f['key'] + "_sigovl"
                    genreg_dict = self._add_overlap(genreg_dict, regions_tmpf, f, toolbox, 3, 7, feattype, False)
                result.append( (genreg_dict, percentiles) )
        else:
            result.append( (genreg_dict, None) )
        assert result, "Result list with tuples (genregions, percentiles) is empty"
        self.logger.debug(dbg.print_mem_io_load())
        return result
    def _perform_feature_computation(self, tuplelist, features, kmers, dbmanager, newdbentry_noid):
        """Tuple list contains tuples of genomic regions and corresponding percentiles (or None)
        Send this to the BEDFeatures object, update DB and write complete datset to disk
        """
        for t in tuplelist:
            self.logger.debug(dbg.print_mem_io_load())
            self.logger.debug('Creating BEDFeatures object')
            bedfeat = bedf.BEDFeatures(features, t[1], kmers)
            cpu_est = osi.estimate_free_cpus()
            self.logger.debug("Control handed over to BEDFeatures object")
            region_lst = bedfeat.execute(t[0], min(cpu_est, len(t[0])))
            del bedfeat
            self.logger.debug("Done for 1 dataset")
            compldat_id = str(dbmanager.get_next_id("completedatasets"))
            new_dbentry = newdbentry_noid
            new_dbentry['key'] = compldat_id
            new_dbentry['filename'] = new_dbentry['filename'] + compldat_id + ".bed"
            dbmanager.update_DBfile_no_id("completedatasets", [new_dbentry])
            self.logger.debug("DB updated...")
            # write complete dataset to file and continue
            good = self._write_complete_dataset(new_dbentry, region_lst)
        return good
    def _add_tfbs_overlap(self, genreg_dict, dbmanager, toolbox):
        """Specialized function to add the overlap with the TFBS to all regions
        To reduce memory requirements, this functions parses on demand the file
        with the TFBS motifs per chrom (that are preprocessed and in binary format)
        """
        # sort regions according to chromosome
        self.logger.debug("Adding TFBS overlap feature")
        self.logger.debug(dbg.print_mem_io_load())
        sort_regions = sorted(genreg_dict.values(), key=lambda d: d['chrom'])
        lastchr = "CHR"
        tfbschrom = None
        for r in sort_regions:
            if r['chrom'] != lastchr:
                self.logger.debug("Next chromosome for assembly {}: {}".format(r['assembly'], r['chrom']))
                lastchr = r['chrom']
                tfbschrom = self._get_tfbsmotifs_chrom(r['assembly'], r['chrom'], dbmanager, toolbox)
                if not tfbschrom:
                    self.logger.debug("No TFBS motifs for {} on chromosome {}".format(r['assembly'], r['chrom']))
            if not tfbschrom: # for the case iff there no motifs for a chrom, though unlikely
                (genreg_dict[r['id']]).update({"tfbsmotifs":[""]})
                continue 
            (genreg_dict[r['id']]).update({"tfbsmotifs":tfbschrom[int(r['start']):int(r['end'])]})
        self.logger.debug(dbg.print_mem_io_load())
        return genreg_dict
    def _get_tfbsmotifs_chrom(self, assembly, chrom, dbmanager, toolbox):
        '''Check if a pre-processed file with TFBS motifs exists already in the
        database; create it if not
        '''
        tfbsmotifs_db = dbmanager.create_DBdict("tfbsmotifs")
        tfbschrom = None
        for entry in tfbsmotifs_db.values():
            if entry['assembly'] == assembly and entry['range_covered'] == chrom \
                and entry['format'] == "BIN":
                tfbschrom = self._parse_tfbsfile_binary(entry)
                self.logger.debug("Found TFBS motifs for {} on chromosome {}".format(assembly, chrom))
                return tfbschrom
        for entry in tfbsmotifs_db.values():
            if entry['assembly'] == assembly and entry['range_covered'] == "genome" and entry['format'] == "BED":
                    self.logger.debug("Need to prepare TFBS file for {} in assembly {}".format(chrom, assembly))
                    tfbschrom = self._create_tfbsfiles_binary(entry, chrom, dbmanager, toolbox)
            if entry['assembly'] == assembly and entry['range_covered'] == "autosomal" and entry['format'] == "BED":
                chrom_matcher = regexpman.RegExpManager().get_chrom_matcher("autosomal")
                if not chrom_matcher(chrom):
                    self.logger.debug("No raw TFBS data for {} in assembly {} available".format(chrom, assembly))
                    return tfbschrom
                self.logger.debug("Need to prepare TFBS file for {} in assembly {}".format(chrom, assembly))
                tfbschrom = self._create_tfbsfiles_binary(entry, chrom, dbmanager, toolbox)
        return tfbschrom # list of strings representing TFBS motifs
    def _create_tfbsfiles_binary(self, dbentry, chrom, dbmanager, toolbox):
        '''The TFBS motif file (BED format) from dbentry is split into chromosomes
        and the appropriate lists with TFBS motif sites are prepared and saved
        as gzips

        The BED file with TFBS motifs has coordinates as left- and right-inclusive
        1-based intervals
        '''
        self.logger.debug("Preparing binary files with TF binding sites...")
        dbchromsizes = dbmanager.create_DBdict("chromsizes")
        chromsizes = dict(toolbox.prs_chrom_sizes(dbchromsizes, dbentry['assembly'], regexpman.RegExpManager().get_chrom_matcher("autosomal")))
        keep_list = None
        # ensure that TFBS is chrom sorted (since needs only to be done once)
        existing_entries = self._check_existing_tfbs(dbentry, dbmanager)
        seen_chroms = set()
        last_chrom = ""
        tfbsmotifs_lst = []
        num_motifs = 0
        debug_count = 0
        with open(os.path.join(dbentry['dir'], dbentry['filename'])) as infile:
            self.logger.debug("Start reading of TFBS file...")
            for line in infile:
                debug_count += 1
                cols = line.strip().split("\t")
                if cols[0] in existing_entries or "#" in line: continue
                if cols[0] != last_chrom and cols[0] in seen_chroms: # unsorted
                    # TODO handle this by sorting the file and restart
                    raise TypeError("Expect TFBS file {} to be sorted by chromosome".format(dbentry))
                elif cols[0] != last_chrom:
                    if last_chrom == chrom: # list is chrom we want, keep it
                        keep_list = tfbsmotifs_lst
                    if tfbsmotifs_lst:
                        _ = self._write_tfbsfile_binary(last_chrom, dbentry, num_motifs, tfbsmotifs_lst, dbmanager)
                        self.logger.debug("Wrote file {} - {}".format(dbentry['assembly'], last_chrom))
                    last_chrom = cols[0]
                    seen_chroms.add(cols[0])
                    tfbsmotifs_lst = [""] * chromsizes[cols[0]] # make large empty list
                    num_motifs = 0
                # add motif to list
                # Since we want to annotate BED regions that are most likely
                # 0-based, convert to 0-based intervals
                # -> -1 for start, end can stay the same (non-inclusive in Python)
                start = int(cols[1])
                end = int(cols[2])
                num_motifs += 1
                motif = " " + cols[4]
                tfbsmotifs_lst[start - 1] += motif
                motif = motif + "_cont"
                for i in range(start, end, 1): # regular start since we set the first occurence manually
                    tfbsmotifs_lst[i] += motif
        if tfbsmotifs_lst:
            _ = self._write_tfbsfile_binary(last_chrom, dbentry, num_motifs, tfbsmotifs_lst, dbmanager)
            self.logger.debug("Wrote file {} - {}".format(dbentry['assembly'], last_chrom))
        if last_chrom == chrom:
            keep_list = tfbsmotifs_lst
        self.logger.debug("Splitting file {} into chromosomes is done - "
                            "read {} lines in total".format(dbentry['filename'], debug_count))
        return keep_list
    def _check_existing_tfbs(self, dbentry, dbmanager):
        '''Create a set of already existing TFBS motif entries
        '''
        existing = set()
        tfbs_db = dbmanager.create_DBdict("tfbsmotifs")
        for entry in tfbs_db.values():
            if entry['assembly'] == dbentry['assembly'] and entry['format'] == "BIN":
                existing.add(entry['range_covered'])
        return existing
    def _write_tfbsfile_binary(self, chrom, dbentry, num_motifs, tfbsmotifs_lst, dbmanager):
        '''Pickle the motif list and dump it to a gzipped file
        Then create appropriate database entry
        '''
        new_dbentry = dbentry
        new_dbentry['key'] = -1
        new_dbentry['format'] = "BIN"
        new_dbentry['range_covered'] = chrom
        new_filename = "_".join([new_dbentry['assembly'], new_dbentry['range_covered'], "TransfacMotifs", "ord0"]) + ".bin"
        new_dbentry['filename'] = new_filename
        new_dbentry['num_tfbs'] = str(num_motifs)
        outpath = os.path.join(new_dbentry['dir'], new_filename)
        with gz.open(outpath, "wb") as outfile:
            pckl.dump(tfbsmotifs_lst, outfile)
            outfile.flush()
        self.logger.debug("Wrote pickled object to file {}".format(outpath))
        dbmanager.update_DBfile_make_id("tfbsmotifs", None, [new_dbentry])
        self.logger.debug("Updated DB file")
        return
    def _parse_tfbsfile_binary(self, dbentry):
        inpath = os.path.join(dbentry['dir'], dbentry['filename'])
        tfbschrom_lst = None
        with gz.open(inpath, "rb") as infile:
            tfbschrom_lst = pckl.load(infile)
        assert tfbschrom_lst, "TFBS list from file {} is empty".format(inpath)
        return tfbschrom_lst
    def _add_overlap_count(self, genreg_dict, regions_tmpf, dbentry, toolbox, key, feattype, issorted):
        """Add info on overlap as mere counting information 0,1,2,3...
        Value is always in last column
        """
        self.logger.debug("Adding overlap count for feature {}".format(feattype))
        isect_outf = toolbox.scl_bed_intersect_count(regions_tmpf.name, os.path.join(dbentry['dir'], dbentry['filename']), "TEMPFILE", issorted)
        store_output = dict()
        with open(isect_outf.name, 'r') as infile:
            for line in infile:
                try:
                    cols = line.strip().split('\t')
                    store_output[cols[key]] = cols[-1]
                except IndexError as ie:
                    self.logger.debug('Intersect file {} is malformed\nFeature type {} and failed line {}'.format(isect_outf.name, feattype, line))
        del isect_outf
        for k in genreg_dict.keys():
            try:
                tmp = int(store_output[k])
                (genreg_dict[k]).update({feattype:"1" if tmp > 0 else "0"})
            except KeyError:
                (genreg_dict[k]).update({feattype:"0"})
        return genreg_dict
    def _add_overlap(self, genreg_dict, regions_tmpf, dbentry, toolbox, key, value, feattype, issorted):
        """Adds the info on overlaps with simple data files
        Intersect with CGI/CRM/signal/TFBS file, parse output and record required overlap info
        """
        self.logger.debug("Adding overlap for feature {}".format(feattype))
        isect_outf = toolbox.scl_bed_intersect_overlap(regions_tmpf.name, os.path.join(dbentry['dir'], dbentry['filename']), "TEMPFILE", issorted)
        store_output = col.defaultdict(list)
        with open(isect_outf.name, 'r') as infile:
            for line in infile:
                try:
                    cols = line.strip().split('\t')
                    store_output[cols[key]].append(cols[value])
                except IndexError as ie:
                    self.logger.debug('Intersect file {} is malformed\nFeature type {} and failed line {}'.format(isect_outf.name, feattype, line))
                    raise IndexError(str(ie))
        del isect_outf
        self.logger.debug(dbg.print_mem_io_load())
        lst = []
        for k in genreg_dict.keys():
            try:
                lst = store_output[k]
                (genreg_dict[k]).update({feattype:lst})
            except KeyError:
                (genreg_dict[k]).update({feattype:[]})
        return genreg_dict
    def _get_signal_intersect(self, fileid, dbmanager, toolbox):
        """Select signal tracks for intersection with
        genomic region file
        Compute percentiles if necessary
        """
        self.logger.debug("Collecting signal tracks for intersection...")
        dbsignaltracks = dbmanager.create_DBdict("signaltracks")
        dbdatasetlinks = dbmanager.create_DBdict("datasetlinks")
        entries = [ d for d in dbdatasetlinks.values() if d['link'] == fileid and d['status'] == "use" ]
        assert entries, "Received file id {} but this file is not in db_datasetlinks".format(fileid)
        result = []
        for e in entries:
            sigtracks = []
            perc = {}
            signallist = list(filter(lambda x: x.startswith('s'), e['intersectfiles'].split('-')))
            assert signallist, 'Did not find any signal files in entry: {}'.format(e)
            signal_files = signallist[0].strip('s()').split(',')
            for s in signal_files:
                try:
                    sigfile = dbsignaltracks[s]
                    if sigfile["p20"] == "-1": # percentiles have not been calculated
                        q = self._compute_signal_percentiles(sigfile, [20, 40, 60, 80, 90, 95])
                        sigfile["p20"] = q[0]
                        sigfile["p40"] = q[1]
                        sigfile["p60"] = q[2]
                        sigfile["p80"] = q[3]
                        sigfile["p90"] = q[4]
                        sigfile["p95"] = q[5]
                        dbmanager.update_DBentry(sigfile, "signaltracks")
                        dbsignaltracks[s] = sigfile
                    sigtracks.append(sigfile)
                    perc[s] = { "p20":float(sigfile["p20"]), "p40":float(sigfile["p40"]), "p60":float(sigfile["p60"]), "p80":float(sigfile["p80"]), "p90":float(sigfile["p90"]), "p95":float(sigfile["p95"]) }
                except KeyError:
                    self.logger.warning("Entry {} contains invalid signal track IDs, check DB".format(s))
                    break
            result.append( (sigtracks, perc) )
        return result
    def _get_peak_intersect(self, fileid, dbmanager, toolbox):
        self.logger.debug('Collecting peak track information')
        dbgenregion = dbmanager.create_DBdict("genomicregions")
        dbdatasetlinks = dbmanager.create_DBdict("datasetlinks")
        entries = [d for d in dbdatasetlinks.values() if d['link'] == fileid and d['status'] == "use"]
        assert entries, "Received file id {} but found no entry in datasetlinks".format(fileid)
        result = []
        for e in entries:
            region_files = list(filter(lambda x: x.startswith('g'), e['intersectfiles'].split('-')))
            assert region_files, 'Did not find any genomic region files in entry: {}'.format(e)
            peak_files = region_files[0].strip('g()').split(',')
            for p in peak_files:
                try:
                    peakfile = dbgenregion[p]
                    result.append(peakfile)
                except KeyError:
                    self.logger.warning("Entry {} contains invalid peak track IDs, check DB".format(e))
        return result
    def _compute_signal_percentiles(self, filedict, percentiles):
        infile = open(os.path.join(filedict['dir'], filedict['filename']), "r")
        self.logger.debug("Computing percentiles for file {}".format(infile.name))
        values = []
        for line in infile:
            try:
                # assumes 4 column signal file
                v = line[line.rfind("\t") :]
                values.append(float(v))
            except ValueError: continue # potential comments are ignored
        assert values, "Computing percentiles for file {} failed".format(infile.name)
        infile.close()
        perc = np.percentile(values, percentiles)
        perc = list(map(lambda x: str(round(x,3)), perc))
        return perc
    def _create_prelim_metadata(self, orig_dbentry, features, table, size):
        """Need timestamp, dir, filename, link, table, features, assembly, type, cell, size
        ID is not yet determined
        """
        metadict = {}
        metadict['key'] = "-1"
        metadict['timestamp'] = dat.date.today().isoformat()
        metadict['dir'] = self.config['outdir']
        metadict['filename'] = "_".join([orig_dbentry['assembly'], orig_dbentry['type'], orig_dbentry['cell'], ""])
        if 'gene' in orig_dbentry['type'] or 'TSS' in orig_dbentry['type']:
            metadict['filename'] = (orig_dbentry['filename']).strip('.bed') + '_'
        metadict['link'] = orig_dbentry['key']
        metadict['table'] = table
        metadict['size'] = str(size)
        metadict['assembly'] = orig_dbentry['assembly']
        metadict['type'] = orig_dbentry['type']
        metadict['cell'] = orig_dbentry['cell']
        metadict['features'] = features
        metadict['status'] = "ignore"
        return metadict
    def _order_keys_compldat(self, prelimorder):
        majors = ["chrom", "start", "end", "id", "length", "assembly", "class", "label", "matching_id", "relax"]
        minors = [ k for k in prelimorder if k not in majors and k != "seq"]
        minors = sorted(minors)
        majors.extend(minors)
        majors.append("seq")
        return majors
    def _order_keys_expression(self, prelimorder):
        majors = ["chrom", "start", "end", "strand", "id", "length", "RPKM", "ensembl_gene_id", "class", "label", "assembly"]
        minors = [ k for k in prelimorder if k not in majors and k != "seq"]
        minors = sorted(minors)
        majors.extend(minors)
        majors.append("seq")
        return majors
    def _write_complete_dataset(self, dbentry, genregions):
        prelimorder = genregions[0].keys()
        order = None
        if 'RPKM' in prelimorder:
            # expression data handled differently
            #TODO: generalize here
            order = self._order_keys_expression(prelimorder)
        else:
            order = self._order_keys_compldat(prelimorder)
        with open(os.path.join(dbentry['dir'], dbentry['filename']), "w") as outfile:
            # write a header
            outfile.write("#" + "\t".join(order) + "\n")
            allregions = [ ("\t".join([str(r[k]) for k in order])) for r in genregions]
            outline = "\n".join(allregions)
            outfile.write(outline)
            outfile.write("\n")
            outfile.flush()
            os.fsync(outfile.fileno())
        return True     
    # ==============================================
    # third use case: predict presence/absence of peaks for regions
    # or predict gene expression
    def predict_peaks(self, dbmanager, toolbox):
        status = True
        self.logger.debug("predict_peaks started")
        # basically, call to external R script
        return status
    def predict_gene_expression(self, dbmanager, toolbox):
        status = True
        self.logger.debug('predict_gene_expression started')
