"""Mapping module
Creates blank mappings
Post-processes blank mappings
Maps signal tracks from source to target assembly
Smoothing of mapped tracks is done using an external R script (lowess smoothing)
"""

import os as os
import logging as log
import collections as col
import io as io
import subprocess as sp
import tempfile as tf
import datetime as dat

import modules.dbmanager as dbm
import modules.regexpmanager as regexpman
import modules.osinterface as osi
from modules.datastructs import RegionCluster

class Mapping(object):
    def __init__(self, configuration):
        self.config = configuration
        self.logger = log.getLogger(__name__)
    # ==============================================
    # first use case: create a blank mapping
    def prepare_mapping(self, dbman, toolbox):
        status = True
        self.logger.debug("prepare_mapping started")
        # all these are lists of length >= 1
        source, target, resolution = self._split_cfgentries()
        chainfiles = dbman.create_DBdict("chainfiles")
        chromsizes = dbman.create_DBdict("chromsizes")
        reman = regexpman.RegExpManager()
        # match_chrom is now local method
        match_chrom = reman.get_chrom_matcher(str(self.config['mapchrom']))
        # check/create output directory for files to be mapped
        osi.check_or_create_dir(self.config['outdir'])
        instr_files = []
        for s in source:
            chromlist = toolbox.prs_chrom_sizes(chromsizes, s, match_chrom)
            if not chromlist:
                status = False
                self.logger.warning("prepare_mapping: list of accepted chromosomes for {} is empty - no DB entry?".format(s))
                continue
            input_files = []
            for r in resolution:
                splitfiles = self._generate_genome_tiling(s, r, chromlist)
                input_files.extend(splitfiles)
            assert input_files, "prepare_mapping: list of input files is empty"
            instr_files.append(self._generate_instruction_file(s, target, input_files, chainfiles))
        self.logger.info("prepare_mapping generated the following instruction files for liftOver mapping:\n" + "\n".join(instr_files))
        return status   
    def _split_cfgentries(self):
        """Split combined entries in configuration file
        into lists of length >= 1
        Returns lists for sources, targets, resolutions
        """
        source = self.config['source'].split("-")
        target = self.config['target'].split("-")
        resolution = self.config['resolution'].split("-")
        # might be a little too much, but doesn't hurt here
        source = [s for s in source if s != ""]
        target = [t for t in target if t != ""]
        resolution = [int(r) for r in resolution if r != ""]
        assert source and target and resolution, "{source}, {target} or {resolution} empty or malformed".format(**self.config)
        self.logger.debug("Returning from normal operation")
        return source, target, resolution
    def _region_generator(self, chrom, start, end, step, species, idcount):
        """Create a generator for genomic regions
        Generates strings of the form chrN\tstart\tend\tid\n
        Function assumes 0-based, half open intervals (compatible to UCSC)
        Regions are enumerated starting at 0
        """
        assert start < end, "_region_generator: defined invalid region {} -> {}".format(start, end)
        current_id = idcount
        for i in range(start, end, step):
            region_id = "@".join([species, chrom, str(current_id)])
            current_id = current_id + 1
            region = "\t".join([chrom, str(i), str(i+step), region_id])
            if i+step > end:
                region = "\t".join([chrom, str(i), str(end), region_id])
            region += "\n"
            yield region
    def _create_splitname(self, species, chrom, counter, resolution):
        filename = "_".join(["split",species,chrom,str(counter),str(resolution)])
        filename += ".bed"
        # outdir created in prepare_mapping, so guaranteed to exist
        filepath = os.path.join(self.config['outdir'],filename)
        return filepath
    def _generate_genome_tiling(self, source, resolution, chromsizes):
        """Function generates regularly spaced genome tiling
        Number of regions per file is set in the configuration file
        in the entry 'splitsize'
        """
        splitfile_list = []
        regions_per_file = int(self.config["splitsize"])
        for chrom_tuple in chromsizes:
            chrom_id = chrom_tuple[0]
            chrom_length = chrom_tuple[1]
            tmp = divmod(chrom_length, resolution * regions_per_file)
            num_of_files = tmp[0] if tmp[1] == 0 else tmp[0] + 1
            region_counter = 0
            # intervals 0-based, half open; UCSC coord. compliant
            regions_start = 0
            regions_end = min(regions_start + regions_per_file * resolution, chrom_length)
            for i in range(num_of_files):
                filebuf = io.StringIO()
                reggen = self._region_generator(chrom_id, regions_start, regions_end, resolution, source, region_counter)
                for line in reggen:
                    filebuf.write(line)
                # write buffer to file
                outfile = self._create_splitname(source, chrom_id, i, resolution)
                with open(outfile, "w") as of:
                    of.write(filebuf.getvalue())
                splitfile_list.append(outfile)
                region_counter += regions_per_file
                regions_start = regions_end
                regions_end = min(regions_start + regions_per_file * resolution, chrom_length)
        assert splitfile_list, "_generate_genome_tiling: list of filenames for genome tiling splits is empty"
        self.logger.debug("returning from normal operation")
        return splitfile_list
    def _create_lift_filename(self, source, target, inputfile, regiontype):
        """_create_lift_filename creates informative filenames
        for the mapped/unmapped output files of liftOver
        """
        old_filename = os.path.basename(inputfile)
        path = os.path.dirname(inputfile)
        assert "split" in old_filename, "_create_lift_filename: did not receive a genome tiling splitfile: {}".format(inputfile)
        parts = old_filename.split("_")
        # the above split at each "_" preserves the file extension .bed
        # filename = "_".join(["split",species,chrom,str(counter),str(resolution)]) this is how the name is generated in _create_splitname
        new_filename = "_".join([regiontype, source, target, parts[2], parts[3], parts[4]])
        return (os.path.join(path, new_filename))
    def _generate_instruction_file(self, source, target, input_files, chainfiles):
        filename = "00_lift_instr_" + source + ".txt"
        instrfile = os.path.join(self.config['outdir'], filename)
        # copy liftOver executable to /scratch
        liftpath = self.config['liftbin']
        liftpath = osi.copy_file_to_dir(liftpath, self.config['outdir'])
        # need [?] to include the "./" since the cwd is not in the PATH of the grid engine
        liftbin = os.path.join(os.path.dirname(liftpath), "./" + os.path.basename(liftpath))
        liftparams = self.config['liftparam']
        with open(instrfile, 'w') as of:
            for t in target:
                if t == source: continue
                key = "-".join([source, t])
                try:
                    DBentry = chainfiles[key]
                    # get chainfile and copy it to /scratch for performance reasons
                    chainfile = os.path.join(DBentry['dir'], DBentry['filename'])
                    chainfile = osi.copy_file_to_dir(chainfile, self.config['outdir'])
                    for infile in input_files:
                        mapfile = self._create_lift_filename(source, t, infile, "map")
                        unmapfile = self._create_lift_filename(source, t, infile, "unmap")
                        line = " ".join([liftbin, liftparams, infile, chainfile, mapfile, unmapfile])
                        of.write(line+"\n")
                except KeyError:
                    self.logger.warning("_generate_instruction_file: No chainfile found for key {}".format(key))
                    continue
        self.logger.debug("_generate_instruction_file: returning from normal operation")
        return instrfile
    # ==================================================
    # next use case: process blank mapping
    def process_mapping(self, dbman, toolbox):
        """process_mapping processes a previously created blank
        mapping with overlapping regions and reduces the mapping
        by selecting the fragments with the highest conservation (fragments are parts of larger clusters).
        The computation is based on the phastCon scores from UCSC, computed
        for the N species multiz alignment
        """
        status = True
        self.logger.debug("process_mapping started")
        mapfiles = dbman.create_DBdict("mapfiles")
        process_files = [mapfd for mapfd in mapfiles.values() if mapfd['level'] == "unprocessed" and mapfd['status'] == "use"]
        if not process_files:
            self.logger.info("process_mapping: no mapfiles to process - if wrong, check DB entries' status")
            return status
        # check if we have the correct cons scores for these mapfiles
        phastcon_files = None
        try:
            phastcon_files = dbman.create_DBdict("phastcons")
        except AssertionError as ae:
            # there are no DB entries for conservation scores
            status = False
            self.logger.error("Please setup DB with files holding appropriate conservation scores and restart\nDBmanager reported {}".format(ae))
            return status
        assert phastcon_files is not None, "process_mapping: No phastcon_files found, problem with DB?"
        dict_pairs = sorted(self._match_mapping_phastcon(process_files, phastcon_files), key=lambda t: t[0]['filename'])
        # dict_pairs is list of 2-tuples of DB dictionaries: tup[0] is phastcons_dict
        last_phastcons_file = ""
        last_phastcons_dict = None
        for p,m in dict_pairs:
            out_buffer = None
            if os.path.join(p['dir'],p['filename']) != last_phastcons_file:
                last_phastcons_file = os.path.join(p['dir'], p['filename'])
                last_phastcons_dict = self._build_phastcons_dict(last_phastcons_file)
            out_buffer = self._reduce_mapping(os.path.join(m['dir'], m['filename']), last_phastcons_dict, m['resolution'], toolbox)
            if out_buffer is None:
                self.logger.error("process_mapping found empty output buffer for mapping file {} and phastcons file {}".format(os.path.join(m['dir'],m['filename']), os.path.join(p['dir'],p['filename']))) 
                status = False
                continue
            self._write_output_buffer(self.config['outdir'], m['filename'], out_buffer)
        toolbox.clean_up()
        self.logger.debug("process_mapping: returning from normal operation")
        return status
    def _match_mapping_phastcon(self, mapfiles_list, phastconfiles_dicts):
        """_match_mapping_phastcon matches each mapping file with the correct
        phastcon file to allow region reduction based on conservation scores.
        Simply returns list of matching tuples of DB dictionaries (mapfile dict, phastconfile dict)
        """
        tuple_list = []
        for mapfdict in mapfiles_list:
            for phcfdict in phastconfiles_dicts.values():
                if mapfdict['source_assembly'] == phcfdict['assembly'] and mapfdict['mapping_type'] == phcfdict['range_covered'] and mapfdict['resolution'] == phcfdict['resolution']:
                    tuple_list.append((phcfdict, mapfdict))
        assert tuple_list, "_match_mapping_phastcon: No matching file with cons. scores for any map file"
        if len(tuple_list) != len(mapfiles_list):
            self.logger.warning("some map files not associated with cons. score files\n{}\n{}".format(mapfiles_list, tuple_list))
            return tuple_list
        self.logger.debug("returning from normal operation")
        return tuple_list
    def _write_output_buffer(self, outdir, outname, outbuffer):
        """ Write the output buffer to the new location (same filename as before)
        and call bedSort on the reduced mapping (sanity check)
        """
        with open(os.path.join(outdir, outname), "w") as outfile:
            outfile.write(outbuffer.getvalue())
            outfile.flush()
            os.fsync(outfile.fileno())
        # DEBUG: use the toolbox here
        #bedsort = self.config['bedsort'
        #p = sp.Popen([bedsort, outfile.name, outfile.name])
        #ret = p.wait()
        return
    def _build_phastcons_dict(self, filename):
        cons_dict = {}
        with open(filename, "r") as infile:
            for line in infile:
                cols = line.strip().split("\t")
                last_at = cols[3].rfind("@")
                region_id = (cols[3])[:last_at]
                cons_value = float((cols[3])[last_at+1:])
                cons_dict[region_id] = cons_value
        return cons_dict
    def _split_self_intersection(self, fobj, phastconsdict):
        good_buffer = io.StringIO()
        ovllist = []
        for line in fobj:
            cols = line.strip().split("\t")
            region_id = cols[3]
            consscore = 0.0
            try:
                consscore = phastconsdict[region_id]
            except KeyError: continue
            ovls = int(cols[4])
            ext_id = region_id + "@" + str(consscore)
            outline = "\t".join([cols[0], cols[1], cols[2], ext_id])
            if ovls == 1: # region only overlaps with itself
                good_buffer.write(outline)
                good_buffer.write("\n")
            else:
                regid_parts = region_id.split("@")
                reg_chrom = regid_parts[1]
                reg_id = int(regid_parts[2])
                region = {'seq_chrom':cols[0], 'seq_start':int(cols[1]), 'seq_end':int(cols[2]), 'reg_chrom':reg_chrom, 'reg_id':reg_id, 'conservation':consscore, 'output_line':outline}
                ovllist.append(region)
            continue
        self.logger.debug("_split_self_intersection: returning from normal operation")
        return good_buffer, ovllist
    def _intersect_files(self, file1, file2):
        """Intersect the mapping file with itself and count the number
        of overlaps for each region. Write the output to a temp file
        that is deleted later
        """
        tmp_out = tf.NamedTemporaryFile(mode="w", buffering=-1, suffix=".tmp", prefix="Py_", dir="/tmp", delete=False)
        isectbed = self.config['intersectbed']
        p = sp.Popen([isectbed, "-c", "-a", file1, "-b", file2], stdout=tmp_out.fileno())
        ret = p.wait()
        tmp_out.close()
        self.logger.debug("_intersect_files: returning from normal operation")
        return tmp_out
    def _reduce_mapping(self, mapfile, phastconsdict, resolution, toolbox):
        """_reduce_mapping reduces a mapfile with overlapping regions based on
        conservation scores (pertain to region with higher conservation)
        Returns file-like object with non-overlapping regions
        """
        self.logger.info("Processing file {}...".format(mapfile))
        # step 1: self-intersection of mapping file to find multiple overlapping regions
        isect_fobj = toolbox.scl_bed_intersect_count(mapfile, mapfile, "TEMPFILE")
        # step 2: split the self-intersection into 1 [self] and n [mult] overlaps
        # and create list of overlapping regions
        good_buffer, ovllist = self._split_self_intersection(open(isect_fobj.name, "r"), phastconsdict)
        isect_fobj.close()
        if not ovllist:
            return good_buffer
        # step 3: iterate through list of overlapping regions and build clusters
        self.logger.debug("Building clusters...")
        cluster_dict = self._build_clusters(ovllist, resolution)
        self.logger.info("Identified {} clusters based on {} overlapping regions".format(len(cluster_dict), len(ovllist)))
        # step 4: determine overlaps among clusters
        self.logger.debug("Computing overlaps among clusters...")
        conslst, cluster_overlap_dict = self._compute_overlaps(cluster_dict)
        self.logger.debug("Starting cluster selection...")
        # step 5: select clusters based on conservation, start with cluster with fewest overlaps
        selected_clusters = self._select_clusters(conslst, cluster_overlap_dict, cluster_dict)
        self.logger.info("Selected {} non-overlapping clusters".format(len(selected_clusters)))
        for memid in selected_clusters:
            cluster = cluster_dict[memid]
            good_buffer.write(cluster.get_cluster_regions())
            good_buffer.write("\n")
        self.logger.debug("Returning from normal operation")
        return good_buffer
    def _build_clusters(self, overlapping_regions, resolution):
        """Iterate through list of overlapping regions
        and build clusters - each region can only be
        part of one cluster
        Also compute cluster conservation
        """
        cluster_dicts = col.defaultdict(dict)
        overlapping_regions = sorted(overlapping_regions, key=lambda s: s['seq_chrom'])
        last_seqchrom = ""
        last_clusterdict = None
        for region in overlapping_regions:
            if region['seq_chrom'] != last_seqchrom:
                if last_clusterdict:
                    self.logger.debug("Building clusters for chromosome {}".format(last_seqchrom))
                    cluster_dicts[last_seqchrom] = last_clusterdict
                last_seqchrom = region['seq_chrom']
                last_clusterdict = cluster_dicts[region['seq_chrom']]
            for cluster in last_clusterdict.values():
                if cluster.is_adjacent(region):
                    cluster.add_region(region)
                    break
            else:
                # inner for loop terminated normally (no break)
                # no adjacency found, create new cluster
                newcluster = RegionCluster(resolution, region)
                last_clusterdict[id(newcluster)] = newcluster
                #cluster_dict[id(newcluster)] = newcluster
        self.logger.debug("Adding last clusters for chromosome {}".format(last_seqchrom))
        cluster_dicts[last_seqchrom] = last_clusterdict
        self.logger.debug("Collapsing all clusters...")
        cluster_dict = {}
        for k in cluster_dicts.keys():
            for a,b in cluster_dicts[k].items():
                cluster_dict[a] = b
        self.logger.debug("All clusters built, compute conservation")
        for c in cluster_dict.values():
            c.compute_conservation()
        self.logger.debug("Returning all clusters")
        return cluster_dict
    def _compute_overlaps(self, cluster_dict):
        """Compute overlaps among all clusters and return sorted list
        containing those clusters with highest conservation first and a dictionary
        of lists storing the overlapping clusters (IDs) for each cluster
        """
        store_ovl = col.defaultdict(list)
        for memid, cluster in cluster_dict.items():
            for memid2, cluster2 in cluster_dict.items():
                if memid == memid2: continue
                if cluster.overlaps(cluster2):
                    store_ovl[memid].append(memid2)
        tmplst = [ ( clust.cluster_cons, clustid ) for (clustid, clust) in cluster_dict.items()  ]
        conslst = sorted(tmplst, reverse=True) # highest conservation first
        return conslst, store_ovl
    def _select_clusters(self, conslst, ovldict, cluster_dict):
        """ Select those cluster with highest conservation score
        Results in set of non-overlapping clusters
        """
        forbidden = set()
        selected = set()
        for cons, memid in conslst:
            if memid in forbidden: continue
            if memid in selected:
                forbidden.update(ovldict[memid]) # just to be sure
                continue
            # in theory, comparing the conservation among
            # all overlapping clusters should not be necessary
            # since we start with the highest - in practice, let's
            # do it the safe way and double-check
            select_cluster = cluster_dict[memid]
            select_id = memid
            for other_id in ovldict[memid]:
                if other_id in forbidden: continue
                if other_id in selected:
                    forbidden.add(memid) # that should always already be in there
                    select_id = other_id
                    break
                checked_cluster = cluster_dict[other_id]
                checked_id = other_id
                if checked_cluster.cluster_cons > select_cluster.cluster_cons:
                    forbidden.add(select_id)
                    select_cluster = checked_cluster
                    select_id = checked_id
                    continue
                else:
                    forbidden.add(checked_id)
                    continue
            selected.add(select_id)
            forbidden.update(ovldict[select_id])
        return selected
    # ==================================
    # next usecase: perform mapping of wig files
    def map_signal_tracks(self, dbman, toolbox):
        """Perform the actual mapping of signal tracks from source to target species
        Trigger is set to convert resulting signal track to regular resolution again
        (i.e. 25, 50, 100 etc.) or to convert to same resolution as input
        """
        status = True
        self.logger.debug("Start mapping of signal tracks")
        wigfiles_dict = dbman.create_DBdict("signaltracks")
        mapfiles_dict = dbman.create_DBdict("mapfiles")
        # TODO: possible change here to avoid changing status - map all files from "species" to "species"
        wigfiles_lst = [ wigfd for wigfd in wigfiles_dict.values() if wigfd['status'] == 'use' ]
        mapfiles_lst = [ mapfd for mapfd in mapfiles_dict.values() if mapfd['level'] == 'processed' and mapfd['status'] == 'use' ]
        if not (wigfiles_lst or mapfiles_lst):
            self.logger.info("No wig files to map or no map files to process - if unexpected, check DB status of files")
            status = False
            return status
        match_list = self._match_mapfile_signalfiles(mapfiles_lst, wigfiles_lst)
        new_wigfiles = []
        for mapfd,wiglst in match_list:
            map_dict = self._build_mapping_dict(mapfd)
            for wigfd in wiglst:
                mapped_file = self._map_signal_track(map_dict, mapfd, wigfd, toolbox)
                new_wigfiles.append(mapped_file)
        self.logger.debug("Updating DB with following entries: {}".format(new_wigfiles))
        dbman.update_DBfile_make_id("signaltracks", wigfiles_dict, new_wigfiles)
        toolbox.clean_up()
        self.logger.debug("Mapping of signal tracks finished")
        return status       
    def _match_mapfile_signalfiles(self, mapfiles_lst, wigfiles_lst):
        """Match each mapfile with corresponding signal tracks
        and return list of tuples
        """
        res_lst = []
        for mapfd in mapfiles_lst:
            matched_wigfiles = []
            map_source = mapfd['source_assembly']
            map_res = mapfd['resolution']
            for wigfd in wigfiles_lst:
                if map_source == wigfd['assembly'] and map_res == wigfd['resolution']:
                    matched_wigfiles.append(wigfd)
            if not matched_wigfiles:
                self.logger.warning("Could not find any matching wiggle tracks for mapping file {}".format(mapfd['filename']))
                continue
            res_lst.append( (mapfd, matched_wigfiles) )
        assert res_lst, "Could not find matching wiggle tracks for any mapping file - check DB entries"
        self.logger.debug("Combined the following files:\n{}".format(res_lst))
        return res_lst
    def _build_mapping_dict(self, mapfd):
        """Builds the dictionary to map tracks from one
        assembly to another. Results in large objects in
        memory (several GB)
        """
        with open(os.path.join(mapfd['dir'], mapfd['filename']), "r") as infile:
            mapdict = {}
            for line in infile:
                right_tab = line.rfind("\t")
                right_at = line.rfind("@")
                # this is chrom\tstart\tend\t
                value = line[:right_tab + 1]
                # this is assembly@chrom@regioncounter
                key = line[right_tab + 1 : right_at]
                mapdict[key] = value
        self.logger.debug("Mapping dictionary {} - {} created with {} entries".format(mapfd['source_assembly'], mapfd['target_assembly'], len(mapdict)))
        return mapdict
    def _compute_consecutive_regions(self, fobj, resolution):
        """Create a buffer with all consecutive regions from
        the mapped signal track to be send to the smoother

        CRITICAL: coordinate system used for interpolation
        UCSC compliant 0-based, half-open BED regions
        Example: region 100-150 has data values at positions 101-150
        """
        output_buffer = io.StringIO()
        with open(fobj.name, "r") as f:
            xvals = []
            yvals = []
            last_chrom = ""
            for line in f:
                chrom, start, end, value = line.strip().split("\t")
                start = int(start)
                end = int(end)
                # input file is sorted, so start of next region is always larger than end of last region
                # Example: start of next 1093 / end of last 987 -> gap too large
                if chrom == last_chrom and start - xvals[-1][1] < resolution:
                    xvals.append( (start, end) )
                    yvals.append( float(value) )
                    continue
                else:
                    output = self._extract_region(last_chrom, xvals, yvals, resolution)
                    output_buffer.write(output)
                    last_chrom = chrom
                    xvals = []
                    yvals = []
                    # CRITICAL: use only positions with data values
                    xvals.append( (start, end) )
                    yvals.append( float(value) )
                    continue
            if xvals:
                output = self._extract_region(last_chrom, xvals, yvals, resolution)
                output_buffer.write(output)
        self.logger.debug("Identified all consecutive regions")
        return output_buffer
    def _extract_region(self, chrom, xvals, yvals, resolution):
        """Extract the well-defined region (proper start - end)
        and return it as a string ready for the smoother
        """
        if chrom == "": return ""
        yval_range = []
        for it in range(len(xvals)):
            # xvals contains tuples (start, end) as 0-based, half-open intervals
            # +1 to simplify calculations, region 100-150 has length 51 here
            length = (xvals[it][1] - xvals[it][0]) + 1
            try:
                if xvals[it][1] == xvals[it + 1][0]:
                    # consecutive, correct for double-counting of boundary position
                    yval_range.extend( [yvals[it]] * (length - 1))
                    continue
                else:
                    yval_range.extend( [yvals[it]] * length)
                    gap_length = (xvals[it + 1][0] - xvals[it][1]) - 1
                    gap_value = (yvals[it] + yvals[it + 1]) / 2.
                    yval_range.extend( [gap_value] * gap_length)
                    continue
            except IndexError:
                # at the last tuple
                yval_range.extend( [yvals[it]] * length)
                continue
        xval_range = range(xvals[0][0], xvals[-1][1] + 1, 1)
        assert len(xval_range) == len(yval_range), "Missed a gap while formatting regions for smoothing: x {} - y {}".format(len(xval_range), len(yval_range))
        startidx = self._find_index(xval_range, resolution, False)
        endidx = self._find_index(xval_range, resolution, True)
        # note here that if startidx and endidx are -1 -> list is empty
        xval_range = xval_range[startidx : (endidx + 1)]
        yval_range = yval_range[startidx : (endidx + 1)]
        if len(xval_range) < resolution: return ""
        xval_range = " ".join(map(str, xval_range))
        yval_range = " ".join(map(str, yval_range))
        s = "\t".join([chrom, xval_range, yval_range])
        return s + "\n"
    def _find_index(self, sequence, modvalue, reverse):
        if reverse:
            for item in reversed(sequence):
                if item % modvalue == 0:
                    return sequence.index(item)
        for item in sequence:
            if item % modvalue == 0:
                return sequence.index(item)
        return -1
    def _run_smoother(self, file_loc, res):
        """Call to R script for lowess smoothing
        """
        free_cpus = osi.estimate_free_cpus()
        self.logger.info("Estimated {} free CPUs".format(free_cpus))
        # make sure everything is str
        args = [str(v) for v in [self.config['smoother'], file_loc, res, free_cpus, self.config['fraction'], self.config['iteration']]]
        self.logger.debug("Calling external smoother with args: {}".format(args))
        p = sp.Popen(args)
        ret = p.communicate()[0]
        self.logger.debug("Smoother returned from operation")
        return ret
    def _create_meta_data(self, mapfd, wigfd, resolution):
        """Create a dictionary of meta data pertaining to the newly created
        signal track - is used to update the DB file
        """
        new_dict = wigfd
        new_dict['timestamp'] = dat.date.today().isoformat()
        new_dict['dir'] = self.config['outdir']
        new_count = None
        if 'MAPDFROM' in wigfd['link']:
            # file has been mapped at least twice
            parts = wigfd['filename'].split("_")
            new_count = str(int(parts[0].strip("x")) + 1) + "x"
        else:
            new_count = "1x"
        filename = "_".join([new_count, mapfd['source_assembly'], mapfd['target_assembly'], wigfd['type'], wigfd['cell'], str(resolution)])
        new_dict['filename'] = filename + ".wig"
        # currently only bedgraph supported
        new_dict['fileformat'] = self.config['outformat']
        # reasonable default?
        new_dict['status'] = 'ignore'
        new_dict['resolution'] = str(resolution)
        new_dict['species'] = mapfd['target_species']
        new_dict['assembly'] = mapfd['target_assembly']
        new_dict['link'] = 'MAPDFROM_' + wigfd['key']
        new_dict['p20'] = "-1"
        new_dict['p40'] = "-1"
        new_dict['p60'] = "-1"
        new_dict['p80'] = "-1"
        new_dict['p90'] = "-1"
        new_dict['p95'] = "-1"
        #s = "{timestamp}\t{dir}\t{filename}\t{fileformat}\t{status}\t{resolution}\t{species}\t{assembly}\t{type}\t{cell}\t{category}\t{lab}\t{source}\t{link}\t{p20}\t{p40}\t{p60}\t{p80}\t{p90}\t{p95}".format(**new_dict)
        self.logger.debug("Created meta information for mapped track")
        return new_dict, os.path.join(new_dict['dir'], new_dict['filename'])
    def _map_signal_track(self, map_dict, mapfd, wigfd, toolbox):
        """Perform mapping and clean-up (sorting, conversion
        to regular resolution, generate info for new file)
        """
        outbuffer = io.StringIO()
        map_count = 0
        line_count = 0
        verrs = 0
        kerrs = 0
        with open(os.path.join(wigfd['dir'], wigfd['filename']), "r") as infile:
            wig_assembly = wigfd['assembly']
            wig_res = int(wigfd['resolution'])
            for line in infile:
                line_count += 1
                try:
                    chrom, start, end, value = line.strip().split()
                    #chrom, start, end, value = line.strip().split("\t")
                    region_count = str(int(start) // wig_res)
                    region_key = "@".join([wig_assembly, chrom, region_count])
                    target_region = map_dict[region_key]
                    outline = target_region + value + "\n"
                    outbuffer.write(outline)
                    map_count += 1
                except ValueError as ve:
                    verrs += 1
                    continue # excepts header lines
                except KeyError:
                    kerrs += 1
                    continue # missing keys may happen
        # file is mapped, next is sorting
        self.logger.info("Mapping of file {} complete, start sorting and signal interpolation...".format(wigfd['filename']))
        self.logger.debug("Mapped a total of {} out of {} regions".format(map_count, line_count))
        self.logger.debug("Recorded {} value errors and {} key errors".format(verrs, kerrs))
        tmp_out = tf.NamedTemporaryFile(mode="r+", buffering=-1, suffix=".tmp", prefix="Py_", dir="/tmp", delete=False)
        tmp_out.write(outbuffer.getvalue())
        tmp_out.flush()
        os.fsync(tmp_out.fileno())
        outbuffer = io.StringIO()
        tmp_out = toolbox.scl_bed_sort(tmp_out, tmp_out)
        out_res = int(wigfd['resolution'] if self.config['interpolate'] == 'input' else self.config['interpolate'])
        outbuffer = self._compute_consecutive_regions(tmp_out, out_res)
        tmp_out = open(tmp_out.name, "w")
        tmp_out.seek(0)
        tmp_out.write(outbuffer.getvalue())
        tmp_out.close()
        # generate the meta data for the mapped track
        new_db_entry, new_file = self._create_meta_data(mapfd, wigfd, out_res)
        # smooth signal
        smooth_ret = self._run_smoother(tmp_out.name, out_res)
        self.logger.info("Smoothing finished with status {}".format(smooth_ret))
        # mv file to new location
        osi.move_file_to_newfile(tmp_out.name, new_file)
        self.logger.debug("Mapping of track {} done".format(wigfd['filename']))
        return new_db_entry
