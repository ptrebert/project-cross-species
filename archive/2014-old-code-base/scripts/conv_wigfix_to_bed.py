#!/usr/bin/python

"""Script to convert fixed step wiggle tracks
into bed files with a fixed resolution.
Currently, it does not make sense to include this
in the cross-species epigenomics framework.

Converter class expects an input file and an output file
(as readable/writeable objects), a resolution

Assumes a step size of 1 for the wiggle input
The wiggle format [coordinates] are 1-based, i.e. a chromosome of length N has positions 1...N
Output BED will be 0-based, half open intervals (UCSC compliant)
"""

import sys as sys
import os as os
import re as re
import cStringIO as sio
import gzip as gz
# OptParse is deprecated in python2.7+
import optparse as opp


class WigFixToBed(object):
    """WigFixToBed conversion - could be derived
    from generic converter class if embedded in framework
    """
    def __init__(self, cfgobj):
        self.inputfile = cfgobj.inputfile
        self.outputfile = cfgobj.outputfile
        self.resolution = cfgobj.resolution
        self.aggregate = cfgobj.aggregate
        self.build = cfgobj.build
        assert self.inputfile is not None and self.outputfile is not None and self.aggregate is not None and self.build is not None, "Converter setup: at least one setting is None"
        assert self.resolution > 0, "Converter setup: resolution is invalid, given as %s" %self.resolution
    def start(self):
        # AVG assumes step size of 1
        agg_dict = { "MIN":min, "MAX":max, "AVG":(lambda lst: sum(lst)/len(lst)) }
        aggregation_method = agg_dict[self.aggregate]
        # check if we have already file-like object
        out_buffer = None
        try:
            self.inputfile.readline()
            self.inputfile.seek(0)
            out_buffer = self._convert(self.inputfile, sio.StringIO(), aggregation_method)
        except AttributeError:
            wigfile = None
            if self.inputfile[-3:] == ".gz":
                wigfile = gz.open(self.inputfile, "rb")
            else:
                wigfile = open(self.inputfile, "r")
            out_buffer = self._convert(wigfile, sio.StringIO(), aggregation_method)
            wigfile.close()
        assert out_buffer is not None, "Output buffer is None"
        ofobj = None
        try:
            self.outputfile.write("")
            self.outputfile.seek(0)
            ofobj = self.outputfile
        except AttributeError:
            ofobj = open(self.outputfile, "w")
        ofobj.write(out_buffer.getvalue())
        ofobj.flush()
        os.fsync(ofobj.fileno())
        ofobj.close()
        return None
    def _convert(self, ifobj, ofobj, aggmethod):
        match_chrom = re.compile("(?<=chrom\=)\w+(?=\s)")
        match_start = re.compile("(?<=start\=)\d+(?=\s)")
        AGGREGATE = False
        position_counter = 0
        start = 0
        current_chrom = ""
        values = []
        for line in ifobj:
            try:
                consvalue = float(line)
                if len(values) == self.resolution:
                    region_count = start / self.resolution
                    aggval = str(aggmethod(values))
                    region_id = "@".join([self.build, current_chrom, str(region_count), aggval])
                    outline = "\t".join([current_chrom, str(start), str(start + self.resolution), region_id])
                    ofobj.write(outline + "\n")
                    values = [consvalue]
                    start = start + self.resolution
                    position_counter += 1
                    continue
                elif AGGREGATE:
                    values.append(consvalue)
                    position_counter += 1
                    continue
                elif position_counter % self.resolution == 0:
                    AGGREGATE = True
                    start = position_counter
                    position_counter += 1
                    continue
                else:
                    position_counter += 1
                    continue
            except ValueError:
                (chrom, block_start) = self._parse_header(line, match_chrom, match_start)
                if chrom == current_chrom and block_start == position_counter:
                    continue
                else:
                    if len(values) == self.resolution:
                        region_count = start / self.resolution
                        aggval = str(aggmethod(values))
                        region_id = "@".join([self.build, current_chrom, str(region_count), aggval])
                        outline = "\t".join([current_chrom, str(start), str(start + self.resolution), region_id])
                        ofobj.write(outline + "\n")
                    current_chrom = chrom
                    values = []
                    if block_start % self.resolution == 0:
                        position_counter = block_start + 1
                        start = block_start
                        AGGREGATE = True
                        next(ifobj)
                        continue
                    elif (block_start - 1) % self.resolution == 0:
                        start = block_start - 1
                        position_counter = block_start
                        AGGREGATE = True
                        continue
                    else:
                        AGGREGATE = False
                        position_counter = block_start
                        continue
        if len(values) == self.resolution:
            region_count = start / self.resolution
            aggval = str(aggmethod(values))
            region_id = "@".join([self.build, current_chrom, str(region_count), aggval])
            outline = "\t".join([current_chrom, str(start), str(start + self.resolution), region_id])
            ofobj.write(outline + "\n")
        return ofobj
    def _parse_header(self, line, match_chrom, match_start):
        assert "fixedStep" in line and "step=1" in line, "Found this header %s - not a fixed step wiggle track with step size 1" %line
        chrom_mobj = match_chrom.search(line)
        start_mobj = match_start.search(line)
        assert chrom_mobj and start_mobj, "Could not find chrom or start in this line %s" %line
        chrom = chrom_mobj.group(0)
        start = int(start_mobj.group(0))
        return (chrom, start)       

if __name__ == "__main__":
    try:
        parser = opp.OptionParser()
        parser.add_option("-i", "--infile", dest="inputfile", type="string", help="Path to input file including full file name (if ending is .gz, gzipped file is assumed)")
        parser.add_option("-r", "--resolution", dest="resolution", type="int", help="Resolution [INTEGER] of the output BED file (in bp)")
        parser.add_option("-o", "--outfile", dest="outputfile", type="string", help="Path to desired location for output file including full file name")
        parser.add_option("-a", "--aggregate", dest="aggregate", type="string", help="Aggregation method (MIN, MAX, AVG) - AVG is average value per bp, i.e. normalized by length of region")
        parser.add_option("-b", "--build", dest="build", type="string", help="Genome build of species, e.g. hg19, mm10, danrer7")
        # Setup ready
        (options, args) = parser.parse_args()
        converter = WigFixToBed(options)
        converter.start()
        sys.stdout.write("\nConversion finished, output should be in %s\nexiting...\n" %options.outputfile)
    except Exception as e:
        sys.stderr.write("\nTerminated with exception %s, exiting...\n" %e)
        sys.exit(1)
    else:
        sys.exit(0)
