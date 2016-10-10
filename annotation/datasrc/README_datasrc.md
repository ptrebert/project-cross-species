# Selected source files

This README merely lists the selection criteria. Not all of the files matching these criteria are actually processed.

## Source: ENCODE portal (encodeproject.org)

Date: 2016-02-11

* Assay
    * ChIP-seq
    * DNase-seq
    * RNA-seq
* Project
    * ENCODE
* Genome assembly
    * hg19
    * mm9
* Biosample type
    * immortalized cell line
    * stem cell
* Available data
    * bigWig
    * bigBed narrowPeak
    * bigBed broadPeak
    * gtf
* Run type
    * single-ended
* Read length (nt)
    * 36
* Replication type
    * isogenic
    

Note 2016-04-29

* Some histone peak files contain spurious peaks (peak size as low as 1bp)
    * Example: ENCFF001MXB [entry: 1 -- chr8 - 11007712 - 11007713]
    * Filter out these peaks, i.e. require peaks to be at least "nucleosome-sized"
    
Note 2016-05-24

* Discovered that mouse (mm9) expression data is mostly listed under assay type "polyA mRNA RNA-seq"
 as opposed to the regular "RNA-seq"; though, in the metadata file, the assay is given as "RNA-seq";
 downloaded additional metadata file for mm9/expression (single ended, all read lengths)
 
Note 2016-05-25

* The expression data for hg19 (SE libraries) are missing quantifications beyond the exon level; downloaded
 additional metadata/files (PE libraries) to get expression estimates on gene level
 
Note 2016-07-13
 
 * changed approach since gene-level expression quantification is anyway required; downloaded all ENCODE
 file metadata and created normalized annotation file from that (encode_metadata_ro.tsv)