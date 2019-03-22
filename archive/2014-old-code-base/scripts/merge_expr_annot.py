#!/usr/bin/python3


import sys as sys
import os as os
import io as io
import collections as col
import fnmatch as fnm
import copy as co
import numpy as np

def parse_promoter_files(lcg_loc, hcg_loc):
    header = ['RefSeq_mRNA', 'chrom', 'strand', 'TSS_loc', 'GC_cont', 'norm_GC_cont', 'CpG_frac', 'length']
    resdict = {}
    with open(lcg_loc, 'r') as infile:
        _ = infile.readline()
        for line in infile:
            cols = line.strip().split(',')
            rid = cols[0]
            tmp = dict(zip(header, cols))
            tmp['promoter_class'] = 'LCG'
            resdict[rid] = tmp
    with open(hcg_loc, 'r') as infile:
        _ = infile.readline()
        for line in infile:
            cols = line.strip().split(',')
            rid = cols[0]
            tmp = dict(zip(header, cols))
            tmp['promoter_class'] = 'HCG'
            resdict[rid] = tmp
    print('Parsed promoter files')
    return resdict

def parse_refseq_ensembl_mapping(file_loc):
    header = ['ensembl_gene_id', 'ensembl_transcript_id', 'chrom', 'strand', 'GO_acc', 'GO_name', 'RefSeq_mRNA']
    resdict = {}
    with open(file_loc, 'r') as infile:
        _ = infile.readline() # skip header
        for line in infile:
            cols = line.strip().split('\t')
            rid = cols[6]
            resdict[rid] = dict(zip(header, cols))
    print('Parsed mapping')
    return resdict

def map_refseq_to_ensembl(lcg, hcg, biom):
    promoter_dict = parse_promoter_files(lcg, hcg)
    biomart_mapping = parse_refseq_ensembl_mapping(biom)
    result = {}
    for k,v in promoter_dict.items():
        if k in biomart_mapping:
            ensdict = biomart_mapping[k]
            ensdict['promoter_class'] = v['promoter_class']
            result[ensdict['ensembl_transcript_id']] = ensdict
    print('Merged annotations')
    return result

def parse_bedfile(file_loc):
    perc_cpg = np.array([], dtype=float)
    beddict = {}
    with open(file_loc, 'r') as infile:
        header = infile.readline().strip('#\n').split('\t')
        for line in infile:
            cols = line.strip().split('\t')
            tmp = dict(zip(header, cols))
            perc_cpg = np.append(perc_cpg, float(tmp['perc_CpG']))
            beddict[tmp['ensembl_transcript_id']] = tmp
    low_cpg, high_cpg = np.percentile(perc_cpg, [30., 70.])
    insertidx = header.index('label')
    header.insert(insertidx, 'promoter_class')
    header.insert(insertidx+1, 'RefSeq_mRNA')
    return beddict, header, low_cpg, high_cpg

def parse_annotation(file_loc):
    genedict = {}
    transcriptdict = {}
    gene_keys = ['ensembl_gene_id', 'gene_type', 'gene_name', 'havana_gene_id']
    gene_values = ['gene_id', 'gene_type', 'gene_name', 'havana_gene']
    transcript_keys = ['ensembl_gene_id', 'ensembl_transcript_id', 'transcript_name', 'transcript_type', 'gene_name', 'havana_gene_id', 'havana_transcript_id']
    transcript_values = ['gene_id', 'transcript_id', 'transcript_name', 'transcript_type', 'gene_name', 'havana_gene', 'havana_transcript']
    with open(file_loc, 'r') as annotfile:
        for line in annotfile:
            cols = line.strip().split("\t")
            kv_part = cols[-1].split(";")
            basic_dict = dict(zip(['chrom', 'start', 'end', 'source', 'element', 'strand'], [cols[0], cols[3], cols[4], cols[1], cols[2], cols[6]]))
            if basic_dict['element'] not in ["gene", "transcript"]: continue
            kv_pairs = col.defaultdict(lambda: "NOTAVAILABLE")
            for kv in kv_part:
                if not kv or 'level' in kv or 'tag' in kv: continue
                k,v = kv.strip().split()
                k = k.strip()
                v = v.strip(' \"')
                kv_pairs[k] = v
            if basic_dict['element'] == "gene" and kv_pairs['gene_type'] == "protein_coding":
                basic_dict.update(dict(zip(gene_keys, [ kv_pairs[v] for v in gene_values ])))
                genedict[basic_dict['ensembl_gene_id']] = basic_dict
            elif basic_dict['element'] == "transcript" and kv_pairs['transcript_type'] == "protein_coding":
                basic_dict.update(dict(zip(transcript_keys, [ kv_pairs[v] for v in transcript_values ])))
                transcriptdict[basic_dict['ensembl_transcript_id']] = basic_dict
    return genedict, transcriptdict

def parse_input_dir(inputdir):
    assembly = "hg19"
    lab = "caltech"
    gencode = "Genc3c"
    all_files = os.listdir(inputdir)
    all_files = fnm.filter(all_files, "*.gtf")
    all_files = [ os.path.join(inputdir, f) for f in all_files ]
    res = []
    for f in all_files:
        celltype = "H1hesc" if "H1hesc" in f else "Gm12878"
        insertsize = "200" if "200" in f else "400"
        element = "TSS" if "TSS" in f else "gene"
        idx = f.find("Rep")
        replicate = f[idx:idx+4]
        fdict = {'assembly':assembly, 'lab':lab, 'gencode':gencode, 'fileloc':f, 'celltype':celltype, 'insertsize':insertsize, 'element':element, 'replicate':replicate}
        res.append(fdict)
    return res

def merge_expr_annot(all_files, genedict, transcriptdict, outdir):
    for fdict in all_files:
        if fdict['element'] == 'gene':
            continue
            print("Processing file {}".format(fdict['fileloc']))
            outbuf, outbuf_ext = quantify_gene(fdict, genedict)
            regfilename = "_".join([fdict['assembly'], fdict['celltype'], fdict['insertsize'], fdict['replicate'], fdict['element'], fdict['gencode']])
            regfilename = os.path.join(outdir, regfilename + ".bed")
            with open(regfilename, "w") as outfile:
                outfile.write(outbuf.getvalue())
                outfile.flush()
            regfilename = "_".join([fdict['assembly'], fdict['celltype'], fdict['insertsize'], fdict['replicate'], fdict['element'], "ext3kb", fdict['gencode']])
            regfilename = os.path.join(outdir, regfilename + ".bed")
            with open(regfilename, "w") as outfile:
                outfile.write(outbuf_ext.getvalue())
                outfile.flush()
        else:
            print("Processing file {}".format(fdict['fileloc']))
            outbuf = quantify_transcript(fdict, transcriptdict)
            regfilename = "_".join([fdict['assembly'], fdict['celltype'], fdict['insertsize'], fdict['replicate'], fdict['element'], 'ext5kb', fdict['gencode']])
            regfilename = os.path.join(outdir, regfilename + ".bed")
            with open(regfilename, "w") as outfile:
                outfile.write(outbuf.getvalue())
                outfile.flush()
    return
            
def quantify_gene(fdict, genedict):
    assertion = ['chrom', 'start', 'end', 'strand']
    header = ['chrom', 'start', 'end', 'strand', 'ensembl_gene_id', 'RPKM', 'gene_name', 'source', 'havana_gene_id', 'element', 'gene_type']
    outbuffer = io.StringIO()
    outbuffer.write("#" + "\t".join(header) + "\n")
    outbuffer_ext = io.StringIO()
    outbuffer_ext.write("#" + "\t".join(header) + "\n")
    with open(fdict['fileloc'], 'r') as exprfile:
        for line in exprfile:
            if not "protein_coding" in line: continue
            cols = line.strip().split("\t")
            kv_part = cols[-1].split(";")
            k,v = kv_part[0].split()
            ens_gene_id = v.strip('\"')
            if ens_gene_id not in genedict: continue
            annotentry = co.deepcopy(genedict[ens_gene_id])
            assert [cols[0], cols[3], cols[4], cols[6]] == [ annotentry[k] for k in assertion ], "\nFound gene mismatch for gene {} and location\n{}\n\n{}".format(ens_gene_id, line.strip(), annotentry)
            annotentry['RPKM'] = cols[5] # expr. value
            outline = "\t".join([annotentry[k] for k in header])
            outbuffer.write(outline + "\n")
            start = int(annotentry['start'])
            start = start - 3000
            assert start > 0, "Gene at the boundary: {}".format(line.strip())
            annotentry['start'] = str(start)
            annotentry['end'] = str(int(annotentry['end']) + 3000)
            outline = "\t".join([annotentry[k] for k in header])
            outbuffer_ext.write(outline + "\n")
    return outbuffer, outbuffer_ext

def quantify_transcript(fdict, transcriptdict):
    assertion = ['chrom', 'start', 'end', 'strand']
    header = ['chrom', 'start', 'end', 'strand', 'ensembl_gene_id', 'ensembl_transcript_id', 'RPKM', 'gene_name', 'source', 'havana_gene_id', 'havana_transcript_id', 'element', 'transcript_type']
    outbuffer = io.StringIO()
    outbuffer.write("#" + "\t".join(header) + "\n")
    with open(fdict['fileloc'], 'r') as exprfile:
        for line in exprfile:
            cols = line.strip().split('\t')
            kv_part = cols[-1].split(';')
            _, ens_id = kv_part[0].split()
            ensembl_gene_id = ens_id.strip('\"')
            _, trans_ids = kv_part[1].split()
            transcript_ids = trans_ids.strip('\"').split(',')
            for tid in transcript_ids:
                if tid not in transcriptdict: continue
                if cols[5] == 'nan': continue
                annotentry = co.deepcopy(transcriptdict[tid])
                annotentry['RPKM'] = cols[5] # expr. value
                start = int(cols[3]) - 3000 # TSS position
                end = int(cols[3]) + 2000
                assert start > 0, 'Transcript at boundary: {}'.format(line.strip())
                annotentry['start'] = str(start)
                annotentry['end'] = str(end)
                outline = '\t'.join([ annotentry[k] for k in header ])
                outbuffer.write(outline + '\n')
    return outbuffer


    
    
def run():
    basedir = '/TL/epigenetics2/work/pebert/data/expression'
    rawdir = '/TL/epigenetics2/work/pebert/data/expression/raw'
    gencode_3c = 'gencode.v3c.annotation.GRCh37.gtf'
    outdir = '/TL/epigenetics2/work/pebert/data/expression/processed'
    try:
        gene_dict, transcript_dict = parse_annotation(os.path.join(basedir, gencode_3c))
        print("Parsed {} genes and {} transcripts of type protein_coding".format(len(gene_dict), len(transcript_dict)))
        exprfiles = parse_input_dir(rawdir)
        _ = merge_expr_annot(exprfiles, gene_dict, transcript_dict, outdir)
    except Exception as e:
        sys.stderr.write('Outmost catch: {}\n'.format(e))
        return 1
    return 0

def run2():
    basedir = '/TL/epigenetics2/work/pebert/data/promoters'
    lcg = os.path.join(basedir, '10310Table6_LCG_promoters.csv')
    hcg = os.path.join(basedir, '10310Table6_HCG_promoters.csv')
    biomart_mapping = os.path.join(basedir, 'ensembl_biomart_export.tsv')
    resdict = map_refseq_to_ensembl(lcg, hcg, biomart_mapping)
    print('Creating header')
    for k,v in resdict.items():
        header = list(v.keys())
        break
    outfile = os.path.join(basedir, 'merged_promoter_classes.tsv')
    print('writing to file')
    with open(outfile, 'w') as outf:
        outf.write('#')
        outf.write('\t'.join(header) + '\n')
        for v in resdict.values():
            outstr = '\t'.join([ v[h] for h in header ])
            outf.write(outstr + '\n')
    return 0

def run3():
    basedir = '/TL/epigenetics2/work/pebert/data/promoters'
    lcg = os.path.join(basedir, '10310Table6_LCG_promoters.csv')
    hcg = os.path.join(basedir, '10310Table6_HCG_promoters.csv')
    biomart_mapping = os.path.join(basedir, 'ensembl_biomart_export.tsv')
    promoter_dict = map_refseq_to_ensembl(lcg, hcg, biomart_mapping)
    basedir = '/TL/epigenetics2/work/pebert/projects/validation/featuredatasets'
    bedfiles = os.listdir(basedir)
    header = ''
    for bedf in bedfiles:
        if '85' in bedf or '86' in bedf:
            print('Working on file {}'.format(bedf))
            beddict, header, low_cpg, high_cpg = parse_bedfile(os.path.join(basedir, bedf))
            for k,v in beddict.items():
                if k in promoter_dict:
                    prom = promoter_dict[k]
                    v['promoter_class'] = ('g' + prom['promoter_class'])
                    v['RefSeq_mRNA'] = prom['RefSeq_mRNA']
                    continue
                if float(v['perc_CpG']) <= low_cpg:
                    v['promoter_class'] = 'dLCG'
                    v['RefSeq_mRNA'] = 'NOTAVAILABLE'
                    continue
                if float(v['perc_CpG']) >= high_cpg:
                    v['promoter_class'] = 'dHCG'
                    v['RefSeq_mRNA'] = 'NOTAVAILABLE'
                    continue
                v['promoter_class'] = 'UNDEF'
                v['RefSeq_mRNA'] = 'NOTAVAILABLE'
            new_file = bedf.replace('ext5kb', 'promCl')
            with open(os.path.join(basedir, new_file), 'w') as outfile:
                outfile.write('#')
                outfile.write('\t'.join(header))
                outfile.write('\n')
                for entry in beddict.values():
                    outline = "\t".join([ entry[h] for h in header ])
                    outfile.write(outline)
                    outfile.write('\n')
            print('File processed: {}'.format(new_file))
    return 0


if __name__ == "__main__":
    sys.exit(run3())
