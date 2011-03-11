#!/usr/bin/python
import _genomicfeatures 
import sys
import os

this_file = os.path.abspath(__file__)
this_dir = os.path.split(this_file)[0]

if sys.argv[1] == 'gtf':
    interval_file = _genomicfeatures.GTFFile(os.path.join(this_dir, gtf_fn))
if sys.argv[1] == 'bed':
    interval_file = _genomicfeatures.BEDFile(os.path.join(this_dir, bed_fn))
if sys.argv[1] == 'bed3':
    interval_file = _genomicfeatures.BEDFile(os.path.join(this_dir, bed_fn+'3'))
if sys.argv[1] == 'sam':
    interval_file = _genomicfeatures.SAMFile(os.path.join(this_dir, 'example.sam'))
if sys.argv[1] == 'bam':
    #fn = '/DATA/analysis/suhw-kc-ripseq/suhw1-ip/suhw1-ip.tophat.filtered.bam'
    #fn = 'example.bam'
    fn = '/DATA/analysis/35nt/35nt-1-ip/35nt-1-ip.tophat.filtered.bam'
    
    interval_file = _genomicfeatures.BAMFile(os.path.join(this_dir, fn))


fout = open('test.bedgraph','w')
fout.write('track type=bedGraph name=score_test\n')
debug = 1
if 1:
    window = _genomicfeatures.Window(interval_file, debug=debug, limit=35, halfwidth=72)
    c = 0
    for w in window:
        c += 1
        score = _genomicfeatures.dups_score_sum(w, window.halfwidth)
        center = w[0].start + window.halfwidth
        if debug:
            print '  %s =========>output:' % c, center, score, 
            print [i.start for i in w]
        fout.write('\t'.join([ w[0].chrom, str(center), str(center+1), str(score)])+'\n')
        fout.flush()
    fout.close()

