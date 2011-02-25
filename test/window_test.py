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
    fn = '/DATA/analysis/shep-ripseq/shep3_hs-ip/shep3_hs-ip.bam'
    #fn 'example.bam'
    interval_file = _genomicfeatures.BAMFile(os.path.join(this_dir, fn))


fout = open('test.bedgraph','w')
fout.write('track type=bedGraph name=score_test\n')
window = _genomicfeatures.Window(interval_file, debug=0, limit=-1, halfwidth=72)
for w in window:
    score = _genomicfeatures.dups_score(w, window.halfwidth)
    center = w[0].start + window.halfwidth
    fout.write('\t'.join([ w[0].chrom, str(center), str(center+1), str(score)])+'\n')
    #fout.flush()
fout.close()

    
