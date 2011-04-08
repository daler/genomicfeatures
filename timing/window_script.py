#!/usr/bin/python
import genomicfeatures 
import sys
import os

this_file = os.path.abspath(__file__)
this_dir = os.path.split(this_file)[0]

fn = 'example.bam'
interval_file = genomicfeatures.BAMFile(os.path.join(this_dir, fn))


fout = open('test.bedgraph','w')
fout.write('track type=bedGraph name=score_test\n')
debug = 1
if 1:
    window = genomicfeatures.Window(interval_file, debug=debug, limit=35, halfwidth=72)
    c = 0
    for w in window:
        c += 1
        score = genomicfeatures.dups_score_sum(w, window.halfwidth)
        center = w[0].start + window.halfwidth
        if debug:
            print '  %s =========>output:' % c, center, score, 
            print [i.start for i in w]
        fout.write('\t'.join([ w[0].chrom, str(center), str(center+1), str(score)])+'\n')
        fout.flush()
    fout.close()

