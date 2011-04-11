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
debug = 0
limit = 100000
if 1:
    #window = genomicfeatures.Window(interval_file, debug=debug, limit=35, halfwidth=72)
    windowsize = 72
    window = genomicfeatures.Window(interval_file, debug=debug, windowsize=windowsize)
    c = 0
    for w in window:
        center, low_reads, high_reads = w

        c += 1
        if limit is not None:
            if c > limit:
                break

        score = genomicfeatures.dups_score(low_reads, high_reads, center, windowsize=windowsize)
        #center = w[0].start + window.halfwidth
        if debug:
            print '  %s =========>output:' % c, center, score, 
            print [i.start for i in reads]
        fout.write('\t'.join([ low_reads[0].chrom, str(center), str(center+1), str(score)])+'\n')
        fout.flush()
    fout.close()

