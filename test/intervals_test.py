#!/usr/bin/python
import _genomicfeatures 
import sys
import pysam
import os

this_file = os.path.abspath(__file__)
this_dir = os.path.split(this_file)[0]

gtf_fn = 'hg19-genes.gtf'
bed_fn = 'hg19-genes.bed'
c = 0
limit = 100
short = False
if sys.argv[1] == 'gtf':
    interval_file = _genomicfeatures.GTFFile(os.path.join(this_dir, gtf_fn))
if sys.argv[1] == 'bed':
    interval_file = _genomicfeatures.BEDFile(os.path.join(this_dir, bed_fn))
if sys.argv[1] == 'bed3':
    interval_file = _genomicfeatures.BEDFile(os.path.join(this_dir, bed_fn+'3'))
if sys.argv[1] == 'sam':
    interval_file = _genomicfeatures.SAMFile(os.path.join(this_dir, 'example.sam'))
if sys.argv[1] == 'bam':
    interval_file = _genomicfeatures.BAMFile(os.path.join(this_dir, 'example.bam'))

for feature in interval_file:
    """
    if c > limit:
        break
    c += 1
    """
    #feature.midpoint()
    pass

print feature
print 'strand:', feature.strand
print 'TSS:',feature.tss()
print 'nfields:', feature.nfields
print 'len:', len(feature)
print 'alignments:', feature.alignments
print 'cigar:', feature.pysam_read.cigar
print 'tags:', feature.tags
