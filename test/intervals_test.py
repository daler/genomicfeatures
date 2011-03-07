#!/usr/bin/python
import _genomicfeatures 
import sys
import pysam
import os

this_file = os.path.abspath(__file__)
this_dir = os.path.split(this_file)[0]

c = 0
limit = 100
short = False

if sys.argv[1] == 'gff':
    interval_file = _genomicfeatures.GFFFile(os.path.join(this_dir, 'dm3.gff'))
if sys.argv[1] == 'gtf':
    interval_file = _genomicfeatures.GTFFile(os.path.join(this_dir, 'hg19-genes.gtf'))
if sys.argv[1] == 'bed':
    interval_file = _genomicfeatures.BEDFile(os.path.join(this_dir, 'hg19-genes.bed'))
if sys.argv[1] == 'bed3':
    interval_file = _genomicfeatures.BEDFile(os.path.join(this_dir, bed_fn+'3'))
if sys.argv[1] == 'sam':
    interval_file = _genomicfeatures.SAMFile(os.path.join(this_dir, 'example.sam'))
if sys.argv[1] == 'bam':
    interval_file = _genomicfeatures.BAMFile(os.path.join(this_dir, 'example.bam'))

for feature in interval_file:
    pass

# generic attributes
print feature
print '\tstrand:', feature.strand
print '\tTSS:',feature.tss()
print '\tnfields:', feature.nfields
print '\tlen:', len(feature)
print '\tstart:', feature.start
print '\tstop:', feature.stop
print '\tmidpoint:', feature.midpoint


if sys.argv[1] in ['gtf','gff']:
    print 
    print 'specific to GTF and GFF'
    print '\tattributes:', feature.attributes
    print '\tmethod:', feature.method
    print '\tscore:', feature.score

if sys.argv[1] in ['bam','sam']:
    print
    print 'specific to BAM and SAM'
    print '\talignments:', feature.alignments
    print '\tcigar:', feature.pysam_read.cigar
    print '\ttags:', feature.tags
