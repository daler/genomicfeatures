
from _Scores import dups_score as score
import genomicfeatures
import os
def test_window():
    #bam_fn = os.path.join(os.path.dirname(__file__), '../timing/window_example.bed')
    bam_fn = os.path.join(os.path.dirname(__file__), '../timing/example.bam')
    iterable = genomicfeatures.BAMFile(bam_fn)
    #iterable = genomicfeatures.BEDFile(bam_fn)
    w = genomicfeatures.Window(iterable, windowsize=100, debug=0)
    windowsize = w.windowsize
    c = 0
    for i in w:
        c += 1
        #if c > 50000:
        #    break
        center, low_reads, high_reads = i
        genomicfeatures.dups_score(low_reads, high_reads, center, windowsize)

