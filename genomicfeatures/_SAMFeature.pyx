from _BaseFeatures cimport Interval, GenericInterval
import pysam
       
cdef class SAMFeature(Interval):

    def __cinit__(self, *args):
        self.nfields = 11

    def __init__(self, pysam_read, pysam_samfile):
        self.pysam_samfile = pysam_samfile
        self.pysam_read = pysam_read

        if self.pysam_read.flag & 0x10:
            self.strand = '-'
        else:
            self.strand = '+'

        self.start = self.pysam_read.pos
        self.stop = self.pysam_read.aend
        self.chrom = self.pysam_samfile.getrname(self.pysam_read.tid)

    property nfields:
        def __get__(self):
            return 11 + len(self.pysam_read.tags)
    
    property alignments:
        def __get__(self):
            alignments = []
            cdef int current_start = self.start
            cdef int current_stop = self.start
            for element in self.pysam_read.cigar:
                operation, bp = element
                if operation not in [0,3,7,8]:
                    lookup = {0: 'M alignment match (can be a sequence match or mismatch)',
                              1: 'I insertion to the reference',
                              2: 'D deletion from the reference',
                              3: 'N skipped region from the reference',
                              4: 'S soft clipping (clipped sequences present in SEQ)',
                              5: 'H hard clipping (clipped sequences NOT present in SEQ)',
                              6: 'P padding (silent deletion from padded reference)',
                              7: '= sequence match',
                              8: 'X sequence mismatch',}
                    raise ValueError, 'CIGAR operation "%s" not yet supported' % lookup[operation]
                if operation in [0,7,8]:
                    current_stop += bp
                    alignments.append(GenericInterval(self.chrom, current_start, current_stop, self.strand))
                    current_start = current_stop+1
                if operation == 3:
                    current_start += bp
            return alignments
                
    property tags:
        def __get__(self):
            return self.pysam_read.tags

    def __repr__(self):
        return str(self.pysam_read)


cdef class SAMFile(object):

    def __cinit__(self,fn):
        self._featureclass = SAMFeature
        self.fn = fn
        
    def __init__(self, str fn):
        self._handle = pysam.Samfile(fn)
    
    def __iter__(self):
        return self

    def __next__(self):
        return self._featureclass(self._handle.next(), self._handle)

    cdef int count(self):
        return int(pysam.view(self.fn, '-c')[0])

cdef class BAMFile(SAMFile):
    def __init__(self, str fn):
        self._handle = pysam.Samfile(fn,'rb')

