from _BaseFeatures import Interval, IntervalFile

cdef class BEDFeature(Interval):
    
    def __cinit__(self, str line):
        # BED features have at least 3 fields, but this can change depending on
        # the parsed line
        self.nfields = 3
        self._line = line

    def __init__(self, str line):
        self.parse_line()
        
    cdef int parse_line(self) except -1:
        """
        Even though we don't *need* to return an int, we do here so that
        exceptions will be passed (methods returning void can't pass exceptions
        to caller in Cython).

        That's also why this is its own function instead of being inside
        __init__...
        """
        L = self._line.strip().split('\t')
        self.nfields = len(L)
        
        if self.nfields == 3:
            self.chrom = L[0]
            self.start = int(L[1])
            self.stop  = int(L[2])
            return 0

        # Optimization note: At first I used try:except blocks and caught
        # IndexError exceptions, but this was really slow.  Then used a chain
        # of "if len(L) > N" statements; but this turned out to be slightly
        # slower than the following and also not as readable. 
        
        if self.nfields == 4:
            self.chrom = L[0]
            self.start = int(L[1])
            self.stop  = int(L[2])
            self.name = L[3]
            return 0

        if self.nfields == 5:
            self.chrom = L[0]
            self.start = int(L[1])
            self.stop  = int(L[2])
            self.name = L[3]
            self.score = float(L[4])
            return 0

        if self.nfields == 6:
            self.chrom = L[0]
            self.start = int(L[1])
            self.stop  = int(L[2])
            self.name = L[3]
            self.score = float(L[4])
            self.strand = L[5]
            return 0

        if self.nfields == 8:
            self.chrom = L[0]
            self.start = int(L[1])
            self.stop  = int(L[2])
            self.name = L[3]
            self.score = float(L[4])
            self.strand = L[5]
            self.thick_start = int(L[6])
            self.thick_end = int(L[7])
            return 0

        if self.nfields == 9:
            self.chrom = L[0]
            self.start = int(L[1])
            self.stop  = int(L[2])
            self.name = L[3]
            self.score = float(L[4])
            self.strand = L[5]
            self.thick_start = int(L[6])
            self.thick_end = int(L[7])
            self.item_rgb = L[8]
            return 0

        if self.nfields == 12:
            self.chrom = L[0]
            self.start = int(L[1])
            self.stop  = int(L[2])
            self.name = L[3]
            self.score = float(L[4])
            self.strand = L[5]
            self.thick_start = int(L[6])
            self.thick_end = int(L[7])
            self.item_rgb = L[8]
            self.block_count = int(L[9])
            self.block_sizes = L[10]
            self.block_starts = L[11]
            return 0

    def __str__(self):
        s = ''
        s += '<%s %s:%s-%s(%s)>' % (self.__class__.__name__, self.chrom,self.start,self.stop,self.strand)
        return s

    def __repr__(self):
        return self._line



cdef class BEDFile(IntervalFile):

    def __cinit__(self,*args):
        self._featureclass = BEDFeature
    
    def __init__(self, str fn):
        self._handle = open(fn)
    
    cdef int is_invalid(self,str line):
        if line[:5] == 'track':
            return 0
        if line[:7] == 'browser':
            return 0
        else:
            return 1

