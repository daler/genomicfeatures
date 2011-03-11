
cdef class IntervalFile(object):

    cdef int is_invalid(self, str line):
        return 1

    def __next__(self):
        line = self._handle.next() 
        while True:
            valid = self.is_invalid(line)
            if valid == 0: 
                line = self._handle.next()
            if valid == 1:
                break
            if valid == -1:
                raise StopIteration

        return self._featureclass(line)
        
    def __iter__(self):
        return self



cdef class CompositeInterval(object):

    def __init__(self, this, other):
        assert isinstance(this, Interval)
        assert isinstance(other, Interval)
        self.featuretypes = (this.__class__, other.__class__)
        self.fields = (this.nfields, other.nfields)
        self.features = (this, other)
        self.nfields = sum(self.fields) 

    cpdef split(self):
        """

        """

cdef class Interval(object):
    """
    Generic interval class.  Has start, stop, and strand as well as a couple
    generic methods
    """

    cpdef int midpoint(self):
        return self.start + (self.stop-self.start)/2
    
    cpdef int tss(self) except -1:
        if self.strand == '+':
            return self.start
        if self.strand == '-':
            return self.stop
        else:
            raise ValueError, 'TSS not implemented for undefined strand "%s"' % self.strand
            return -1
    
    def __len__(self):
        return self.stop - self.start

    def __repr__(self):
        return self._line

    def __str__(self):
        return '<%s %s:%s-%s(%s)>' % (self.__class__.__name__, self.chrom, self.start, self.stop, self.strand)

cdef class GenericInterval(Interval):
    def __init__(self, chrom,start,stop,strand,other_attributes=""):
        self.chrom = chrom
        self.start = start
        self.stop = stop
        self.strand = strand
        self.other_attributes = other_attributes
        self._line = '%s' % ((chrom,start,stop,strand,other_attributes),)
