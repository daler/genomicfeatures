cdef class CompositeInterval(object):
    # attributes
    cdef public int nfields
    cdef public object fields
    cdef public object featuretypes
    cdef public object features
    
    # methods
    cpdef split(self)

cdef class Interval(object):
    # attributes
    cdef public int start, stop
    cdef public str chrom, strand
    cdef public str _line
    cdef object _split_line
    cdef public int nfields
    cdef public str other_attributes

    # methods
    cpdef int midpoint(self)
    cpdef int tss(self)

cdef class GenericInterval(Interval):
    pass


cdef class IntervalFile(object):
    # attributes
    cdef type _featureclass
    cdef object _handle

    # methods
    cdef int is_invalid(self, str line)
