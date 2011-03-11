from _BaseFeatures cimport Interval
cdef class SAMFeature(Interval):
    cdef public object pysam_read, pysam_samfile

cdef class SAMFile(object):
    cdef type _featureclass
    cdef object _handle
    cdef public str fn

    cdef int count(self)

cdef class BAMFile(SAMFile): pass
