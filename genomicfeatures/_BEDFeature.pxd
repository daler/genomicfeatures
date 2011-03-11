from _BaseFeatures cimport Interval, IntervalFile

cdef class BEDFeature(Interval):
    cdef public str name, item_rgb
    cdef public float score
    cdef public int thick_start, thick_end, block_count
    cdef public str block_sizes, block_starts
    
    cdef int parse_line(self)

cdef class BEDFile(IntervalFile):
    pass
