from _BaseFeatures cimport Interval, IntervalFile

cdef class GFeature(Interval):
    cdef public float score
    cdef public str featuretype, method
    cdef public str phase
    cdef dict _attributes
    cdef str _strattributes
    cdef int _attrs_parsed
    cdef str _attribute_delimiter
    cdef str _field_sep

    cdef int parse_line(self)
    cpdef int add_attributes(self,dict d)
    cdef int _parse_attributes(self)

cdef class GTFFile(IntervalFile): pass
cdef class GFFFile(IntervalFile): pass


