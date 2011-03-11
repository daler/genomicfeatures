from _BaseFeatures cimport Interval, IntervalFile

cdef class GFeature(Interval):
    """
    Base class for GTF and GFF files, which pretty much only differ in their
    attributes.
    """

    cdef int parse_line(self) except -1:
        """
        Simply split the line and put things where they belong, converting as
        necessary
        """
        cdef object L
        L = self._line.split('\t')
        self.chrom = intern(L[0])
        self.start = int(L[3])
        self.stop = int(L[4])
        try:
            self.score = float(L[5])
        except ValueError:
            self.score = 0
        self.strand = L[6]
        self.phase = L[7]
        self.method = L[1]
        self.featuretype = L[2]

        # here we're saving just the string of attributes which will be parsed
        # only when asked for
        self._strattributes = L[8]
        
        # dictionary, where attributes will eventually go
        self._attributes = {}

    # Using the property mechanism, we can postpone parsing the attributes
    # until we actually need them . . .
    property attributes:
        def __get__(self):
            if self._attrs_parsed == 0:
                self._parse_attributes()
            return self._attributes

        def __set__(self, value):
            if not isinstance(value, dict):
                raise ValueError, 'attributes must be a dictionary'
            self._attributes = value
    
    cpdef int add_attributes(self,dict d):
        """
        Add dict *d* of field:values to attributes
        """
        str_to_add = self._attribute_delimiter.join(['%s%s%s'%(key, self._field_sep, value) for key,value in d.items()])
        if self._strattributes[-1] != self._attribute_delimiter:
            self._strattributes += self._attribute_delimiter
        self._strattributes += str_to_add
        
        # tell the attributes method that _strattributes has changed and needs
        # to be re-parsed next time attributes are asked for
        self._attrs_parsed = 0
        return 0

    cdef int _parse_attributes(self) except -1:
        """
        Parse the attributes stored in self._strattributes and store 'em in a
        dictionary, self._attributes.
        """
        cdef dict attrs
        cdef str field, value
        cdef list items
        attrs = {}
        try:
            items = self._strattributes.strip().split(';')
        except AttributeError:
            print self._strattributes
        for item in items:
            if len(item) == 0:
                continue
            field, value = item.strip().split(self._field_sep)
            value = value.replace('"','')
            attrs[field] = value
        self._attributes = attrs 

        # This lets the getter know that we've already parsed; from now on it
        # should just look at self._attributes instead of parsing
        # self._strattributes again.
        self._attrs_parsed = 1


cdef class GFFFeature(GFeature):
    def __init__(self, str line):
        self._line = line
        self._attrs_parsed = 0
        self._attribute_delimiter = ';'
        self._field_sep = '='
        self.nfields = 8
        self.parse_line()

cdef class GTFFeature(GFeature):
    def __init__(self, str line):
        self._line = line
        self._attrs_parsed = 0
        self._attribute_delimiter = ';'
        self._field_sep = ' '
        self.nfields = 8
        self.parse_line()



cdef class GTFFile(IntervalFile):

    def __cinit__(self,*args):
        self._featureclass = GTFFeature
    
    def __init__(self, str fn):
        self._handle = open(fn)
    
    cdef int is_invalid(self,str line):
        if line[0] == '#':
            return 0
        else:
            return 1

cdef class GFFFile(IntervalFile):

    def __cinit__(self,*args):
        self._featureclass = GFFFeature
    
    def __init__(self, str fn):
        self._handle = open(fn)
    
    cdef int is_invalid(self,str line):
        # complete stop if we hit the FASTA section of a GFF file
        if line[0] == '>':
            return -1
        if line[0] == '#':
            return 0
        else:
            return 1


