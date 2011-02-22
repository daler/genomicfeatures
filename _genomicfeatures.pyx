#!/usr/bin/python
# encoding: utf-8
 
# To build for testing:
# python setup.py build_ext --inplace

# Check your optimization in HTML, right from vim....
# setlocal makeprg=cython\ -a\ %\ &&\ google-chrome\ '%<.html' 
# :make

# DONE: separate BEDFeature subclasses for each field length (3-6, 8, 9, 12)
# DONE: BEDFile should check first valid line of file and then dispatch  
#       the right feature type . . . or is it better to have different parser
#       methods instead of entirely new classes?
# DONE: SAMFile, SAMFeature -- wrap pysam's AlignedRead class
# DONE: BAMFile, BAMFeature
# TODO: GFFFile, GFFFeature
# TODO: VCFFile, VCFFeature
# TODO: composite features
# TODO: format auto-detect

import pysam


cdef class AutoDetect(object):
    cdef public str line

    def __init__(self,line):
        raise NotImplementedError, 'this is just a placeholder for now...'
        self.line = line
        
    
    def inspect_fields(self):
        """
        identify the number fields and the string fields.

        Also look for things like cigar strings, seq and qual fields with same
        length, and other format-specific field properties
        """

cdef class Interval(object):
    """
    Generic interval class.  Has start, stop, and strand as well as a couple
    generic methods
    """
    cdef public int start, stop
    cdef public str chrom, strand
    cdef str _line
    cdef object _split_line
    cdef public int nfields

    cpdef int midpoint(self):
        return self.start + (self.stop-self.start)/2
    
    cpdef int tss(self) except -1:
        if self.strand == '+':
            return self.start
        if self.strand == '-':
            return self.stop
        else:
            raise ValueError, 'Undefined strand "%s"' % self.strand
            return -1
    
    def __len__(self):
        return self.stop - self.start

    def __repr__(self):
        return self._line

    def __str__(self):
        return '<%s %s:%s-%s(%s)>' % (self.__class__.__name__, self.chrom, self.start, self.stop, self.strand)

cdef class GenericInterval(Interval):
    def __init__(self, chrom,start,stop,strand):
        self.chrom = chrom
        self.start = start
        self.stop = stop
        self.strand = strand
        self._line = '%s' % ((chrom,start,stop,strand),)

cdef class SAMFeature(Interval):
    # wrapper around pysam AlignedRead object
    cdef public object pysam_read, pysam_samfile


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

cdef class BEDFeature(Interval):
    cdef public str name, item_rgb
    cdef public float score
    cdef public int thick_start, thick_end, block_count
    cdef public str block_sizes, block_starts
    
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


cdef class GTFFeature(Interval):

    # extra attributes that GTFFeatures have over plain ol' Intervals
    cdef public float score
    cdef public str featuretype, method
    cdef public str phase
    cdef dict _attributes
    cdef str _strattributes
    cdef int _attrs_parsed
    
    def __cinit__(self, line):
        self._line = line
        self._attrs_parsed = 0
        
        # GTFs always have exactly 8 fields
        self.nfields = 8

    def __init__(self, str line):
        self.parse_line()

    cdef void parse_line(self):
        """
        Simply split the line and put things where they belong, converting as
        necessary
        """
        cdef object L
        L = self._line.split('\t')
        self.chrom = intern(L[0])
        self.start = int(L[3])
        self.stop = int(L[4])
        self.score = float(L[5])
        self.strand = L[6]
        self.phase = L[7]
        self.method = L[1]
        self.featuretype = L[2]

        # here we're saving just the string of attributes which will be parsed
        # only when asked for
        self._strattributes = L[8]

    # using the property mechanism, we can postpone parsing the attributes
    # until we actually need them . . .
    property attributes:
        def __get__(self):
            if self._attrs_parsed == 0:
                self._parse_attributes()
            return self._attributes
    
    cpdef add_attributes(self,dict d):
        """
        Add dict *d* of field:values to attributes
        """
        self._attributes.update(d)

    cdef void _parse_attributes(self):
        """
        Parse the attributes stored in self._strattributes and store 'em in a
        dictionary, self._attributes.
        """
        cdef dict attrs
        cdef str field, value
        cdef list items
        attrs = {}
        items = self._strattributes.strip().split(';')
        for item in items:
            if len(item) == 0:
                continue
            field, value = item.strip().split()
            value = value.replace('"','')
            attrs[field] = value
        self._attributes = attrs 

        # This lets the getter know that we've already parsed and that from now
        # on it should just look at self._attributes instead of parsing
        # self._strattributes again.
        self._attrs_parsed = 1


cdef class IntervalFile(object):
    cdef type _featureclass
    cdef object _handle

    cdef int is_invalid(self, str line):
        return 1

    def __next__(self):
        line = self._handle.next() 
        while self.is_invalid(line) == 0: 
            line = self._handle.next()
        return self._featureclass(line)
        
    def __iter__(self):
        return self


cdef class SAMFile(object):
    cdef type _featureclass
    cdef object _handle

    def __cinit__(self,*args):
        self._featureclass = SAMFeature
        
    def __init__(self, str fn):
        self._handle = pysam.Samfile(fn)
    
    def __iter__(self):
        return self

    def __next__(self):
        return self._featureclass(self._handle.next(), self._handle)

cdef class BAMFile(SAMFile):
    def __init__(self, str fn):
        self._handle = pysam.Samfile(fn,'rb')

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

