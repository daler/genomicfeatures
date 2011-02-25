# cython: profile=True

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
from collections import deque
import numpy as np
cimport numpy as np


cpdef float dups_score(object x, int halfwidth, float scalar=1) except -1:
    """
    Returns the score at the centerpoint for a window of features, *x*.  *x*
    should be a list or deque.

    *halfwidth* is what the halfwidth of the window should be.

    *scalar* is what each score should be mutliplied by (e.g., a scale factor
    for reads per million mapped)

    Score is calculated by taking the sum of reads in the center and iding
    by average number of of non-zero and non-center duplicates.

    Window size is assumed to be (x[0].start + halfwidth).
    
    The window [0, 0, 0, 4, 4, 5, 5, 5, 5, 5, 5, 8] with a halfwidth=5 would
    be counted like this::

                  _
                  _
                  _
        _         _
        _       _ _
        _       _ _     _
        0 1 2 3 4 5 6 7 8 9
        ^         ^
        |         |
        start    center  

        number of center duplicates = 6

        avg nonzero duplicates = (3 reads at pos 0) + (2 reads at pos 4) + (1 read at pos 8)  
                                 ___________________________________________________________  = 2.0
                                                  3 positions with >0 reads
                                                  
        score at pos 5 = 6 / 2.0 = 3.0


        NEXT WINDOW:
          _
          _
          _
          _
        _ _
        _ _     _
        4 5 6 7 8 9 10 11 12 13
                  ^
                  |
                center

       center dups = 0
       score at center = 0

    """
    cdef np.ndarray[np.int_t, ndim=1] starts = np.array([j.start for j in x])
    
    # subtract the left edge, so now the window goes from 0 to 2*halfwidth
    # instead of being in genomic coords
    starts -= starts[0]

    # "The output, b[i], represents the number of times that i is found in
    # `x`."
    #
    # So, for example, 
    #
    # >>> bincount([0, 0, 0, 4, 4, 5, 5, 5, 5, 5, 5, 8])
    # 
    # index: 0  1  2  3  4  5  6  7  8
    #        |  |  |  |  |  |  |  |  |
    # array([3, 0, 0, 0, 2, 6, 0, 0, 1])
    #        ^                       ^
    #    '0' was found 3 times       '8' was found once
    
    cdef np.ndarray[np.int_t, ndim=1] dup_counts = np.bincount(starts)

    # if we can't even get a duplicate count on the centerpoint, then the score
    # is defined to be zero.
    cdef float dc_len 
    dc_len = float(len(dup_counts))
    if dc_len <= halfwidth:
        return 0 

    # Imagine for the example above that halfwidth is 6 
    # There are (2*halfwidth - len(dup_counts) ) zeros that could be appended on the end.
    #
    # dup_counts = [3, 0, 0, 0, 2, 6, 0, 0, 1]

    # What it would look like if we took the time to actually pad it:
    # padded     = [3, 0, 0, 0, 2, 6, 0, 0, 1, 0, 0, 0]
    
    # But instead let's just get how many extra zeros there should be
    padding = 2*halfwidth - dc_len
    
    # To avoid divide-by-zero, let's add 1 to everything.  A side effect is
    # that the lower nonzero bound on a score will be 1 / (2*halfwidth), for
    # the case of a single read all by its lonesome, with no other reads in the
    # window

    # add 1 to everything we have a value for so far
    dup_counts += 1

    # number of actual center duplicates plus 1
    center_dups = dup_counts[halfwidth] 
    
    # the padding is added cause it's like we're adding padding * 1 reads to
    # the window total.
    non_center_dups = dup_counts.sum() + padding - center_dups

    avg_dups = non_center_dups / (dc_len - 1) 
    score = center_dups / avg_dups * scalar
    return score
    

cdef class Window(object):
    cdef public object iterable
    cdef public object features
    cdef public int pos
    cdef public int halfwidth
    cdef public int limit
    cdef public int debug
    cdef object buffered_feature
    cdef int counter
    cdef int TRIM
    cdef int left_edge
    cdef int right_edge
    cdef str chrom
    cdef int READY
    cdef int STARTING

    def __init__(self, iterable, halfwidth=100, limit=10, debug=0):
        """
        *iterable* is a an IntervalFile subclass instance.

        *halfwidth* is 0.5 * how big you want the window to be.
        """
        
        self.debug = debug
        self.READY = 0
        self.TRIM = 0
        self.limit = limit
        self.iterable = iterable
        self.features = deque()
        self.pos = 0
        self.halfwidth = halfwidth
        self.counter = 0

        first_read = self.iterable.next()
        self.left_edge = first_read.start
        self.right_edge = self.left_edge + 2*self.halfwidth
        self.chrom = first_read.chrom
        self.STARTING = 1

    cdef int within_current_window(self, str chrom, int start):
        """
        if self.features[0].start < feature.start < right_edge, then we're
        good. Otherwise store the feature.
        """
        if self.debug:
            print '\tchecking if %s:%s within %s-%s...' % (chrom,start, self.left_edge, self.right_edge),
        if self.STARTING:
            self.STARTING = 0
            if self.debug:
                print '(yes, starting)',
            return 1
        if start <= self.right_edge:
            if chrom == self.chrom:
                if self.debug:
                    print '(yes)',
                return 1
        if self.debug:
            print '(no)',
        return 0

    cdef int accumulate_features(self) except -1:
        """
        Grabs new features as long as they're within range of the window
        """
        while True:
            self.buffered_feature = self.iterable.next()
             
            if self.within_current_window(self.buffered_feature.chrom, self.buffered_feature.start) == 1:
                self.features.append(self.buffered_feature)
                if self.debug:
                    print 'added buffered read to window'
            
            # if not in range then set the TRIM flag, which __next__ will see
            else:
                if self.debug:
                    print 'buffered read not in range, saving for later'
                self.READY = 1
                self.TRIM = 1
                break
        if len(self.features) == 0:
            self.READY = 0 
            self.TRIM = 0

    cdef int trim(self) except -1:
        """
        Trims duplicates off of self.features.  Duplicates are defined as
        having identical chrom and start positions.  The currently buffered
        read will be appended to the new window, if it fits.
        """
        if len(self.features) == 0:
            return 0
        first_item = self.features[0]

        while len(self.features) > 0:
            if self.debug:
                print 'inside trimming loop...',
            # if they're the same, then popleft 
            if (first_item.start == self.features[0].start) and (first_item.chrom == self.features[0].chrom):
                _ = self.features.popleft()
                if self.debug:
                    print 'popping off leftmost read', _.start
                continue

            # if not, then break!
            else:
                if self.debug:
                    print self.features[0].start, 'not a duplicate of', first_item.start, 'so done trimming!'
                break
        
        # check to see if the buffered read will fit
        if len(self.features) > 0:
            self.set_window_edges()
        
        if self.within_current_window(self.buffered_feature.chrom, self.buffered_feature.start):

            if self.debug:
                print 'appending buffered read within trim()'
            self.features.append(self.buffered_feature)
            self.TRIM = 0
            self.READY = 0
        
        else:
            if self.debug:
                print 'buffered feature won\'t fit new window, will need trimming next time'
            # buffered read won't fit, so we return the trimmed window and
            # prepare for trimming next time
            self.TRIM = 1
            self.READY = 1
        
        return 0

        
    cdef int set_window_edges(self):

        self.left_edge = self.features[0].start
        self.right_edge = self.left_edge + 2*self.halfwidth
        self.chrom = self.features[0].chrom
        if self.debug:
            print 'new edges: %s:%s-%s' % (self.chrom, self.left_edge, self.right_edge)
        return 0

    def __next__(self):
        """
        Returns deques of all reads whose start positions fall within
        2*self.halfwidth bp of the first read in the window.
        """
         
        if self.limit > 0:
            if self.counter > self.limit:
                raise StopIteration
        if len(self.features) > 0:
            self.set_window_edges()
        self.READY = 0
        while self.READY == 0:
            
            # for stopping early...
            self.counter += 1
            
            # "self.TRIM = 0" means get as many reads as can fit in the window
            if self.TRIM == 0:    
                if self.debug:
                    print 'entering accumulation'
                self.accumulate_features()
                if self.debug:
                    print 'done accumulation, features now:', [i.start for i in self.features]
                if len(self.features) > 0:
                    self.set_window_edges()
                    break

            # "self.TRIM = 1" means pop off duplicates and see if the buffered
            # read will fit.

            # NOTE: might need to test for !READY here.
            elif self.TRIM == 1 and self.READY == 0: 

                self.trim()
                if len(self.features) == 0:
                    if self.debug:
                        print 'features is empty -- adding buffered read in trim()'
                    self.features.append(self.buffered_feature)

            if len(self.features) == 0:
                self.READY = 0

        return self.features
            
    def __iter__(self):
        return self


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
    cdef public str fn

    def __cinit__(self,fn):
        self._featureclass = SAMFeature
        self.fn = fn
        
    def __init__(self, str fn):
        self._handle = pysam.Samfile(fn)
    
    def __iter__(self):
        return self

    def __next__(self):
        return self._featureclass(self._handle.next(), self._handle)

    def count(self):
        return int(pysam.view(self.fn, '-c')[0])

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

def _test_windows():
    fn = 'test/example.bam'
    window = Window(BAMFile(fn), debug=0, limit=-1, halfwidth=1000)
    c = 0

    for w in window:
        c += 1
        if c > 50000:
            break
        pass
        #score = dups_score(w, window.halfwidth)

