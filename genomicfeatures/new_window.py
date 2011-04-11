import numpy as np
import genomicfeatures
import os
from collections import deque

## {{{ http://code.activestate.com/recipes/576611/ (r11)
from operator import itemgetter
from heapq import nlargest
from itertools import repeat, ifilter

class Counter(dict):
    '''Dict subclass for counting hashable objects.  Sometimes called a bag
    or multiset.  Elements are stored as dictionary keys and their counts
    are stored as dictionary values.

    >>> Counter('zyzygy')
    Counter({'y': 3, 'z': 2, 'g': 1})

    '''

    def __init__(self, iterable=None, **kwds):
        '''Create a new, empty Counter object.  And if given, count elements
        from an input iterable.  Or, initialize the count from another mapping
        of elements to their counts.

        >>> c = Counter()                           # a new, empty counter
        >>> c = Counter('gallahad')                 # a new counter from an iterable
        >>> c = Counter({'a': 4, 'b': 2})           # a new counter from a mapping
        >>> c = Counter(a=4, b=2)                   # a new counter from keyword args

        '''        
        self.update(iterable, **kwds)

    def __missing__(self, key):
        return 0

    def most_common(self, n=None):
        '''List the n most common elements and their counts from the most
        common to the least.  If n is None, then list all element counts.

        >>> Counter('abracadabra').most_common(3)
        [('a', 5), ('r', 2), ('b', 2)]

        '''        
        if n is None:
            return sorted(self.iteritems(), key=itemgetter(1), reverse=True)
        return nlargest(n, self.iteritems(), key=itemgetter(1))

    def elements(self):
        '''Iterator over elements repeating each as many times as its count.

        >>> c = Counter('ABCABC')
        >>> sorted(c.elements())
        ['A', 'A', 'B', 'B', 'C', 'C']

        If an element's count has been set to zero or is a negative number,
        elements() will ignore it.

        '''
        for elem, count in self.iteritems():
            for _ in repeat(None, count):
                yield elem

    # Override dict methods where the meaning changes for Counter objects.

    @classmethod
    def fromkeys(cls, iterable, v=None):
        raise NotImplementedError(
            'Counter.fromkeys() is undefined.  Use Counter(iterable) instead.')

    def update(self, iterable=None, **kwds):
        '''Like dict.update() but add counts instead of replacing them.

        Source can be an iterable, a dictionary, or another Counter instance.

        >>> c = Counter('which')
        >>> c.update('witch')           # add elements from another iterable
        >>> d = Counter('watch')
        >>> c.update(d)                 # add elements from another counter
        >>> c['h']                      # four 'h' in which, witch, and watch
        4

        '''        
        if iterable is not None:
            if hasattr(iterable, 'iteritems'):
                if self:
                    self_get = self.get
                    for elem, count in iterable.iteritems():
                        self[elem] = self_get(elem, 0) + count
                else:
                    dict.update(self, iterable) # fast path when counter is empty
            else:
                self_get = self.get
                for elem in iterable:
                    self[elem] = self_get(elem, 0) + 1
        if kwds:
            self.update(kwds)

    def copy(self):
        'Like dict.copy() but returns a Counter instance instead of a dict.'
        return Counter(self)

    def __delitem__(self, elem):
        'Like dict.__delitem__() but does not raise KeyError for missing values.'
        if elem in self:
            dict.__delitem__(self, elem)

    def __repr__(self):
        if not self:
            return '%s()' % self.__class__.__name__
        items = ', '.join(map('%r: %r'.__mod__, self.most_common()))
        return '%s({%s})' % (self.__class__.__name__, items)

    # Multiset-style mathematical operations discussed in:
    #       Knuth TAOCP Volume II section 4.6.3 exercise 19
    #       and at http://en.wikipedia.org/wiki/Multiset
    #
    # Outputs guaranteed to only include positive counts.
    #
    # To strip negative and zero counts, add-in an empty counter:
    #       c += Counter()

    def __add__(self, other):
        '''Add counts from two counters.

        >>> Counter('abbb') + Counter('bcc')
        Counter({'b': 4, 'c': 2, 'a': 1})


        '''
        if not isinstance(other, Counter):
            return NotImplemented
        result = Counter()
        for elem in set(self) | set(other):
            newcount = self[elem] + other[elem]
            if newcount > 0:
                result[elem] = newcount
        return result

    def __sub__(self, other):
        ''' Subtract count, but keep only results with positive counts.

        >>> Counter('abbbc') - Counter('bccd')
        Counter({'b': 2, 'a': 1})

        '''
        if not isinstance(other, Counter):
            return NotImplemented
        result = Counter()
        for elem in set(self) | set(other):
            newcount = self[elem] - other[elem]
            if newcount > 0:
                result[elem] = newcount
        return result

    def __or__(self, other):
        '''Union is the maximum of value in either of the input counters.

        >>> Counter('abbb') | Counter('bcc')
        Counter({'b': 3, 'c': 2, 'a': 1})

        '''
        if not isinstance(other, Counter):
            return NotImplemented
        _max = max
        result = Counter()
        for elem in set(self) | set(other):
            newcount = _max(self[elem], other[elem])
            if newcount > 0:
                result[elem] = newcount
        return result

    def __and__(self, other):
        ''' Intersection is the minimum of corresponding counts.

        >>> Counter('abbb') & Counter('bcc')
        Counter({'b': 1})

        '''
        if not isinstance(other, Counter):
            return NotImplemented
        _min = min
        result = Counter()
        if len(self) < len(other):
            self, other = other, self
        for elem in ifilter(self.__contains__, other):
            newcount = _min(self[elem], other[elem])
            if newcount > 0:
                result[elem] = newcount
        return result

class Window(object):
    def __init__(self, iterable, windowsize=100, debug=0):
        """
        Moving window over an *iterable* of features (e.g., BAMFile(bamfn)) of
        size *windowsize*.  Use *debug=1* to see all sorts of output for
        double-checking.

        The resulting Window instance can be iterated over.  Each iteration
        returns a tuple of::

            (center, low_reads, high_reads)

        where *center* is the current center of the window; *low_reads* is a
        deque of reads that includes the center and everything lower than it
        that fits within the window; and *high_reads* is a deque of reads that
        includes everything higher within the window.

        The strategy is to hold one read as the "centered read", which is
        currently in focus.  Reads are checked to see if they fit within the
        window centered on this read.  There is always a buffered read, which
        is last read taken from the iterable.  If the buffered read doesn't fit
        in the window, it remains the buffered read until the current window is
        returned.  It will continue to remain the buffered read (and no more
        reads will be taken from the iterable) until it fits within the current
        window.

        The window is implemented in two parts, a low_reads and a high_reads part.

        The next window's center jumps to the next available read position.
        This will typically be the first item in the high_reads deque.

        """
        self.iterable = iterable
        self.windowsize = windowsize
        self.left_edge = 0
        self.right_edge = 0
        self.debug = debug

        # Here we pull the first thing from the iterable to set up the various
        # attributes
        first_read = self.iterable.next()
        self.chrom = first_read.chrom
        first_start_pos = first_read.start
        self.left_edge = first_start_pos - self.windowsize/2
        self.right_edge = self.left_edge + self.windowsize
        self.center = first_start_pos
        self.buffered_read = first_read
        self.high_reads = deque()
        self.low_reads = deque([self.buffered_read])
        self.START = 1

    def accumulate_reads(self):
        """
        Fill up the window surrounding the currently-centered read.
        """
        if self.debug:
            print 'appending:\n\t',

        while True:


            # Need to short-circuit if starting, cause we've already filled
            # buffered_read
            if self.START:
                self.START = 0
                self.buffered_read = self.iterable.next()
                continue

            if self.buffered_read.chrom != self.chrom:
                if self.debug:
                    print 'new chrom -- %s' % self.buffered_read.chrom
                break

            # While accumulating, the only time low_reads will fill up is if
            # they are duplicates of the currently-centered read
            if self.buffered_read.start == self.center:

                if self.debug:
                    print self.buffered_read.start,

                self.low_reads.append(self.buffered_read)

            # Otherwise, if it's within the window then it's added to
            # high_reads.
            elif self.buffered_read.start < self.right_edge:

                if self.debug:
                    print  self.buffered_read.start,

                self.high_reads.append(self.buffered_read)

            else:
                break

            # The positioning of this is important -- we only get a new
            # buffered read if the last buffered read has been treated --
            # either added to low_reads or high_reads
            self.buffered_read = self.iterable.next()

        if self.debug:
            print

    def trim(self):
        """
        Trims reads off window edges, which is basically just shifting the
        window.
        """

        # If there is nothing in the high reads, then use the current buffered
        # read as the center.
        if len(self.high_reads) == 0:
            self.center = self.buffered_read.start
            self.chrom = self.buffered_read.chrom
            self.left_edge = self.center - self.windowsize/2
            self.right_edge = self.center + self.windowsize/2

        # Otherwise, use the next read in the high_reads deque
        else:
            self.chrom = self.high_reads[0].chrom
            self.center = self.high_reads[0].start
            self.left_edge = self.center - self.windowsize/2
            self.right_edge = self.center + self.windowsize/2

        # Now that the center point has been reset, remove reads from low_reads
        # list that no longer fit in the window
        if self.debug:
            print 'removed:', 
        while True:

            # Must be a better way to do this other than popping it off and
            # then back on again if it's in range, though the appendleft will
            # only happen during one (i.e. the last) time through the loop
            try:
                popped = self.low_reads.popleft()
                if (popped.start < self.left_edge) or (popped.chrom != self.buffered_read.chrom):
                    if self.debug:
                        print popped.start,
                    continue
                else:
                    self.low_reads.appendleft(popped)
                    break

            # If there's nothing left in the low_reads, then stop removing
            except IndexError:
                break

        # Next we remove any additional reads that are duplicates of the
        # centered read and add these to low_reads
        while True: 
            try:
                popped = self.high_reads.popleft()
                if popped.start == self.center:
                    self.low_reads.append(popped)
                else:
                    self.high_reads.appendleft(popped)
                    break
            except IndexError:
                break

        # Run accumulator again to see if we can add the current buffered read
        # and/or any additional reads to the window.
        #self.accumulate_reads()
        if self.debug:
            print 

    def __iter__(self):
        return self

    def next(self):

        if not self.START:
            # This moves the window...
            self.trim()

        # First we accumulate reads
        self.accumulate_reads()

        if self.debug:
            print 'chrom        :', self.chrom
            print 'left         :', self.left_edge
            print 'center       :', self.center
            print 'right        :', self.right_edge
            print 'low contents :', [i.start for i in self.low_reads]
            print 'high contents:', [i.start for i in self.high_reads]
            print 'buffer       :', self.buffered_read.start

        return self.center, self.low_reads, self.high_reads

def score(x, center, windowsize):
    """
    Computes the number of duplicates in the center and the average number
    within the window of size *windowsize*
    """
    c = Counter([i.start for i in x])
    center_count = float(c.pop(center))
    total = sum(c.values())
    num = len(c)
    if num == 0:
        return center_count
    mn = total / (num)
    return center_count / (mn)

if __name__ == "__main__":
    bam_fn = os.path.join(os.path.dirname(__file__), '../timing/window_example.bed')
    #bam_fn = os.path.join(os.path.dirname(__file__), '../timing/example.bam')
    #iterable = genomicfeatures.BAMFile(bam_fn)
    iterable = genomicfeatures.BEDFile(bam_fn)
    w = Window(iterable,windowsize=100, debug=True)
    halfwidth = w.windowsize/2
    windowsize = w.windowsize
    for i in w:

        center = i[0]
        reads = list(i[1])
        reads.extend(i[2])

        print
        print '\t',
        print [i.start for i in reads]
        print '\t',  center, score(reads, center=center, windowsize=windowsize)
        print




