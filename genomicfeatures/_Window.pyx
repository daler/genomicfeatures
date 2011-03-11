from collections import deque


cdef class Window(object):
    cdef public object iterable
    cdef public object features
    cdef public int pos
    cdef public int halfwidth
    cdef public int limit
    cdef public int debug
    cdef object buffered_feature, _first_read
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

        If *debug* = 0, then lots of debugging info will be printed.  You
        probably will want *limit* to be set!
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

        self._first_read = self.iterable.next()
        self.left_edge = self._first_read.start
        self.right_edge = self.left_edge + 2*self.halfwidth
        self.chrom = self._first_read.chrom
        self.STARTING = 1

    cdef int within_current_window(self, str chrom, int start):
        """
        if self.features[0].start < feature.start < right_edge, then we're
        good. Otherwise store the feature.
        """
        if self.debug:
            print '\tchecking if %s:%s within %s-%s...' % (chrom,start, self.left_edge, self.right_edge),
        if self.STARTING:
            self.features.append(self._first_read)
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

    cdef int set_window_edges(self) except -1:
        try:
            self.left_edge = self.features[0].start
            self.chrom = self.features[0].chrom
        except IndexError:
            pass
        self.right_edge = self.left_edge + 2*self.halfwidth
        if self.debug:
            print 'new edges: %s:%s-%s' % (self.chrom, self.left_edge, self.right_edge)
        return 0

    def __next__(self):
        """
        Returns deques of all reads whose start positions fall within
        2*self.halfwidth bp of the first read in the window.
        """
         
        # For debugging
        if self.limit > 0:
            if self.counter > self.limit:
                raise StopIteration
        
        # Reset the window edges; not sure if this needs to be done here.
        if len(self.features) > 0:
            self.set_window_edges()

        # Every "next()" call we assume it's not ready, so the loop below will
        # run at least once. 
        #
        # self.READY will be set to 1 when:
        #
        #   * self.accumulate_features() hits a read that won't fit, so the
        #     current window is ready
        #
        #   * self.trim() pops a batch of duplicates off the beginning of the
        #     window, so now it's ready
        #
        self.READY = 0

        while self.READY == 0:
            
            # For stopping early...
            self.counter += 1
            
            # "self.TRIM = 0" means we're not trimming, we're accumulating.  So
            # get as many reads as can fit in the window
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
            # read will fit.  This is an "else" because we don't want to run it
            # till the next call.  self.TRIM will be 1 if:
            #
            #   * self.accumulate_features() finishes running
            #
            #   * self.trim() does a round of trimming, but the buffered read
            #     is still out-of-range even from the newly trimmed window

            elif self.TRIM == 1 and self.READY == 0: 
                self.trim()
                if self.debug:
                    print 'done trimming.  features now:', [i.start for i in self.features]
                if len(self.features) == 0:
                    if self.debug:
                        print 'features is empty after trimming -- appending buffered feature and continuing the while-loop in next();',
                    self.features.append(self.buffered_feature)
                    if self.debug:
                        print 'features now:', [i.start for i in self.features]
                    self.set_window_edges()
                    self.TRIM = 0
                    self.READY = 0

            # If we got through all of that, but still nothing in the window,
            # then we're still not ready.  Head around back to the beginning of
            # the while-loop
            if len(self.features) == 0: 
                if self.debug:
                    print 'still nothing!  moving on to next read'
                self.READY = 0

        return self.features
            
    def __iter__(self):
        return self


