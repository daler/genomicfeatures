import numpy as np
cimport numpy as np

cpdef float dups_score_sum(object x, int halfwidth, float scalar=1) except -1:
    """
    Returns the score at the centerpoint for a window of features, *x*.  *x*
    should be a list or deque.

    In contrast to dups_score(), which computes the *ratio* of center to
    average window, this function takes the *difference*.
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
    score = (center_dups - avg_dups) * scalar
    return score

cpdef float dups_score(object x, int halfwidth, float scalar=1) except -1:
    """
    Returns the score at the centerpoint for a window of features, *x*.  *x*
    should be a list or deque.

    *halfwidth* is what the halfwidth of the window should be.

    *scalar* is what each score should be mutliplied by (e.g., a scale factor
    for reads per million mapped)

    Score is calculated by taking the sum of reads in the center and dividing
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
    #
    # This is the number of duplicates at each location.
    

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
    # Here's what dup_counts might look like if there are only reads 9
    # positions into the window:
    #
    # dup_counts = [3, 0, 0, 0, 2, 6, 0, 0, 1]
    #
    # What it would look like if we took the time to actually pad it:
    # padded     = [3, 0, 0, 0, 2, 6, 0, 0, 1, 0, 0, 0]
    #
    # But instead let's just get how many extra zeros there should be
    padding = 2*halfwidth - dc_len
    
    # To avoid divide-by-zero, let's add 1 to everything.  A side effect is
    # that the lower nonzero bound on a score will be 1 / (2*halfwidth), for
    # the case of a single read all by its lonesome, with no other reads in the
    # window

    # Originally, I was adding 1 to everything we have a value for so far:
    #
    # dup_counts += 1
    #
    # but really only the window average should be +1

    # number of actual center duplicates
    center_dups = dup_counts[halfwidth]
    
    # the padding is added cause it's like we're adding padding * 1 reads to
    # the window total.

    # set the center to be zero
    dup_counts[halfwidth] = 0

    avg_dups = dup_counts.mean() 
    score = center_dups / (avg_dups+1.0) * scalar
    return score
    

