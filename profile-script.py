import pstats, cProfile

import pyximport
pyximport.install()

import _genomicfeatures

cProfile.runctx('_genomicfeatures._test_windows()', globals(), locals(), 'Profile.prof')
s = pstats.Stats('Profile.prof')
s.strip_dirs().sort_stats('time').print_stats()

