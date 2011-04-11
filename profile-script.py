import pstats, cProfile

#import pyximport
#pyximport.install()

import genomicfeatures

cProfile.runctx('genomicfeatures.test_window()', globals(), locals(), 'Profile.prof')
s = pstats.Stats('Profile.prof')
s.strip_dirs().sort_stats('time').print_stats()

