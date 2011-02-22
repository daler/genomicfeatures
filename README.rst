genomicfeatures
===============
Classes and iterators for files of genomic features like BED, GTF, GFF,
SAM, BAM, and VCF.

The classes are written in Cython to try to gain some speed, and any
instance attributes that aren't super fast to create (e.g., parsing GTF
attributes into string, floats, etc as needed) are done only when asked for
("lazily").

SAM and BAM support require pysam to be installed.  The ``SAMFeature``
class is a loose wrapper around the ``pysam.AlignedRead`` class in order to
make it look more like the BED and GFF features in ``genomicfeatures`` and
take advantage of any new functionality added to to the ``Interval`` base
class.
