from _BaseFeatures import GenericInterval
from _GFeatures import GFFFeature, GTFFeature, GFFFile, GTFFile
from _BEDFeature import BEDFeature, BEDFile
from _SAMFeature import SAMFeature, SAMFile, BAMFile
from _Window import Window
from _Scores import dups_score, dups_score_sum


def inspect_fields(line):
    """
    Given a *line*, try to figure out what kind of feature it is, then return
    that feature.

    If it can't, it will return a GenericInterval by making some assumptions
    about chrom, start, stop, and strand
    """
    fields = line.split('\t')
    
    nfields = len(fields)
    # only BED files (out of currently supported formats) can have fewer
    # than 8 fields
    if nfields < 8:
        return BEDFeature(line)

    if nfields == 11:
        return SAMFeature(line)

    # could be GFF or GTF.  Check the format of the last field...
    if nfields == 8:
        attributes = fields[-1]
        if len(attributes) == 0:
            raise ValueError, 'ambiguous feature type -- GFF or GTF; last field empty'

        # does it have about as many '='s as ';'s?
        num_eq = attributes.count('=')
        num_semi = attributes.count(';')
        if (num_eq == num_semi) or (num_eq == num_semi+1):
            return GFFFeature(line)
        elif ('gene_id' in attributes) and ('transcript_id' in attributes):
            return GTFFeature(line)

    # If you got here, you can't figure out what it is.  Make some
    # assumptions, like chrom being the first field that can't be turned
    # into an int, or strand is the first thing that's  either a '+','-',
    # or '.'
    start = None
    stop = None
    chrom = None
    strand = None
    for i,field in enumerate(fields):
        try:
            number = int(field)
            if start is None:
                start = number
            elif stop is None:
                stop = number
        except ValueError:
            if (strand is None) and ((field == '+') or (field == '.') or (field == '-') ):
                strand = field
            elif chrom is None:
                chrom = field
   
    return GenericInterval(chrom,start,stop,strand)



