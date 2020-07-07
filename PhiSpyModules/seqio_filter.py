import types
from Bio import SeqFeature
import copy

from .log_and_message import log_and_message
from Bio.SeqFeature import FeatureLocation, CompoundLocation

class SeqioFilter( list ):
    """This class is to allow filtering of the Biopython SeqIO record

    SeqIO.parse returns a generator object so anytime you want to perform
    an action on it, you must iterate through the entire list. This
    class adds the ability to filter and return only a subset of the
    features. Those are split into separate features or joined based
    on the distance between the parts and then sorted based on the
    start to end location.

    Note:
        To use simply pass a SeqIO.parse object to it and then when
        the object is called a keyword is passed to it and only those
        features matching the keyword are returned.
    Example:
        record = SeqioFilter(SeqIO.parse(infile))):
        #no change to standard SeqIO calls
        for entry in record:
            print(entry.id, entry.seq)
        #now we can get only certain features
        for cds in record.get_features('CDS'):
            print(cds)

    """

    def __init__( self, content ):
        self.n = 0
        self.get_n = dict()
        for n, item in enumerate(content):
            self.attach_methods(item)
            self.get_n[item.id] = n
            self.append(item)

    def __iter__(self):
        return self

    def __next__(self):
        try:
            item = self[self.n]
        except IndexError:
            self.n = 0
            raise StopIteration()
        self.n += 1
        return item

    def __call__(self, keyword=''):
        pass

    def get_entry(self, id):
        return self[self.get_n[id]]

    def merge_or_split(self, seq, feature, mindistance = 100):
        """
        Merge or split sequence features with compound locations
        :param seq: the sequence object
        :param feature: The feature with a compound location to merge or split
        :param mindistance: the distance with which they will be merged/split
        :return:
        """
        thisid = " ".join(feature.qualifiers.get('locus_tag', [str(feature.location)]))
        # find the index of this in the list. We could also pass this as an arg
        idx = seq.features.index(feature)

        # do we need to do anything
        if type(feature.location) != CompoundLocation:
            log_and_message(f"Error {thisid} does not appear to be a compound location\n",
                            c="RED", stderr=True, loglevel="WARNING")
            return

        if not feature.strand:
            # per the biopython docs, a compound feature with some parts on one strand
            # and some parts on the other are given a strand designation of None
            # https://biopython.org/DIST/docs/api/Bio.SeqFeature.CompoundLocation-class.html#__init__
            log_and_message(f"Error {thisid} compound location seems to be on both strands. We do not know how to handle this!\n",
                            c="RED", stderr=True, loglevel="WARNING")
            return

        log_and_message(f"merging/splitting {thisid} original location: {feature.location}")

        if 'product' in feature.qualifiers:
            feature.qualifiers['product'][0] = 'Merged-or-split: ' + feature.qualifiers['product'][0]
        else:
            feature.qualifiers['product'] = ['Merged-or-split: not assigned']


        # simplify our coding
        loc = feature.location
        # first we find the strand
        strand = loc.parts[0].strand

        # test that all locations are on the same strand
        # record an error if not
        for p in loc.parts:
            if p.strand != strand:
                msg = "Error: We can not handle compound locations on different strands. For {thisid} we have {loc}"
                log_and_message(msg, c="RED", stderr=True, loglevel="WARNING")
                return

        all_locs = []
        merged = loc.parts[0]
        # handle features on the + strand
        if strand > 0:
            for p in loc.parts[1:]:
                if p.start < merged.start:
                    # this feature spans a break
                    msg = (f"Feature {thisid} spans the origin: {loc}\n"
                           f"WARNING: THIS IS AN UNTESTED FEATURE!\n"
                           f"We have not thoroughly tested conditions where an ORF on the +ve strand appears to cross "
                           f"the origin of the contig. We would appreciate you posting an issue on GitHub and sending "
                           f"Rob a copy of your genome to test!\n")
                    log_and_message(msg, c="YELLOW", stderr=True, loglevel='WARNING')

                    all_locs.append(merged)
                    merged = p
                    continue
                if p.start > merged.start and p.start < merged.end:
                    merged = FeatureLocation(merged.start, p.end, strand)
                elif p.start > merged.end:
                    if p.start - merged.end > mindistance:
                        all_locs.append(merged)
                        merged = FeatureLocation(p.start - 1, p.end, strand)
                    else:
                        merged = FeatureLocation(merged.start, p.end, strand)
                else:
                    all_locs.append(merged)
                    merged = FeatureLocation(p.start, p.end, strand)
        # handle features on the -ve strand
        else:
            for p in loc.parts[1:]:
                if merged.start < p.start:
                    # this feature spans a break
                    msg = (f"Feature {thisid} spans the origin: {loc}\n"
                           f"WARNING: THIS IS AN UNTESTED FEATURE!\n"
                           f"We have not thoroughly tested conditions where an ORF on the -ve strand appears to cross "
                           f"the origin of the contig. We would appreciate you posting an issue on GitHub and sending "
                           f"Rob a copy of your genome to test!\n")
                    log_and_message(msg, c="YELLOW", stderr=True, loglevel='WARNING')

                    all_locs.append(merged)
                    merged = p
                    continue
                if p.end > merged.start and p.end < merged.end:
                    # trivial case, the ORFs overlap
                    merged = FeatureLocation(p.start, merged.end, strand)
                elif p.end < merged.start:
                    # more complex case, there is a gap
                    if merged.start - p.end > mindistance:
                        all_locs.append(merged)
                        merged = FeatureLocation(p.start, p.end - 1, strand)
                    else:
                        merged = FeatureLocation(p.start, merged.end, strand)
                else:
                    all_locs.append(merged)
                    merged = FeatureLocation(p.start, p.end, strand)

        if all_locs:
            # we have multiple features, so we need to add features
            # make sure we add the last feature
            all_locs.append(merged)
            log_and_message(f"We could not join the whole {feature.type} feature into a single new feature.")
            # replace the existing compound feature with the first one of the split features
            newfeat = feature
            newfeat.location = all_locs[0]
            seq.features[idx] = newfeat
            # append the other features
            for f in all_locs[1:]:
                newfeat = copy.deepcopy(feature)
                newfeat.location = f
                seq.features.append(newfeat)
                log_and_message(f"Appended part of a multiple {newfeat.type} feature {thisid} loc: {f}\n")
        else:
            # we just replace the old feature with the new
            # update the location
            feature.location = merged
            # add the new feature
            seq.features[idx] = feature
            log_and_message(f"Created a single {feature.type} feature: {thisid} loc: {merged}\n")

    def attach_methods(self, target):
        """This method allows attaching new methods to the SeqIO entry object

           Args:
               target: is the SeqIO object that will be attaching a method to
        """

        # resolve compound locations!
        for feat in target.features:
            if type(feat.location) == CompoundLocation:
                self.merge_or_split(target, feat)

        # sort the features based on their location in the genome
        target.features.sort( key = lambda feature : tuple([min(feature.location.start, feature.location.end),
                                                            max(feature.location.start, feature.location.end)]) )

        def get_features(target, feature_type):
            for feature in target.features:
                feature.id       = " ".join(feature.qualifiers.get('locus_tag', [str(feature.location)]))
                feature.function = " ".join(feature.qualifiers.get('product',['unknown']))
                # set start and stop for the feature and swap for direction if needed
                feature.start    = int(feature.location.start) + 1
                feature.stop     = int(feature.location.end)
                feature.phmm     = [1.0] if 'phmm' not in feature.qualifiers else [float(x.split(':')[1]) for x in feature.qualifiers.get('phmm')]
                if not feature.strand:
                    # silently ignore compound features on two strands. Hopefully the error is recorded
                    # in split and merge!
                    continue
                if feature.strand < 0:
                    feature.start, feature.stop = feature.stop, feature.start
                # if a feature type was provided, only return those features
                if not feature_type or feature.type == feature_type:
                    yield feature
        target.get_features = types.MethodType(get_features,target)

        def append_feature(target, feature):
            target.features.append(feature)

        target.append_feature = types.MethodType(append_feature, target)
