import types
from Bio import SeqFeature
from copy import copy


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

    def distance_between(self, locations):
        return -(locations[0].end - locations[1].start)

    def merge_or_split(self, feature):
        cutoff_distance           = 100
        feature.location_operator = None
        last_part_start           = feature.location.parts[0].start
        last_part_end             = feature.location.parts[0].end
        mor_features              = []

        if 'product' in feature.qualifiers:
            feature.qualifiers['product'][0] = 'Merged-or-split: ' + feature.qualifiers['product'][0]
        else:
            feature.qualifiers['product'] = ['Merged-or-split: not assigned']

        # in case parts are not written in increasing order
        tmp = []
        for i in range(len(feature.location.parts)):
            tmp.append((feature.location.parts[i].start, feature.location.parts[i]))
        feature.location.parts = []
        for t in sorted(tmp):
            feature.location.parts.append(t[1])

        for i in range(len(feature.location.parts) - 1):
            if self.distance_between(feature.location.parts[i : i+2]) < cutoff_distance:
                # do merging stuff here
                last_part_end = feature.location.parts[i + 1].end
            else:
                # do splitting stuff here
                if last_part_end == feature.location.parts[i].end:
                    nf = copy(feature)
                    nf.location = SeqFeature.FeatureLocation(last_part_start, last_part_end, feature.strand)
                    mor_features.append(nf)
                    last_part_start = feature.location.parts[i + 1].start
                else:
                    last_part_start = feature.location.parts[i].start
                last_part_end = feature.location.parts[i + 1].end
        nf = copy(feature)
        nf.location = SeqFeature.FeatureLocation(last_part_start, last_part_end, feature.strand)
        mor_features.append(nf)

        return mor_features

    def attach_methods(self, target):
        """This method allows attaching new methods to the SeqIO entry object

           Args:
               target: is the SeqIO object that will be attaching a method to
        """

        # if feature has a complex location then merge or split it
        new_features    = []
        joined_features = []
        for feature in target.features:
            if feature.location_operator == 'join':
                new_features.extend(self.merge_or_split(feature))
                joined_features.append(feature)
        for feature in joined_features:
            target.features.remove(feature)
        target.features.extend(new_features)

        # sort the features based on their location in the genome
        target.features.sort( key = lambda feature : tuple([min(feature.location.start, feature.location.end),  max(feature.location.start, feature.location.end)]) )

        def get_features(target, feature_type):
            for feature in target.features:
                feature.id       = " ".join(feature.qualifiers.get('locus_tag', [str(feature.location)]))
                feature.function = " ".join(feature.qualifiers.get('product',['unknown']))
                # set start and stop for the feature and swap for direction if needed
                feature.start    = int(feature.location.start) + 1
                feature.stop     = int(feature.location.end)
                feature.phmm     = [1.0] if 'phmm' not in feature.qualifiers else [float(x.split(':')[1]) for x in feature.qualifiers.get('phmm')]
                if feature.strand < 0:
                    feature.start, feature.stop = feature.stop, feature.start
                # if a feature type was provided, only return those features
                if not feature_type or feature.type == feature_type:
                    yield feature
        target.get_features = types.MethodType(get_features,target)

        def append_feature(target, feature):
            target.features.append(feature)

        target.append_feature = types.MethodType(append_feature, target)
