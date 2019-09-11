import types


class SeqioFilter( list ):
    """This class is to allow filtering of the Biopython SeqIO record

    SeqIO.parse returns a generator object so anytime you want to perform
    an action on it, you must iterate through the entire list. This
    class add the ability to filter and return only a subset of the
    features.

    Note:
        To use simply pass a SeqIO.parse object to it and then when
        the object is called a keyword is passed to it and only those
        features matching the keyword are returned.
    Example:
        record = SeqioFilter(SeqIO.parse(infile)):
        #no change to standard SeqIO calls
        for entry in record:
            print(entry.id, entry.seq)
        #now we can get only certain features
        for cds in record.get_feature('CDS'):
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

    def __call__( self, keyword='' ):
        pass

    def get_entry(self, id):
        return self[self.get_n[id]]
 
    def attach_methods(self, target):
        """This method allows attaching new methods to the SeqIO entry object

           Args:
               target: is the SeqIO object that will be attaching a method to
        """
        def get_features(target, feature_type):
            for feature in target.features:
                feature.id       = " ".join(feature.qualifiers.get('locus_tag', [str(feature.location)]))
                feature.function = " ".join(feature.qualifiers.get('product',['unknown']))
                feature.start    = int(feature.location.start) + 1
                feature.stop     = int(feature.location.end)
                if feature.strand < 0:
                    feature.start, feature.stop = feature.stop, feature.start
                if not feature_type or feature.type == feature_type:
                    yield feature
        target.get_features = types.MethodType(get_features,target)

