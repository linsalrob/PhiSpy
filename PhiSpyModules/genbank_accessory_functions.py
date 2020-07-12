"""
These functions are to help parse genbank records.
"""

def feature_id(seq, feat):
    """
    Choose the appropriate id for the feature
    :param feat: the feature
    :return: the id
    """

    if 'locus_tag' in feat.qualifiers:
        return "|".join(feat.qualifiers['locus_tag'])
    elif 'protein_id' in feat.qualifiers:
        return '|'.join(feat.qualifiers['protein_id'])
    elif 'db_xref' in feat.qualifiers:
        return '|'.join(feat.qualifiers['db_xref'])
    else:
        return seq.id + "." + str(feat.location)
