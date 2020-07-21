import re
import os


def is_phage_func(func):
    func = func.lower()
    func = func.replace('-', ' ')
    func = func.replace(',', ' ')
    func = func.replace('/', ' ')
    a = re.split(' ', func)
    if (
            'integrase'     in a or
            'phage'         in a or
            'lysin'         in a or
            'endolysin'     in a or
            'holin'         in a or
            'capsid'        in a or
            'tail'          in a or
            'bacteriophage' in a or
            'prophage'      in a or
            'portal'        in a or
            'terminase'     in a or
            'tapemeasure'   in a or
            'baseplate'     in a or
            'virion'        in a or
            'antirepressor' in a or
            'excisionase'   in a or
            re.search(r"\b%s\b" % "tape measure", func) or
            re.search(r"\b%s\b" % "Cro-like repressor", func) or
            re.search(r"\b%s\b" % "CI-like repressor", func) or
            re.search(r"\b%s\b" % "rIIA lysis", func) or
            re.search(r"\b%s\b" % "rI lysis", func) or
            re.search(r"\b%s\b" % "rIIB lysis", func) or
            re.search(r"\b%s\b" % "base plate", func) or
            ("head" in a and "decoration" in a) or
            ("helix" in a and "turn" in a) or
            "HNH endonuclease" in func or
            "single-stranded DNA-binding protein" in func
    ):
        return True
    return False


def is_not_phage_func(x):
    """
    These are some annotations that are definitely not phages, but often look like them.

    :param x: the function to test
    :return: True if it is NOT a phage function.
        Note that False suggests it maybe a phage function, but may not be (we are not sure)!
    """

    x = x.lower()
    if (
            ('phage' in x and 'shock' in x) or
            ('trna' in x and 'synthase' in x) or
            "conjugal transfer" in x or
            "conjugative" in x or
            "flagella" in x or
            "flagellar" in x or
            "flagellin" in x or
            "flagellum" in x or
            "ribosomal protein" in x or
            "translation elongation factor" in x or
            "secy" in x or
            "summary phrase" in x or
            "dna binding domain" in x or
            "abortive infection bacteriophage resistance protein" in x
    ):
        return True
    return False


def is_unknown_func(x):
    x_lower = x.lower()
    if (
            (len(x) == 0) or
            # ('hypoth' in x_lower) or
            'mobile element protein' == x_lower or
            ('hypothetical' in x_lower) or
            ('conserved protein' in x_lower) or
            ('gene product' in x_lower) or
            ('interpro' in x_lower) or
            ('uncharacterized' in x_lower) or
            ('pseudogene' in x_lower) or
            ('similar to' in x_lower) or
            ('similarity' in x_lower) or
            ('glimmer' in x_lower) or
            ('unknown' in x_lower) or
            ('complete' in x_lower) or
            ('ensang' in x_lower) or
            ('unnamed' in x_lower) or
            ('expressed' in x_lower) or
            ('similar to' in x_lower) or
            (' identi' in x_lower) or
            ('ortholog of' in x_lower) or
            ('structural feature' in x_lower) or
            ('cds_' in x_lower) or
            ('predicted by psort' in x_lower) or
            ('AGR_' in x) or
            ('EG:' in x) or
            ('RIKEN' in x) or
            re.search('lmo\d+ protein', x_lower) or
            re.search('lmo\d+protein', x_lower) or
            re.search('B[sl][lr]\d', x_lower) or
            re.search('^U\d', x) or
            re.search('[a-zA-Z]{2,3}\|', x) or
            re.search('orf\d+', x_lower) or
            re.match('orf[^_]', x_lower) or
            re.match('predicted', x_lower) or
            re.match('bh\d+', x_lower) or
            re.match('y[a-z]{2,4}\\b', x) or
            re.match('[a-z]{2,3}\d+[^:\+\-0-9]', x_lower)
    ):
        return True
    return False


def downweighting_unknown_functions(self):
    for d in self.initial_tbl:
        if d[8] > 0 and is_unknown_func(d[1]):
            d[8] = 0.5
    return self.initial_tbl