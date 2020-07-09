"""
An alternative search phmms method that does not require a separate reading of the file. We parse
the data record we already have.

Need to adjust the location of the phmms call in PhiSpy.py
"""

import os
import sys

from io import StringIO
import subprocess
from argparse import Namespace
from tempfile import NamedTemporaryFile
from .genbank_accessory_functions import feature_id

from .log_and_message import log_and_message

from Bio import SearchIO

def search_phmms(**kwargs):
    """
    This is an alternate implementation of the search phmms
    but we write to a temp file and then parse that.

    See hmms/run_hmmer.py in EdwardsLab git repo for timing tests
    of different ways of doing this!
    :param kwargs:
    :return:
    """

    self = Namespace(**kwargs)
    log_and_message(f"Running HMM profiles against {self.phmms}", c="GREEN", stderr=True, quiet=self.quiet)
    # write the amino acids to a named temporary file that we will unlink later
    aaout = NamedTemporaryFile(mode='w+t', dir=self.output_dir, delete=False)
    log_and_message(f"hmmsearch: writing the amino acids to temporary file {aaout.name}\n", c="GREEN", stderr=True, quiet=self.quiet)
    aaout.seek(0)
    all_features = {}
    for seq in self.record:
        for feat in seq.get_features('CDS'):
            aa = ""
            if 'translation' in feat.qualifiers:
                aa = feat.qualifiers['translation'][0]
            else:
                aa = str(feat.extract(seq).translate().seq)
                feat.qualifiers['translation'] = [aa]
            myid = feature_id(seq, feat)
            if myid in all_features:
                log_and_message(f"FATAL: {myid} is not a unique id. We need unique protein IDs to run hmmsearch\n", c="RED",
                                stderr=True, quiet=False)
                sys.exit(-1)
            aaout.write(f">{myid}\n{aa}\n")
            all_features[myid] = feat
    aaout.close()

    try:
        search = subprocess.Popen(["hmmsearch", '--cpu', '6', '-E', '1e-10', '--domE', '1e-5', '--noali', self.phmms, aaout.name], stdout=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        sys.stderr.write(f"Error running hmmscan:\n{e}\n")
        sys.exit(-1)

    hmmresult = search.communicate()[0]
    results = SearchIO.parse(StringIO(hmmresult.decode()), 'hmmer3-text')

    for res in results:
        for hit in res:
            if hit.id not in all_features:
                log_and_message(f"Found an HMM hit to {hit.id} but we don't have a record of that feature", c="RED", stderr=True,
                                quiet=False)
                continue
            if 'phmm' not in all_features[hit.id].qualifiers:
                all_features[hit.id].qualifiers['phmm'] = []
            all_features[hit.id].qualifiers['phmm'].append(f"{res.id}:{hit.evalue}")

    # remove the amino acids
    os.unlink(os.path.join(self.output_dir, aaout.name))
    log_and_message(f"Completed running HMM profiles against {self.phmms}", c="GREEN", stderr=True, quiet=self.quiet)
    return self.record