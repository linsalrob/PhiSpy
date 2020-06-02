# read and set the version

from .seqio_filter import SeqioFilter
from .makeTest import make_set_train, measure_features
from .classification import call_randomforest, make_initial_tbl
from .protein_functions import downweighting_unknown_functions, is_phage_func, is_unknown_func
from .evaluation import fixing_start_end
from .helper_functions import get_args, print_list, is_gzip_file
from .pathtype import PathType
from .writers import write_gff3, write_genbank, write_phage_and_bact, write_prophage_tbl, write_prophage_tsv
from .writers import prophage_measurements_to_tbl
from .writers import log_and_message
from .search_phmms import search_phmms
from .formatting import message
from .errors import ColorNotFoundError, NoBasesCounted

from .version import __version__

__all__ = ['SeqioFilter',
           'make_set_train', 'measure_features',
           'call_randomforest', 'make_initial_tbl',
           'downweighting_unknown_functions', 'is_phage_func', 'is_unknown_func',
           'fixing_start_end',
           'get_args', 'print_list', 'is_gzip_file',
           'PathType',
           'write_gff3', 'write_genbank', 'write_phage_and_bact', 'write_prophage_tbl', 'write_prophage_tsv',
           'prophage_measurements_to_tbl', 'log_and_message',
           'search_phmms',
           'message',
            'ColorNotFoundError', 'NoBasesCounted'
           ]
