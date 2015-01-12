import argparse
from utils.path import make_full_path
import utils.parameters
import os


def make_full_paths(args):
    args.ind = make_full_path(args.data_folder, args.ind)
    args.sample = make_full_path(args.data_folder, args.sample)
    args.loc = make_full_path(args.data_folder, args.loc)

    if not os.path.exists(args.analysis_folder):
        os.makedirs(args.analysis_folder)

    if args.output_folder is None:
        args.output_folder = make_full_path(args.analysis_folder,
                                            'output')
    if not os.path.exists(args.output_folder):
        os.makedirs(args.output_folder)
    for nd in args.nDemes:
        if not os.path.exists("%s/%s" % (args.output_folder, nd)):
            os.makedirs("%s/%s" % (args.output_folder, nd))

    if args.input_folder is None:
        args.input_folder = make_full_path(args.analysis_folder,
                                           'input')
    if not os.path.exists(args.input_folder):
        os.makedirs(args.input_folder)

    if args.tmp_folder is None:
        args.tmp_folder = make_full_path(args.analysis_folder,
                                         'tmp')
    if not os.path.exists(args.tmp_folder):
        os.makedirs(args.tmp_folder)

    args.proj = make_full_path(args.input_folder, args.proj)

    global TMP_PLINK
    TMP_PLINK = make_full_path(args.tmp_folder, 'tmp_plink')


class Parameters(utils.parameters.Parameters):
    """ class that handles parameters I/O

    the main idea is that this class unifies the different ways a project can
    be loaded. Current options are:
        - a dict like object, passed to init
        - from an input file
        - from the command line
        
    the main input format is from the command line
    """

    def __init__(self, **kwargs):
        Parameters.create_parser()
        for k, v in kwargs.items():
            setattr(self, k, v)

    @staticmethod
    def create_from_dict(d, defaults=True):
        """create_from_dict
        
        allows creation of a set of parameters from a dictionary.
        
        
        Parameters
        ----------
        d : dict
            Dictionary with format 'args' => val
        defaults : bool
            if True, then argparse is run once to get default arguments
        
        Returns
        -------
        params : Parameters
            the parameters object
        """

        params = Parameters(**d)
        if defaults:
            defaults = Parameters.parser.parse_args("")
            for k, v in defaults.__dict__.iteritems():
                if not hasattr(params, k):
                    setattr(params, k, v)

        params.postprocess_args()
        return params

    @staticmethod
    def create_parser():
        """
        generates ArgumentParser and reads options from CLI

        use -h flag for details
        """

        if hasattr(Parameters, "parser"):
            return
        parser = argparse.ArgumentParser("eems_pipeline", 
                                         fromfile_prefix_chars='@')
        parser.add_argument('--dry', default=False, action='store_true',
                            help="if set, files are created, but eems is not ran")
        parser.add_argument('--bed', '--bfile',
                            help="The common name of the bed file to be used")
        parser.add_argument('--diffs', default=None,
                            help='''A diffs file from a previous bed2diffs run,
                            without extension. If not None, the diffs file is
                            filtered rather than a new diffs file created from the
                            bed file.''')

        parser.add_argument('--loc', default=None,
                            help=""" File with location information.
                            Should have a column named `pop` and columns named
                            `latitude` and `longitude` or a variation thereof.""")
        parser.add_argument('--loc-has-no-header', 
                            action='store_false', dest="location_header", help="""add this flag if the
                            location file (--loc) doesn't have a header line. In
                            this case, the first three columns are assumed to be
                            the
                            population id, latitude and longitude, respectively."""
                            )
        parser.add_argument('--loc-has-header', default=True,
                            action='store_true', dest="location_header", help="""add this flag if the
                            location file (--loc) does have a header line. In
                            this case, the first three columns are assumed to be
                            the
                            population id, latitude and longitude, respectively."""
                            )

        parser.add_argument('--sample-has-no-header', 
                            action='store_false', dest="sample_header", help="""add this flag if the
                            sample file (--sample) doesn't have a header line. In
                            this case, the first three columns are assumed to be
                            the
                            sample id and population id , respectively."""
                            )
        parser.add_argument('--sample-has-header', default=True,
                            action='store_true', dest="sample_header", help="""add this flag if the
                            sample file (--sample) does have a header line. In
                            this case, the first three columns are assumed to be
                            the
                            sample id and population id, latitude and longitude, respectively."""
                            )
        parser.add_argument('--sample', default=None,
                            help=""" File with sample information.
                            Should have a column named `sample` and columns named
                            `pop` or a variation thereof.""")

        parser.add_argument('--ind', default=None,
                            help=""" File with individual based 
                            location information.
                            Should have a column named `sample` and columns named
                            `latitude` and `longitude` or a variation thereof.""")
        parser.add_argument('--ind-has-no-header', 
                            action='store_false', dest="ind_header", help="""add this flag if the
                            individual file (--ind) doesn't have a header line. In
                            this case, the first three columns are assumed to be
                            the
                            individual id, latitude and longitude, respectively."""
                            )
        parser.add_argument('--ind-has-header', default=True,
                            action='store_true', dest="ind_header", help="""add this flag if the
                            ind file (--ind) does have a header line."""
                            )

        parser.add_argument('--eems', '--eems-folder',
                            default="/data/eems-project/eems/",
                            help="The base folder of your local eems installation")
        parser.add_argument('--polygon',
                            default="/data/eems-project/human_origins/america.polygon",
                            help="""A file with the polygon describing the region to
                            run eeems on.
                            """)
        parser.add_argument('--bed2diffs',
                            default=None,
                            help="""The bed2diffs executable, defaults to
                            bed2diffs/bed2diffs in the directory of --eems""")
        parser.add_argument('--eems_snps',
                            default=None,
                            help="""The eems_snps executable, defaults to
                            runeems_snps/runeems_snps in the directory of --eems"""
                            )
        parser.add_argument('--input-folder', '--input_folder',
                            default=None,
                            help="""the folder where all eems input files will
                            be stored. (default: ./input )"""
                            )
        parser.add_argument('--output-folder', '--output_folder',
                            default=None,
                            help="""the folder where all eems output files will
                            be stored. (default: ./output )"""
                            )
        parser.add_argument('--tmp-folder', '--tmp_folder',
                            default=None,
                            help="""the folder where tempory output files will
                            be stored. (default: ./tmp )"""
                            )
        parser.add_argument('--analysis-folder', '--analysis_folder',
                            default='./',
                            help="""the folder where all analysis will be stored
                             (default: . )"""
                            )
        parser.add_argument('--data-folder', '--data_folder',
                            default=".",
                            help="""the folder where all data files are read from
                            (default: . )"""
                            )
        parser.add_argument('--proj', '--proj-name', '--proj_name',
                            default='eems_proj',
                            help="""The name of the output files""")

        parser.add_argument('--wrap', '--wrap_america',
                            default=True,
                            help="""Should all corinates be wrapped s.t. the
                            americas
                            appear in the east?
                            """)

        parser.add_argument('--nDemes',
                            default=[100], nargs='*',
                            help="""eems arg: number of demes in the model
                            (default: 100 )"""
                            )
        parser.add_argument('--diploid',
                            default='true',
                            help="""eems arg: is diploid? (true/false)
                            (default: true )"""
                            )
        parser.add_argument('--numMCMCIter', '--n_mcmc', '--n-mcmc',
                            default=20000,
                            help="""eems arg: number of MCMC iterations
                            (default: 20000 )"""
                            )
        parser.add_argument('--numBurnIter', '--n_burn', '--n-burn',
                            default=10000,
                            help="""eems arg: number of burn in steps
                            (default: 10000 )"""
                            )
        parser.add_argument('--numThinIter', '--n_thin', '--n-thin',
                            default=99,
                            help="""eems arg: thinning interval
                            (default: 99 )"""
                            )

        Parameters.parser = parser

    @staticmethod
    def from_command_line():
        """loads arguments from command line
        
        Returns
        -------
        p : Parameters
            the parameters object read from the command line
        """

        Parameters.create_parser()
        parser = Parameters.parser
        params = parser.parse_args()

        p = Parameters(**params.__dict__)
        p.postprocess_args()
        return p

    def postprocess_args(p):
        """
        processes some args in p
        """
        p.eems_args = dict()
        p.eems_args['diploid'] = p.diploid.lower()
        p.eems_args['numMCMCIter'] = p.numMCMCIter
        p.eems_args['numBurnIter'] = p.numBurnIter
        p.eems_args['numThinIter'] = p.numThinIter
        p.eems_args['nDemes'] = p.nDemes

        make_full_paths(p)

        if p.eems_snps is None:
            p.eems_snps = "%s/runeems_snps/src/runeems_snps" % p.eems
        if p.bed2diffs is None:
            p.bed2diffs = "%s/bed2diffs/bed2diffs" % p.eems
