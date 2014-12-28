"""
Simple computing of EEMS surfaces from a bed file and files with coordinate
info.

The possible location options are:

a) location file and sample file seperate:
    - in this case, we have one file that maps individuals to populations
      and another file that maps populations to geographic locations

b) ind file:
    - this file is assumed to already contain the location for each individual
"""

import numpy as np
import pandas as pd
import os
import argparse
# from geopy.geocoders import Nominatim
from shapely.geometry.polygon import Polygon
from shapely.geometry import Point

from bedmerger.bedmerger.utils import run_plink

TMP_PLINK = 'tmp_plink'


def parse_cmdline():
    parser = argparse.ArgumentParser("eems_pipeline")
    parser.add_argument('--bed', '--bfile',
                        help="The common name of the bed file to be used")

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
                        default=100,
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

    args = parser.parse_args()
    args.eems_args = dict()
    args.eems_args['diploid'] = args.diploid.lower()
    args.eems_args['numMCMCIter'] = args.numMCMCIter
    args.eems_args['numBurnIter'] = args.numBurnIter
    args.eems_args['numThinIter'] = args.numThinIter
    args.eems_args['nDemes'] = args.nDemes

    make_full_paths(args)

    if args.eems_snps is None:
        args.eems_snps = "%s/runeems_snps/src/runeems_snps" % args.eems
    if args.bed2diffs is None:
        args.bed2diffs = "%s/bed2diffs/bed2diffs" % args.eems
    return args


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


def make_full_path(path, fi):
    if fi is None:
        return fi
    if os.path.isabs(fi):
        return fi
    else:
        return os.path.sep.join((path, fi))


def load_location_file(location_file, location_has_header=True, format=None,
                       column_names=None, wrap=True,
                       **keywords):
    """load_location_file
    
    this function loads the location data using lat, long and pop_id. 
        The approach here is that I'll draw them from a list of possible ids,
        but the resulting data frame will always have three cols named POP, LAT, LONG.

    Parameters
    ----------
    location_file : path
        the file with the path
    location_has_header : bool
        does the location file has a header row?
    format : str
        the format of the file, allowed formats are xls, xlsx, or csv. If not
        given, it is guessed based on extension. Otherwise, read_table is used.
    wrap : bool
        should coordinates mapping to the americas be wrapped around?
    column_names : (str, str, str)
        the names of the columns to be used for population id, latitude and
        longitude, respectively.
    
    Returns
    -------
    location_data : pandas.DataFrame
        a table containing the location information
    """

    read_function = get_read_fun_from_extension(location_file, format)

    if location_has_header:
        location_data = read_function(location_file, **keywords)
    else:
        location_data = read_function(location_file, header=None, **keywords)
        headers = list(location_data.columns)
        headers[:3] = ['POP', 'LATITUDE', 'LONGITUDE']
        location_data.columns = headers



    header = location_data.columns.values

    
    if column_names is None:
        # possible ids for population, possibly add more
        allowed_pop_names = ['POP', 'ID', 'POP_ID', 'POP_NAME', 'POP_ORIG',
                             'POPULATION', 'ECOTYPE_ID',
                             "verbose Population ID"]
    
        # possible ids for latitude & longitude
        allowed_lat_names = ['LAT', 'LATITUDE', 'Y', 'LAT-ITUDE']
        allowed_long_names = ['LONG', 'LONGITUDE', 'X', 'LON-GI-TUDE']
    
        POP = get_first_id(header, allowed_pop_names)
        LAT = get_first_id(header, allowed_lat_names)
        LONG = get_first_id(header, allowed_long_names)

    else:
        POP, LAT, LONG = column_names

    location_data = location_data[[POP, LAT, LONG]]
    location_data.columns = ['POP', 'LAT', 'LONG']
    
    if wrap:
        location_data['LONG'] = wrap_america(location_data['LONG'])

    return location_data.drop_duplicates()


def get_read_fun_from_extension(file_name, format):
    """
    gets the proper pd function to read a data set
    """
    if format is None:
        format = os.path.splitext(file_name)[1]

    if format.startswith("."):
        format = format[1:]
    
    if format == "xlsx" or format == "xls":
        read_function = pd.read_excel
    elif format == "csv":
        read_function = pd.read_csv
    else:
        read_function = pd.read_table

    return read_function


def get_first_id(header, allowed_names):
    """private function that gets the first match from the possible list of matches"""
    for name in allowed_names:
        for i, h in enumerate( header ):
            if str(name).upper() == str(h).upper():
                return h
    raise ValueError('did not find correct header for %s.\
                     Please adjust file/ script'%allowed_names)


def load_sample_file(sample_file, sample_has_header=True, format=None,
                     column_names=None, **kwargs):
    """load_sample_file
    
    this file loads the sample ids that are used when analyzing VCF/bed files

    Parameters
    ----------
    sample_file : file or filename
        The file to assign samples to populations
    sample_has_header : bool
        does the sample have a header row
    format : str
        the format of the file, allowed formats are xls, xlsx, or csv. If not
        given, it is guessed based on extension. Otherwise, read_table is used.
    column_names : (str, str)
        the names of the columns to be used for sample id and population id,
        respectively.
    
    Returns
    -------

    sample_data : pd.DataFrame
        a data frame with columns 'SAMPLE' and `POP` that assigns samples
        to populations
    """
    read_function = get_read_fun_from_extension(sample_file, format)

    if sample_has_header:
        sample_data = read_function(sample_file, **kwargs)
    else:
        sample_data = read_function(sample_file, header=None, **kwargs)
        sample_data.columns[:2] = ['SAMPLE', 'POP']

    header = sample_data.columns.values
    
    if column_names is None:
        # possible ids for population
        allowed_pop_names = ['POP', 'POP_ID', 'POP_NAME', 'POP_ORIG',
                             'POPULATION', 'SIMPLE_MAJORITY']
        
        # possible ids for individuals
        allowed_ind_names = ['SAMPLE', 'IID', 'IND', 'INDIVIDUAL',
                             'INDIVIDUALS', 'ID', 'SAMPLE']
    
        POP = get_first_id(header, allowed_pop_names)
        IND = get_first_id(header, allowed_ind_names)
    else:
        POP, IND = column_names

    sample_data = sample_data[[IND, POP]]
    sample_data.columns = ['IND', 'POP']

    return sample_data


def get_country_from_coords(coords):
    """
    uses the openstreetmap api to look up the country of a point

    
    Parameters
    ----------
    coords : Point or (float, float)
        a set of coordinates (either a tuple or point object)
    
    Returns
    -------
    the country code of the location, None if not on land
    """
    geolocator = Nominatim()
    location = geolocator.reverse(coords, timeout=10)

    address = location.address
    if address is None:
        return None

    country = address['country_code']
    return country


def get_polygon(polygon, expand=0, wrap=True):
    """get_polygon
    gets the boundary polygon tu run eems on
    
    Parameters
    ----------
    polygon : int or str
        if str, file with polygon in eems format
        if set to `0`, get a  convex hull around the sample points
        if set to `1`, use the countries the points are in
        if set to `2`, get continent data
    expand : float
        how much the polygon should be expanded
    wrap : bool
        should coordinates in the americas be wrapped s.t they appear in the
        east?
    
    Returns
    -------
    polygon: shapely Polygon object
        the polygon to be used in the eems analysis
    """
    if polygon == 0:
        pass
    elif polygon == 1:
        raise NotImplementedError("country mode nyi")
    elif polygon == 2:
        raise NotImplementedError("continent mode nyi")
    else:
        polygon = np.loadtxt(polygon)
        polygon[:, 1] = wrap_america(polygon[:, 1])
        return Polygon(polygon)


def create_points(meta_data):
    """ adds a point object for meta data """
    if not 'POINTS' in meta_data:
        meta_data['POINTS'] = [Point(a[1]['LAT'], a[1]['LONG']) for a in
                               meta_data.iterrows()]
    return meta_data


def filter_individuals_based_on_location(meta_data, polygon):
    """filter_individuals
    only retains individuals that are inside the polygon
    
    Parameters
    ----------
    meta_data : pd.data_frame
        data frame with individuals and their location
    polygon : list of (int, int)
        a polygon, output from get_polygon
    
    Returns
    -------
    filtered_meta_data : pd.DataFrmae
        a data frame only with individuals used for analysis
    """
    create_points(meta_data)
    to_keep = [polygon.contains(pt) for pt in meta_data['POINTS']]
        
    return meta_data[to_keep]


def filter_data(meta_data, bedfile, plink="plink"):
    """filter_data
    filters bed file to only keep the individuals in meta_data, uses plink
    
    Parameters
    ----------
    meta_data : pd.DataFrame
        pandas data frame with individuals to keep
    bedfile : path
        the bedfile to be filtered
    plink : path
        the plink executable to be used
    
    Returns
    -------
    """
    include_name = '%s.incl' % TMP_PLINK

    fam = pd.read_table("%s.fam" % bedfile, header=None,
                        skipinitialspace=True, sep=" ")
    fam.columns = ['FAM', 'IND', 'a', 'b', 'c', 'd']

    extract_data = meta_data.merge(fam, on='IND', how='left')
    extract_data.to_csv(include_name, sep=' ',
                        columns=('FAM', 'IND'),
                        header=None, index=None)

    flags = dict()
    flags['make-bed'] = ''
    flags['bfile'] = bedfile
    flags['out'] = TMP_PLINK
    flags['keep'] = include_name
    flags['indiv-sort'] = 'f %s' % include_name

    run_plink(plink, flags)


def create_eems_files(bedfile, meta_data, polygon, bed2diffs="bed2diffs",
                      outname="eems_input", destination="eems_output",
                      **kwargs):
    """create_eems_files 
    
    Parameters
    ----------
    bedfile : path
        Path to the bed files to be converted
    meta_data : pd.DataFrame
        the location data for all individuals
    polygon : list of (str, str)
        the polygon to be written
    bed2diffs : path
        the bed2diffs executable to be used
    
    """
    create_diffs_file(bedfile, bed2diffs, outname)
    create_polygon_file(polygon, outname)
    create_sample_file(meta_data, outname)

    nSites = sum(1 for line in open("%s.bim" % bedfile))
    create_ini_file(TMP_PLINK, destination, outname, meta_data,
                    nSites=nSites, **kwargs)


def create_polygon_file(polygon, outname):
    """create_polygon_file
    writes a polygon to outname, for input in eems
    
    Parameters
    ----------
    polygon : shapely Polygon
        the polygon to be written
        
    outname : str
        file name of the polygon
    
    Returns
    -------
    """
    points = list(polygon.exterior.coords)
    np.savetxt("%s.outer" % outname, points, fmt="%f")


def create_diffs_file(bedfile, bed2diffs, outname, nthreads=4):
    """create a file with pairwise differences using the bed2diffs executable
        
    Parameters
    ----------
    bedfile : path
        the bed file to run bed2diffs on
    bed2diffs : path
        path to the bed2diffs executable
    outname : str
        output file name
    
    """
    tpl = bed2diffs, nthreads, bedfile
    s = "%s --nthreads %d --bfile %s" % (tpl)
    s += " && mv %s.order %s.order " % (bedfile, outname)
    s += " && mv %s.diffs %s.diffs " % (bedfile, outname)
    print s
    os.system(s)


def create_sample_file(meta_data, outname):
    """create_sample_file
    
    Parameters
    ----------
    meta_data : type
        Desc
    outname : type
        Desc
     
    Returns
    -------
    """
    out = "%s.coord" % outname
    meta_data.to_csv(out, sep=" ", header=False, index=False,
                     columns=('LAT', 'LONG'))


def create_ini_file(ini_name, mcmcpath, datapath, meta_data, **kwargs):
                    kwargs['nSites'] = kwargs.get('nSites', 10000)
                    kwargs['nDemes'] = kwargs.get('nDemes', 200)
                    kwargs['nIndiv'] = len(meta_data)
                    kwargs['diploid'] = kwargs.get('diploid', 'True')
                    kwargs['numMCMCIter'] = kwargs.get('numMCMCIter', 20000)
                    kwargs['numBurnIter'] = kwargs.get('numBurnIter', 10000)
                    kwargs['numThinIter'] = kwargs.get('numThinIter', 99)
                    kwargs['datapath'] = datapath
                    kwargs['mcmcpath'] = mcmcpath

                    with open("%s.ini" % ini_name, 'w') as f:
                        for k, v in kwargs.iteritems():
                            f.write("%s = %s\n" % (k, v))


def filter_diffs_file(diffs_file, individuals, bed2diffs, outname):
    """Filters an existing diffs an order file combo to only retain the
    individuals in `individuals`
    
    Parameters
    ----------
    diffs_file : path
        path to a diffs file
    individuals : list of strings
        the individuals to be retianed
    bed2diffs : path
        path to the bed2diffs executable
    outname : str
        output file name
    
    Returns
    -------
    """
    pass


def create_plots(plot_scripts_folder, out_path, in_path, tmp_script, wrap):
    plot_script = "%s/sed-plots.R" % plot_scripts_folder
    tpl = out_path, in_path, plot_scripts_folder, plot_script, tmp_script
    s = 'sed -e "s!MCMCPATH!%s!g; s!PLOTPATH!%s!g; s!EEMSWD!%s!g;" %s > %s' % (tpl)
    os.system(s)

    if wrap:
        print in_path
        coords_file = '%s.coord' %in_path
    data = np.loadtxt(coords_file)
    data[:, 1] = unwrap_america(data[:, 1])
    np.savetxt(coords_file, data, fmt="%f")
    os.system('Rscript %s' % tmp_script)


def run_eems(eems_exe, ini_file):
    tpl = eems_exe, ini_file
    s = ("%s --params %s.ini" % tpl)
    print "EEMS-command:", s
    os.system(s)

def wrap_america(data):
    return [i + 360 if i < -45 else i for i in data]


def unwrap_america(data):
    return [i - 360 if i >= 180 else i for i in data]

if __name__ == "__main__":
    args = parse_cmdline()
    location_data = load_location_file(args.loc, args.location_header,
                                       wrap=args.wrap)
    sample_data = load_sample_file(args.sample, args.sample_header)
    meta_data = sample_data.merge(location_data)
    
    polygon = get_polygon(args.polygon, wrap=args.wrap)
    meta_data = filter_individuals_based_on_location(meta_data, polygon)
    filter_data(meta_data, args.bed)
    #create_eems_files(TMP_PLINK, meta_data, polygon, args.bed2diffs,
    #                  outname=args.proj,
    #                  destination=args.output_folder,
    #                  **args.eems_args)

    run_eems(args.eems_snps, ini_file=TMP_PLINK)
    create_plots(plot_scripts_folder=args.eems + '/runeems_snps/plot/',
                 out_path=args.output_folder, in_path=args.input_folder,
                 tmp_script=args.tmp_folder + '/plot.r',
                 wrap=args.wrap)
