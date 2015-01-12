from shapely.geometry.polygon import Polygon
from shapely.geometry import Point
import numpy as np
from load import wrap_america, unwrap_america
from utils.plink import run_plink
import pandas as pd
import os

TMP_PLINK = 'tmp_plink'


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


def create_points(meta_data):
    """ adds a point object for meta data """
    if 'POINTS' not in meta_data:
        meta_data['POINTS'] = [Point(a[1]['LAT'], a[1]['LONG']) for a in
                               meta_data.iterrows()]
    return meta_data


def filter_data(meta_data, bedfile, plink="plink"):
    """filter_data
    filters bed file to only keep the individuals in meta_data, uses plink
    
    the data is read from bedfile, and written to the file in the TMP_PLINK
    global variable
    
    Parameters
    ----------
    meta_data : pd.DataFrame
        pandas data frame with individuals to keep
    bedfile : path
        the bedfile to be filtered
    plink : path
        the plink executable to be used
    
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


def create_eems_files(meta_data=None, polygon=None, bed2diffs="bed2diffs",
                      bedfile=None, diffs_file=None,
                      eems_input_name="eems_input", 
                      eems_output_name="eems_output",
                      **kwargs):
    """create_eems_files

    creates the eems files
        - diff file with pairwise differences (using bed2diffs)
        - polygon file with the location to run run on
        - sample file with sampling locations
        - ini file with all remaining parameters

    if there is more than one argument passed to nDemes, a different ini file
    is created for each of them.
    
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
    if diffs_file is not None:
        filter_diffs_file(diffs_file, meta_data.IND, eems_input_name)
    elif bedfile is not None:
        create_diffs_file(bedfile=bedfile, bed2diffs=bed2diffs,
                          outname=eems_input_name)
    else:
        raise ValueError("either bedfile or diffs_file needs to be specified")

    create_polygon_file(polygon, eems_input_name)
    create_sample_file(meta_data, eems_input_name, order_file=eems_input_name)

    nSites = sum(1 for line in open("%s.bim" % bedfile))
    nDemes = kwargs['nDemes']
    del kwargs['nDemes']
    for nd in nDemes:
        ini_name = "%s_%s" % (eems_input_name, nd)
        out_name = "%s/%s" % (eems_output_name, nd)
        create_ini_file(ini_name, datapath=eems_input_name,
                        mcmcpath=out_name,
                        meta_data=meta_data,
                        nSites=nSites, n_demes=nd, **kwargs)

    kwargs['nDemes'] = nDemes


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


def filter_diffs_file(diffs_file, individuals, outname):
    """Filters an existing diffs an order file combo to only retain the
    individuals in `individuals`
    
    Parameters
    ----------
    diffs_file : path
        path to a diffs/order file (created by bed2diffs)
    individuals : list of strings
        the individuals to be retained
    outname : str
        output file name
    
    """
    diff_matrix = np.loadtxt("%s.diffs" % diffs_file)
    sample_order = pd.read_table("%s.order" % diffs_file, header=None, sep=" ")
    individuals = pd.DataFrame(individuals)
    individuals.columns = [0]

    to_keep = [np.where(i == sample_order[1])[0][0] for i in individuals[0]]
    new_diff_matrix = diff_matrix[to_keep][:, to_keep]
    individuals = individuals.merge(sample_order, how="left",
                                    left_on=0, right_on=1)

    np.savetxt("%s.diffs" % outname, new_diff_matrix, fmt="%f")
    individuals.to_csv("%s.order" % outname, header=False, index=False,
                       columns=('0_y', '0_x'), sep=" ")


def create_ini_file(ini_name, mcmcpath, datapath, meta_data, n_demes,
                    **kwargs):
    """create_ini_file
    
    Parameters
    ----------
    ini_name : path
        the name of the ini file to be created
    mcmcpath : path
        path for the mcmc output of eems
    datapath : path
        path where the input data is loaded from
    meta_data : pd.DataFrame
        data structure with sample meta data
    
    """
    kwargs['nSites'] = kwargs.get('nSites', 10000)
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
        f.write("%s = %s\n" % ('nDemes', n_demes))


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


def create_sample_file(meta_data, outname, order_file=None):
    """create_sample_file
    
    Parameters
    ----------
    meta_data : pd.DataFrame
        data frame of the meta data
    outname : path
        the file of the outputed sample file (without extension)
    order_file : path
        a file with the sample ordering (without extension .order).
        This file is created by bed2diffs. If present, the sample
        coordinates are aligned with this file
    """
    out = "%s.coord" % outname
    if order_file is not None:
        sample_order = pd.read_table("%s.order" % order_file, header=None,
                                     sep=" ")
        meta_data = sample_order.merge(meta_data, how="left", left_on=1,
                                       right_on='IND')

    meta_data.to_csv(out, sep=" ", header=False, index=False,
                     columns=('LAT', 'LONG'))


def write_all_files(args, meta_data, polygon):
    if args.diffs is None:
        print 'de novo'
        filter_data(meta_data, args.bed)
        create_eems_files(meta_data=meta_data, polygon=polygon,
                          bed2diffs=args.bed2diffs,
                          bedfile=TMP_PLINK,
                          eems_input_name=args.proj,
                          eems_output_name=args.output_folder,
                          **args.eems_args)
    else:
        print 'filtering'
        filter_data(meta_data, args.bed)

        create_eems_files(meta_data=meta_data, polygon=polygon,
                          bedfile=TMP_PLINK,
                          bed2diffs=args.bed2diffs,
                          diffs_file=args.diffs,
                          eems_input_name=args.proj,
                          eems_output_name=args.output_folder,
                          **args.eems_args)


def run_all(args):
    if not args.dry:
        run_eems(args.eems_snps, ini_file=args.proj, n_demes=args.nDemes)
        create_plots(plot_scripts_folder=args.eems + '/runeems_snps/plot/',
                     out_path=args.output_folder, in_path=args.proj,
                     tmp_script=args.tmp_folder + '/plot.r',
                     wrap=args.wrap)


def create_plots(plot_scripts_folder, out_path, in_path, tmp_script, wrap):
    plot_script = "%s/sed-plots.R" % plot_scripts_folder
    tpl = out_path, in_path, plot_scripts_folder, plot_script, tmp_script
    s = 'sed -e "s!MCMCPATH!%s!g; s!PLOTPATH!%s!g; s!EEMSWD!%s!g;" %s > %s' % (tpl)
    os.system(s)

    if wrap:
        print in_path
        coords_file = '%s.coord' % in_path
    data = np.loadtxt(coords_file)
    data[:, 1] = unwrap_america(data[:, 1])
    np.savetxt(coords_file, data, fmt="%f")
    os.system('Rscript %s' % tmp_script)


def run_eems(eems_exe, ini_file, n_demes, dry=False):
    s = ""
    for nd in n_demes:
        tpl = eems_exe, ini_file, nd
        s += "%s --params %s_%s.ini &\n" % tpl
    s += 'wait'
    print "EEMS-command:", s
    os.system(s)
