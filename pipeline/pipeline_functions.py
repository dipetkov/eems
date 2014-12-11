import pandas as pd
import numpy as np
import wormtable as wormtable
import os

from pipeline_load import *




def genotype2allele( genotype ):
    """ formats the genotype string into allele counts
    
    Parameters
    ----------
    genotype : string
        the genotype to format. In vcf format, e.g 0|1
    
    Returns
    -------
    out : int
        the number of non-reference alleles
    """
            
    try:
        return ( int(genotype[0]) + int(genotype[2]) , 2)
    except:
        return 0 , 0

def list_genotype2allele( genotype_list ):
    return np.array([ genotype2allele( ind ) for ind in genotype_list ]).transpose()

def genetic_distance( allele_counts ):
    """ calculates genetic distance according to eems paper"""
    return np.subtract.outer( allele_counts, allele_counts ) ** 2

    
def write_eems_input_files( name, location_data, pw_dist, n_samp_snp, 
        eems_options):
    """writes the eems input files.
    
    Parameters
    ----------
    name : string,
        the common name of all files
    location_data : array_like, with headers 'LAT' and 'LONG'
        the data object with `n` locations
    pw_dist : array_like, `n` x `n` matrix
        the pairwise distance matrix used by eems
    n_samp_snp : tuple (int, int)
        tuple with the number of samples and the number of snp
    eems_options : dict[str] => var 
        tuple with the number of samples and the number of snp
    """
    np.savetxt(name+'.diffs', pw_dist )
    location_data.to_csv(name+'.coord', sep=' ', header=False, columns=['LONG',
        'LAT'], index=False)
    
    long_pc, lat_pc = eems_options['long_add'],eems_options['lat_add']
    long_d, lat_d = location_data['LONG'], location_data['LAT']

    #the absolute amount to be added
    long_add = ( max( long_d ) - min( long_d ) ) * long_pc
    lat_add = ( max( lat_d ) - min( lat_d ) ) * lat_pc

    long_limits = min( long_d ) - long_add,   max( long_d ) + long_add
    lat_limits  = min( lat_d  ) - lat_add ,   max( lat_d  ) + lat_add

    with open( name + '.dimns', 'w' ) as settings_file:
        settings_file.write( "%f %f\n"%( long_limits ) )
        settings_file.write( "%f %f\n"%( lat_limits ) )
        settings_file.write( "%f %f\n"%(  n_samp_snp  ) )


def get_pw_dist ( wt_file, individuals=[]):
    """ calculates the pairwise distance matrix of all genotypes in a wt
    file for individuals individuals
    
    Parameters
    ----------
    wt_file : wormtable table
        the data structure containing the genetic data, see load_wt_data
    individuals : array_like, string
        the ids of the individuals for which we calculate the pairwise distance
    
    Returns
    -------
    dist_matrix : array_like, float
        the distance matrix calculated
    n_inds_snps : (int, int)
        the number of individuals and the max number of SNP comparisons
    """
    n_inds = len( individuals )
    dist_matrix = np.zeros( (n_inds, n_inds), dtype="int")
    n_snps_matrix = np.zeros_like( dist_matrix, dtype="int" )
     
    tc, individuals = load_wt_file( wt_file, individuals=individuals )
    for i, row in enumerate( tc ):
        if i % 1000 == 0: print i

        genotypes, counts = list_genotype2allele( row )
        tmp_dist_matrix = genetic_distance( genotypes )
        n_snps_matrix += np.outer( counts, counts )
        dist_matrix += tmp_dist_matrix
    
    dist_matrix = np.array( dist_matrix, dtype="float")
    dist_matrix /= ( n_snps_matrix / 4.)
        
    return dist_matrix, (n_inds, np.max( n_snps_matrix ))

def test ( wt_file, individuals=[]):
    """ calculates the pairwise distance matrix of all genotypes in a wt
    file for individuals individuals
    
    Parameters
    ----------
    wt_file : wormtable table
        the data structure containing the genetic data, see load_wt_data
    individuals : array_like, string
        the ids of the individuals for which we calculate the pairwise distance
    
    Returns
    -------
    dist_matrix : array_like, float
        the distance matrix calculated
    n_inds_snps : (int, int)
        the number of individuals and the max number of SNP comparisons
    """
    n_inds = len( individuals )
    dist_matrix = np.zeros( (n_inds, n_inds), dtype="int")
    n_snps_matrix = np.zeros_like( dist_matrix, dtype="int" )
     
    tc, individuals = load_wt_file( wt_file, individuals=individuals )
    for i, row in enumerate( tc ):
        if i % 1000 == 0: print i

        genotypes, counts = list_genotype2allele( row )
    
        
