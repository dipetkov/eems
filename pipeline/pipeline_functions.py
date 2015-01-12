import pandas as pd
import numpy as np
import wormtable as wormtable


def load_location_file( location_file, location_dir="./", **keywords ):
    """ this function loads the location data using lat, long and pop_id. 
        The approach here is that I'll draw them from a list of possible ids,
        but the resulting data frame will always have three cols named POP, LAT, LONG.
    """
    location_data = pd.read_table( location_file, **keywords)
    header = location_data.columns.values
    
    #possible ids for population, possibly add more depending on what people come up with
    allowed_pop_names = [ 'ID', 'POP_ID', 'POP_NAME', 'POP_ORIG', 'POPULATION', 'ECOTYPE_ID' ] 
    
    #possible ids for latitude & longitude
    allowed_lat_names = [ 'LAT', 'LATITUDE', 'Y' ]
    allowed_long_names = [ 'LONG', 'LONGITUDE', 'X' ]    
    
    def get_first_id( header, allowed_names):
        """private function that gets the first match from the possible list of matches"""
        for name in allowed_names:
            for i, h in enumerate( header ):
                if name.upper() == h.upper():
                    return h
        raise ValueError('did not find correct header for %s.\
                         Please adjust file/ script'%allowed_names)
    
    POP = get_first_id( header, allowed_pop_names )
    LAT = get_first_id( header, allowed_lat_names )
    LONG = get_first_id( header, allowed_long_names )

    location_data = location_data[ [POP, LAT, LONG] ]
    location_data.columns = [ 'POP', 'LAT', 'LONG' ]
    
    return location_data

def load_sample_file( sample_file, sample_dir="./", **kwargs):
    """ 
    this file loads the sample ids that are used when analyzing VCF/bed files
    """
    sample_data = pd.read_table( sample_file, **kwargs )
    header = sample_data.columns.values
    
    #possible ids for population, possibly add more depending on what people come up with
    allowed_pop_names = [ 'POP_ID', 'POP_NAME', 'POP_ORIG', 'POPULATION', 'SIMPLE_MAJORITY', 'SAMPLE' ] 
    
    #possible ids for individuals & families (optional)
    allowed_ind_names = [ 'IID', 'IND', 'INDIVIDUAL', 'INDIVIDUALS', 'ID', 'SAMPLE' ]
    allowed_family_names = [ 'FID', 'FAM', 'FAMILY', 'FAMILIES' ]    
    
    def get_first_id( header, allowed_names):
        """private function that gets the first match from the possible list of matches"""
        for name in allowed_names:
            for i, h in enumerate( header ):
                if name.upper() == h.upper():
                    return h
        raise ValueError('did not find correct header for %s.\
                         Please adjust file/ script'%allowed_names)
    
    POP = get_first_id( header, allowed_pop_names )
    IND = get_first_id( header, allowed_ind_names )
    #FAMILY = get_first_id( header, allowed_family_names )

    sample_data = sample_data[ [IND, POP] ]
    sample_data.columns = [ 'IND', 'POP' ]


    return sample_data

def vcf2wt( vcf_file, wt_file, forced=True ):
    if forced:
        wormtable.os.system( 'vcf2wt -f %s %s'%( vcf_file, wt_file )  )
    else:
        wormtable.os.system( 'vcf2wt %s %s'%( vcf_file, wt_file )  )

def load_vcf_file ( vcf_file, wt_file = 'data.wt', individuals=[]):
    """ loads vcf """
    return load_wt_file ( wt_file, individuals )
    
def load_wt_file( wt_file, individuals=[] ):
    """ loads wormtable file
    """
    genotypes = [ind+'.GT' for ind in individuals ]
    table = wormtable.open_table( wt_file )
    tc = table.cursor( genotypes )
    return tc, individuals


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

    long_add = ( max( long_d ) - min( long_d ) ) * long_pc
    lat_add = ( max( lat_d ) - min( lat_d ) ) * lat_pc

    long_limits = min( long_d ) * (1 - long_pc), max( long_d ) + (1 + long_pc)
    lat_limits = min( lat_d ) * (1 + lat_pc), max( lat_d ) + (1 + lat_pc)

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
        if i % 10000 == 0: print i

        genotypes, counts = list_genotype2allele( row )
        tmp_dist_matrix = genetic_distance( genotypes )
        n_snps_matrix += np.outer( counts, counts )
        dist_matrix += tmp_dist_matrix
    
    dist_matrix = np.array( dist_matrix, dtype="float")
    dist_matrix /= ( n_snps_matrix / 4.)
        
    return dist_matrix, (n_inds, np.max( n_snps_matrix ))


"""
import cyvcf
from collections import Countern

def read_vcf_slow( data_file ):
    vcf_handle = vcf.VCFReader( data_file )
    n_processes
    samples = vcf_handle.samples
    samples_to_keep = [  i for i, s in enumerate( samples )     if s in merged_data[ 'IND' ].tolist()  ]
     
    def get_snp_allele_freq( gt ):
        if gt[0] == '.' or gt[2] == '.':
            return np.nan
        return int( gt[0] ) + int( gt[2] )
     
    
    n_samples = len(samples_to_keep)
    pairs = []
    for i in xrange( n_samples ):
            for j in xrange( n_samples):
                pairs.append(  ( i, j )  )
                    
                
    pw_comps = np.zeros( (n_samples, n_samples ) )
    pw_diff = np.zeros( (n_samples, n_samples ) )
        
    p = Pool( n_processes )
    for line,snp in enumerate(vcf_handle):
        current_snp = [ get_snp_allele_freq( snp.genos[i][0] ) for i in samples_to_keep ]
            
                
        
        def pw_fun( pair ):
            i, j = pair
            valid_comparison = not np.isnan( current_snp[i] - current_snp[j] )
            if valid_comparison:
                pw_comps[i,j] += 1
                pw_diff[i,j] += abs( current_snp[i] - current_snp[j] )
            
        #p.map( pw_fun, pairs )
        map( pw_fun, pairs )
        if line % 100 == 0:
            print line
def read_vcf( data_file):
    vcf_handle = cyvcf.Reader( open(data_file, 'r' ) )
    samples = vcf_handle.samples
    
    samples_to_keep =  merged_data[ 'IND' ].tolist()  
    n_samples = len( samples_to_keep )
        
    c_missing = Counter()
    c_diff1 = Counter()
    c_diff2 = Counter()
        
     
    
    for i,line in enumerate(vcf_handle):
        hom_ref = line.get_hom_refs()
        hom_alt = line.get_hom_alts()
        hets = line.get_hets()
        no_call = line.get_unknowns()
            
        for s1 in hom_ref:
            if s1.sample not in samples_to_keep: continue
            for s2 in hets:
                if s2.sample not in samples_to_keep: continue
                c_diff1.update( [(s1.sample, s2.sample)] )
            for s2 in no_call:
                if s2.sample not in samples_to_keep: continue
                c_missing.update( [(s1.sample, s2.sample)] )
            for s2 in hom_alt:
                if s2.sample not in samples_to_keep: continue
                c_diff2.update( [(s1.sample, s2.sample)] )
                    
        for s1 in hom_alt:
            if s1.sample not in samples_to_keep: continue
        
            for s2 in hets:
                if s2.sample not in samples_to_keep: continue
                c_diff1.update( [(s1.sample, s2.sample)] )
            for s2 in no_call:
                if s2.sample not in samples_to_keep: continue
                c_missing.update( [(s1.sample, s2.sample)] )
            
        for k,s1 in enumerate( no_call) :
            if s1.sample not in samples_to_keep: continue
        
            for s2 in hets:
                if s2.sample not in samples_to_keep: continue
                c_missing.update( [(s1.sample, s2.sample)] )
            for s2 in no_call[k+1:]:
                if s2.sample not in samples_to_keep: continue
                c_missing.update( [(s1.sample, s2.sample)] )     
                    
        for k,s1 in enumerate(hets):
            if s1.sample not in samples_to_keep: continue
        
            for s2 in hets:
                if s2.sample not in samples_to_keep: continue
                c_diff1.update( [(s1.sample, s2.sample)] )
            
        
        if i % 10000 == 0:
            print i
                
    
    n_snp = i
        
    pw_dist = np.zeros( (n_samples, n_samples) )
    
    for i, s1 in enumerate( samples_to_keep ):
        for j, s2 in enumerate( samples_to_keep) :
            pw_dist[i, j] = 2 *c_diff2[s1, s2] + 2 * c_diff2[s2, s1]
            pw_dist[i, j] += c_diff1[s1,s2] + c_diff1[s2, s1]
            pw_dist[i, j] /= ( n_snp - c_missing[s1, s2] - c_missing[s2, s1] ) 
         
    return( pw_dist, ( n_samples, n_snp )  )
"""    

                            
           
