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

