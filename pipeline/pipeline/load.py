import pandas as pd
import os


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
        for i, h in enumerate(header):
            if str(name).upper() == str(h).upper():
                return h
    raise ValueError('did not find correct header for %s.\
                     Please adjust file/ script'%allowed_names)


def wrap_america(data):
    return [i + 360 if i < -45 else i for i in data]


def unwrap_america(data):
    return [i - 360 if i >= 180 else i for i in data]
