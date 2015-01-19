from parameters import Parameters
from load import load_location_file, load_sample_file
from polygon import filter_individuals_based_on_location, get_polygon
from eems import write_all_files, run_all, filter_data


def run(params):
    location_data = load_location_file(params.loc, params.location_header)
    sample_data = load_sample_file(params.sample, params.sample_header)

    meta_data = sample_data.merge(location_data)
    polygon = get_polygon(params.polygon,
                          wrap=params.wrap)
    meta_data = filter_individuals_based_on_location(meta_data, polygon)
        
    filter_data(meta_data, params.bed)
    write_all_files(params, meta_data, polygon)
    run_all(params)


def run_cli():
    """run_cli

    runs the pipeline from the command line interface.
    see -h flag for arguments
    """
    params = Parameters.from_command_line()
    run(params)
