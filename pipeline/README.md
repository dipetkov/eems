
# EEMS: python-pipeline

this is a basic pipeline to facilitate running eems. It mainly facilitates
running eems on subsets of a data set and setting parameters.

run `python pipeline_cli.py -h` to get an overview of all the options. A
possibly not up-to-date version of the output is replicated at the end of this
readme.

There are two ways to run this: 
 - an ipython notebook (eems_pipeline_new.ipynb), which allows to run steps
   interactively
 - a python command line interface (pipeline_cli.py)


The main input files are:
 - loc:
    a file that assigns each population a physical location
 - sample:
    assign each individual a population
 - bed:
    a bed/bim/fam file combination that contains data for all individuals
 - polygon:
    a boundary polygon for the area to analyze. Only individuals inside the 
    polygon are considered for the analysis.

###Running notebook:
simply copy the notebook (to keep the original), and adjust the parameters in 
the second input cell.

###Running from cli:
parameters can be set either directly on cli, or using arguments files. If read 
from a file, they need to be prepended by an "@" symbol. For the example, you
can run

    python pipeline_cli.py @example/default.args @example/example.args

which will read arguments from the files in `example/default.args` and 
`example/example.args`. This hopefully facilitates reproduceability.

###Install:

I am still experimenting with modules, but in general, I require the following 
python modules:
- numpy
- pandas
- shapely

to install the notebook, simply
```
sudo apt-get install ipython-notebook python-numpy python-matplotlib python-pandas python-pygments
sudo apt-get install pandoc
ipython notebook
``` 
then  view the notebook in a webbrowser

### Usage (2015-01-12)


    usage: eems_pipeline [-h] [--dry] [--bed BED] [--diffs DIFFS] [--loc LOC]
                         [--loc-has-no-header] [--loc-has-header]
                         [--sample-has-no-header] [--sample-has-header]
                         [--sample SAMPLE] [--ind IND] [--ind-has-no-header]
                         [--ind-has-header] [--eems EEMS] [--polygon POLYGON]
                         [--bed2diffs BED2DIFFS] [--eems_snps EEMS_SNPS]
                         [--input-folder INPUT_FOLDER]
                         [--output-folder OUTPUT_FOLDER] [--tmp-folder TMP_FOLDER]
                         [--analysis-folder ANALYSIS_FOLDER]
                         [--data-folder DATA_FOLDER] [--proj PROJ] [--wrap WRAP]
                         [--nDemes [NDEMES [NDEMES ...]]] [--diploid DIPLOID]
                         [--numMCMCIter NUMMCMCITER] [--numBurnIter NUMBURNITER]
                         [--numThinIter NUMTHINITER]

    optional arguments:
      -h, --help            show this help message and exit
      --dry                 if set, files are created, but eems is not ran
      --bed BED, --bfile BED
                            The common name of the bed file to be used
      --diffs DIFFS         A diffs file from a previous bed2diffs run, without
                            extension. If not None, the diffs file is filtered
                            rather than a new diffs file created from the bed
                            file.
      --loc LOC             File with location information. Should have a column
                            named `pop` and columns named `latitude` and
                            `longitude` or a variation thereof.
      --loc-has-no-header   add this flag if the location file (--loc) doesn't
                            have a header line. In this case, the first three
                            columns are assumed to be the population id, latitude
                            and longitude, respectively.
      --loc-has-header      add this flag if the location file (--loc) does have a
                            header line. In this case, the first three columns are
                            assumed to be the population id, latitude and
                            longitude, respectively.
      --sample-has-no-header
                            add this flag if the sample file (--sample) doesn't
                            have a header line. In this case, the first three
                            columns are assumed to be the sample id and population
                            id , respectively.
      --sample-has-header   add this flag if the sample file (--sample) does have
                            a header line. In this case, the first three columns
                            are assumed to be the sample id and population id,
                            latitude and longitude, respectively.
      --sample SAMPLE       File with sample information. Should have a column
                            named `sample` and columns named `pop` or a variation
                            thereof.
      --ind IND             File with individual based location information.
                            Should have a column named `sample` and columns named
                            `latitude` and `longitude` or a variation thereof.
      --ind-has-no-header   add this flag if the individual file (--ind) doesn't
                            have a header line. In this case, the first three
                            columns are assumed to be the individual id, latitude
                            and longitude, respectively.
      --ind-has-header      add this flag if the ind file (--ind) does have a
                            header line.
      --eems EEMS, --eems-folder EEMS
                            The base folder of your local eems installation
      --polygon POLYGON     A file with the polygon describing the region to run
                            eeems on.
      --bed2diffs BED2DIFFS
                            The bed2diffs executable, defaults to
                            bed2diffs/bed2diffs in the directory of --eems
      --eems_snps EEMS_SNPS
                            The eems_snps executable, defaults to
                            runeems_snps/runeems_snps in the directory of --eems
      --input-folder INPUT_FOLDER, --input_folder INPUT_FOLDER
                            the folder where all eems input files will be stored.
                            (default: ./input )
      --output-folder OUTPUT_FOLDER, --output_folder OUTPUT_FOLDER
                            the folder where all eems output files will be stored.
                            (default: ./output )
      --tmp-folder TMP_FOLDER, --tmp_folder TMP_FOLDER
                            the folder where tempory output files will be stored.
                            (default: ./tmp )
      --analysis-folder ANALYSIS_FOLDER, --analysis_folder ANALYSIS_FOLDER
                            the folder where all analysis will be stored (default:
                            . )
      --data-folder DATA_FOLDER, --data_folder DATA_FOLDER
                            the folder where all data files are read from
                            (default: . )
      --proj PROJ, --proj-name PROJ, --proj_name PROJ
                            The name of the output files
      --wrap WRAP, --wrap_america WRAP
                            Should all corinates be wrapped s.t. the americas
                            appear in the east?
      --nDemes [NDEMES [NDEMES ...]]
                            eems arg: number of demes in the model (default: 100 )
      --diploid DIPLOID     eems arg: is diploid? (true/false) (default: true )
      --numMCMCIter NUMMCMCITER, --n_mcmc NUMMCMCITER, --n-mcmc NUMMCMCITER
                            eems arg: number of MCMC iterations (default: 20000 )
      --numBurnIter NUMBURNITER, --n_burn NUMBURNITER, --n-burn NUMBURNITER
                            eems arg: number of burn in steps (default: 10000 )
      --numThinIter NUMTHINITER, --n_thin NUMTHINITER, --n-thin NUMTHINITER
                            eems arg: thinning interval (default: 99 )
