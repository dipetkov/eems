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

from pipeline import run_cli

if __name__ == "__main__":
    run_cli()
