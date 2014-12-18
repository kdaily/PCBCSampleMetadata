"""Use a bootstrapped set of annotations to update all Synapse file annotations.

Currently for PCBC only.

This updates the experimental annotations for mRNA, miRNA, and methylation data.

These annotations should be provided by default in the future for any files added for mRNA,
miRNA, or methylation in the PCBC project.

"""

import re

# import pandas as pd
import csv

import synapseclient

import UpdateAnnotations
from UpdateAnnotations import UpdatePCBCAnnotations

def main():
    import argparse
    import sys
    
    parser = argparse.ArgumentParser()

    parser.add_argument('inputfile', nargs='?',
                        type=argparse.FileType('r'),
                        default=sys.stdout)

    parser.add_argument('--overwrite', action="store_true",
                        default=False,
                        help="Overwrite existing annotations [default: %(default)s]")
    
    parser.add_argument('--dryrun', action="store_true",
                        default=False,
                        help="Run without making changes [default: %(default)s]")

    parser.add_argument('--verbose', action="store_true",
                        default=False,
                        help="Output status [default: %(default)s]")

    args = parser.parse_args()

    syn = synapseclient.Synapse(debug=False, skip_checks=True)
    syn.login(silent=True)

    bootstrapped_data = list(csv.DictReader(args.inputfile, delimiter='\t'))
    
    updater = UpdatePCBCAnnotations(syn, bootstrapped_data)
    
    updater.update_annotations(overwrite=args.overwrite,
                               dryrun=args.dryrun, verbose=args.verbose)
    
if __name__ == "__main__":
    main()
