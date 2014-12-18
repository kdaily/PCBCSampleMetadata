"""Use existing annotations to make a table of experimental metadata.

Works specifically for PCBC data.

"""

import re
import sys
import tempfile

import pandas as pd
import synapseclient
import PCBCAnnotations

import logging

def main():
    import argparse
    import sys
    
    parser = argparse.ArgumentParser()

    # parser.add_argument('outfile', nargs='?',
    #                     type=argparse.FileType('w'),
    #                     default=sys.stdout)

    parser.add_argument('--empty',
                        type=str,
                        default=None)
    
    # parser.add_argument('--storetable',
    #                     action="store_true",
    #                     default=False)

    # parser.add_argument('--projectId', type=str)

    args = parser.parse_args()

    syn = synapseclient.login(silent=True)

    a = PCBCAnnotations.PCBCAnnotationChecker(syn)

    annots = a.update_annots_synapse()

    if args.empty:
        annots = annots[annots[args.empty] == ""]

    annots.to_csv("", sep="\t", index=False)

if __name__ == "__main__":
    main()
