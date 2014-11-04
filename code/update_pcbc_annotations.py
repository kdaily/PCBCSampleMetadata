"""Use existing annotations to make a table of experimental metadata.

Works specifically for PCBC data - mRNA, miRNA, and methylation.

"""

import re
import sys
import tempfile

import pandas as pd
import synapseclient
import PCBCAnnotations
import synapseHelpers

import logging

def main():
    import argparse
    import sys
    
    parser = argparse.ArgumentParser()

    parser.add_argument('--dryrun', action="store_true",
                        default=False,
                        help="Run without making changes [default: %(default)s]")

    parser.add_argument('--storetable', action="store_true", default=False)
    
    parser.add_argument('--projectId', type=str)


    args = parser.parse_args()
    syn = synapseclient.login(silent=True)
    
    a = PCBCAnnotations.PCBCAnnotationTableUpdate(syn)

    my_synapse_id = None

    if not args.dryrun:
        my_synapse_id = synapseHelpers.thisCodeInSynapse(parentId='syn2758110', syn=syn)

    a.update_annots_synapse(executed=my_synapse_id, dryrun=args.dryrun)

    if args.storetable:
        tbl = a.update_annots_table_table(projectId=args.projectId, dryrun=args.dryrun)



if __name__ == "__main__":
    main()
