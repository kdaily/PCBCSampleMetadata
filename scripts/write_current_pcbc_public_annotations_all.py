"""Use existing annotations to make a table of experimental metadata.

Works specifically for PCBC data - mRNA, miRNA, and methylation.

"""

import re
import sys
import tempfile

import pandas as pd
import synapseclient
# from synapseclient import client
import PCBCAnnotations
import synapseHelpers

import logging

def main():
    import argparse
    import sys
    
    parser = argparse.ArgumentParser()

    parser.add_argument("--configPath", default=synapseclient.client.CONFIG_FILE)
    parser.add_argument('--dryrun', action="store_true",
                        default=False,
                        help="Run without making changes [default: %(default)s]")

    parser.add_argument('--parentId', type=str)
    parser.add_argument('--storetable', action="store_true", default=False)
    
    parser.add_argument('--projectId', type=str)


    args = parser.parse_args()

    syn=synapseclient.Synapse(skip_checks=True, configPath=args.configPath)
    syn.login(silent=True)
    
    annots = PCBCAnnotations.PCBCAllAnnotations(syn)
    a = PCBCAnnotations.PCBCAllAnnotationTableUpdate(syn, annots)
    
    my_synapse_id = None
    
    # if not args.dryrun:
    #     my_synapse_id = synapseHelpers.thisCodeInSynapse(parentId='syn2758110', syn=syn)
    
    a.update_annots_synapse(parentId=args.parentId, executed=my_synapse_id, dryrun=args.dryrun)
    
    if args.storetable:
        tbl = a.update_annots_table_synapse(projectId=args.projectId, dryrun=args.dryrun)

if __name__ == "__main__":
    main()
