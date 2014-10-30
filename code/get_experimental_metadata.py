import re

import pandas as pd
import synapseclient
import ExperimentalMetadata

import logging
logging.basicConfig(level=logging.DEBUG)


def main():
    import argparse
    import sys
    
    parser = argparse.ArgumentParser()

    parser.add_argument('outfile', nargs='?',
                        type=argparse.FileType('w'),
                        default=sys.stdout)

    parser.add_argument('--storetable',
                        action="store_true",
                        default=False)

    parser.add_argument('--projectId', type=str)

    args = parser.parse_args()

    syn = synapseclient.login(silent=True)

    # pcbcmeta = ExperimentalMetadata.AllExperimentalMetadata(syn)
    # pcbc_sample_metadata = pcbcmeta()

    pcbcmeta = ExperimentalMetadata.MRNAMetadata(syn)
    pcbc_mrna = pcbcmeta()

    pcbcmeta = ExperimentalMetadata.MIRNAMetadata(syn)
    pcbc_mirna = pcbcmeta()

    pcbcmeta = ExperimentalMetadata.MethylMetadata(syn)
    pcbc_methyl = pcbcmeta()

    pcbc_sample_metadata = pd.concat([pcbc_mrna, pcbc_mirna, pcbc_methyl],
                                     ignore_index=True)

    _keepcols = ["name", "id", "UID", "C4_Cell_Line_ID",
                 "Diffname short", "dataType", "fileType"]

    pcbc_sample_metadata = pcbc_sample_metadata[_keepcols]

    pcbc_sample_metadata.to_csv(args.outfile, sep="\t", index=False)
    
    if args.storetable:
        tbl = ExperimentalMetadata.to_table(syn, pcbc_sample_metadata,
                                            args.projectId)

    # pcbcmeta.to_csv(args.outfile)
        
if __name__ == "__main__":
    main()
