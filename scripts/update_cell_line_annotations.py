"""Add cell line annotations to files related to PCBC analyses from mRNA, miRNA, and methylation assays.

Requires that the files have an annotation called `C4_Cell_Line_ID` which will be used to merge and pull the cell line level annotations from a Synapse table.

"""

import tempfile
import csv

import pandas as pd
import synapseclient

import UpdateAnnotations
import PandasTable

def main():
    import argparse
    
    parser = argparse.ArgumentParser()

    parser.add_argument('--exptMetaId', type=str, required=True,
                        help="Synapse ID of experimental metadata table")

    parser.add_argument('--cellMetaId', type=str, required=True,
                        help="Synapse ID of cell line metadata table")
    
    parser.add_argument('--overwrite', action="store_true",
                        default=False,
                        help="Overwrite existing annotations [default: %(default)s]")

    parser.add_argument('--dryrun', action="store_true",
                        default=False,
                        help="Run without making changes [default: %(default)s]")

    parser.add_argument('-v', '--verbose', action="store_true",
                        default=False,
                        help="Verbose output [default: %(default)s]")
    
    parser.add_argument('--storeTableId', type=str, default=None,
                        help="Synapse project ID to store the data as a table.")

    args = parser.parse_args()
    
    syn = synapseclient.login(silent=True)
    
    cellline_metadata = PandasTable.PandasTable(syn=syn,
                                                table_id=args.cellMetaId)

    pcbc_experimental_metadata = PandasTable.PandasTable(syn=syn,
                                                         table_id=args.exptMetaId)

    # Remove anything without C4_Cell_Line_ID set
    pcbc_experimental_metadata = pcbc_experimental_metadata[pcbc_experimental_metadata["C4_Cell_Line_ID"].notnull()]

    merged_metadata = pd.merge(pcbc_experimental_metadata,
                               cellline_metadata,
                               on="C4_Cell_Line_ID", how="left")

    # Required for updater, expects list of dicts
    # Can't figure out how to properly do mapping on pandas.DataFrame
    tmpfile = tempfile.NamedTemporaryFile(delete=True)
    merged_metadata.to_csv(tmpfile.name, sep="\t", index=False)
    bootstrapped_data = list(csv.DictReader(tmpfile, delimiter='\t'))
    tmpfile.close()   
    
    updater = UpdateAnnotations.UpdatePCBCCellLineAnnotations(syn, bootstrapped_data)

    updater.update_annotations(overwrite=args.overwrite, dryrun=args.dryrun,
                               verbose=args.verbose)

    if args.storeTableId is not None and not args.dryrun:
        tbl = updater.update_annots_table_synapse(projectId=args.storeTableId,
                                                  dryrun=args.dryrun)

if __name__ == "__main__":
    main()
