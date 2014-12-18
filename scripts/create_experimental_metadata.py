"""Use existing annotations and static files to make a table of experimental metadata.

Works specifically for PCBC data - mRNA, miRNA, and methylation.

"""


import pandas as pd
import synapseclient
import BootstrapExperimentalMetadata

import logging

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
    
    ## mRNA metadata from bam, bed, fastq, and fpkm files
    pcbcmeta = BootstrapExperimentalMetadata.MRNABamMetadata(syn)
    pcbc_mrna_bam = pcbcmeta()

    sys.stderr.write("Getting mRNA fastq data\n")
    pcbcmeta = BootstrapExperimentalMetadata.MRNAFastqMetadata(syn)
    pcbc_mrna_fastq = pcbcmeta()

    pcbcmeta = BootstrapExperimentalMetadata.MRNABedMetadata(syn)
    pcbc_mrna_bed = pcbcmeta()
    
    pcbcmeta = BootstrapExperimentalMetadata.MRNAFpkmMetadata(syn)
    pcbc_mrna_fpkm = pcbcmeta()
    
    ## miRNA metadata from Fastq and Expr files
    pcbcmeta = BootstrapExperimentalMetadata.MIRNAFastqMetadata(syn)
    pcbc_mirna_fastq = pcbcmeta()
    
    pcbcmeta = BootstrapExperimentalMetadata.MIRNAExprMetadata(syn)
    pcbc_mirna_expr = pcbcmeta()
    
    pcbcmeta = BootstrapExperimentalMetadata.MethylMetadata(syn)
    pcbc_methyl = pcbcmeta()

    pcbc_sample_metadata = pd.concat([pcbc_mrna_bam, pcbc_mrna_bed,
                                      pcbc_mrna_fpkm,
                                      pcbc_mrna_fastq,
                                      pcbc_mirna_fastq, pcbc_mirna_expr,
                                      pcbc_methyl],
                                     ignore_index=True)

    _keepcols = ["name", "id", "UID", "C4_Cell_Line_ID",
                 "Diffname_short", "dataType", "fileType"]

    pcbc_sample_metadata = pcbc_sample_metadata[_keepcols]

    pcbc_sample_metadata.to_csv(args.outfile, sep="\t", index=False)
    
    if args.storetable:
        tbl = BootstrapExperimentalMetadata.to_table(syn, pcbc_sample_metadata,
                                                     args.projectId)

    # pcbcmeta.to_csv(args.outfile)
        
if __name__ == "__main__":
    main()
