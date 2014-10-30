import re

import synapseclient
import numpy as np
import pandas as pd

import logging
logging.basicConfig(level=logging.DEBUG)

syn = synapseclient.login(silent=True)

_PCBC_IDS = (1773109, 2249208)

_PCBC_PRIVATE_ID = 2249208
_PCBC_PUBLIC_ID = 1773109

_MRNA_QUERY = "select sampleName,dataType,id,parentId,benefactorId,fileType,name from entity where fileType=='bam' AND dataType=='mRNA' AND bamType=='mapped'"

_MIRNA_QUERY = "select sampleName,dataType,id,parentId,benefactorId,fileType,name from entity where fileType=='fastq' AND dataType=='miRNA'"

# One file should not be included - it's the array description file: syn2654330
_METHYL_QUERY = "select sampleName,dataType,id,parentId,benefactorId,fileType,name from entity where parentId=='syn2653626'"

_METHYL_FILE_ID = "syn2677043"
_MRNA_MANUAL_META_ID = "syn2278178"
_MIRNA_MANUAL_META_ID = "syn2278179"

# possible differentation state suffixes
_DIFF_STATES = ('DE', 'EB', 'MESO', 'ECTO')

def fixRow(row):
    """Get the first item from column values with more than one.

    """

    for k, v in row.iteritems():
        if type(v) == list:
            row[k] = v[0]

    return row

def splitSampleName(s):
    """Split the sample name into separate components.

    These are a unique identifier, differentiation state, run, lane, and index.

    """

    (line, run, lane, index) = s.split(".")

    ## Remove and save the differentiation state
    diffState = re.findall("|".join(_DIFF_STATES), line)
    diffState = diffState[0] if len(diffState) > 0 else ''

    ## Remove and save the unique id
    uid = line
    uid = uid.split(diffState)[0] if diffState else uid
    uid = uid.split("_")[0]
    uid = uid.rstrip('-:')
    
    return pd.Series([uid, diffState, run, lane, index])

def getC4ID(name):
    """Clean up the sample name to get the C4 cell line ID.

    """

    c4id = None
    techrep = None
    biorep = None
    
    if name.startswith("SC"):
        # Matches a pattern SCNN-NNN
        c4idre = re.findall("SC[0-9][0-9]-[0-9][0-9][0-9][AB]{0,1}", name)
        c4id = c4idre[0]

        # Biological replicate possibilities
        if c4id.endswith("A") or c4id.endswith("B"):
            biorep = c4id[-1]
            c4id = c4id[:-1]
            
    elif name.startswith("H9"):
        c4id = "H9"
    elif name.startswith("IPS18"):
        c4id = "IPS18"
    else:
        c4id = name

    return c4id

def cleanColumnNames(x):
    if x.startswith("entity"):
        cleancol = x.split("entity.")[1]
    else:
        cleancol = x

    return cleancol

def getAnnotations(query):
    pcbcmeta_qry = syn.query(query)

    pcbcmeta = map(lambda x: fixRow(x),
                   pcbcmeta_qry['results'])
    
    pcbcmeta = pd.DataFrame(pcbcmeta)

    # Remove entity prefix from column names
    pcbcmeta.rename(columns=cleanColumnNames,
                    inplace=True)

    # Only those that are in the PCBC project
    pcbcmeta = pcbcmeta[pcbcmeta['benefactorId'].isin(_PCBC_IDS)]

    # # Split the sample name into separate columns
    # pcbc_sample_names = pcbcmeta['sampleName'].apply(splitSampleName)

    # pcbc_sample_names.rename(columns={0: 'UID',
    #                                   1: "Differentiation_State",
    #                                   2: 'Run', 3: 'Lane', 4: 'Index'},
    #                          inplace=True)

    # pcbc_sample_names['C4_Cell_Line_ID'] = pcbc_sample_names['UID'].apply(getC4ID)
    
    # pcbcmeta = pd.merge(pcbcmeta, pcbc_sample_names,
    #                left_index=True, right_index=True)

    return pcbcmeta


def main():
    import argparse
    import sys
    
    parser = argparse.ArgumentParser()

    parser.add_argument('outfile', nargs='?',
                        type=argparse.FileType('w'),
                        default=sys.stdout)
    
    args = parser.parse_args()

    # Fetch methyl data (no annotations)
    methyl_manual_meta = pd.read_csv(syn.get(_METHYL_FILE_ID).path, sep='\t')
    methyl_manual_meta['UID'] = methyl_manual_meta['Sample']
    methyl_manual_meta.rename(columns={"C4 Cell Line ID": "C4_Cell_Line_ID"},
                              inplace=True)
    methyl_manual_meta.set_index(['File'], inplace=True)

    
    # Fetch manually curated metadata
    mrna_manual_meta = pd.read_csv(syn.get(_MRNA_MANUAL_META_ID).path, sep='\t')
    mrna_manual_meta['UID'] = mrna_manual_meta['Decorated Name']
    mrna_manual_meta.set_index(['Decorated Name'], inplace=True)
    mrna_manual_meta.rename(columns={"CellLine": "C4_Cell_Line_ID"},
                            inplace=True)

    
    mirna_manual_meta = pd.read_csv(syn.get(_MIRNA_MANUAL_META_ID).path, sep='\t')
    mirna_manual_meta['UID'] = mirna_manual_meta['sample']
    mirna_manual_meta.set_index('sample', inplace=True)
    mirna_manual_meta.rename(columns={"C4 Cell Line ID": "C4_Cell_Line_ID"},
                             inplace=True)
    
    pcbc_mrna = getAnnotations(_MRNA_QUERY)
    pcbc_mrna.set_index(['sampleName'], inplace=True)

    pcbc_mrna = pd.merge(pcbc_mrna, mrna_manual_meta,
                         left_index=True, right_index=True)

    pcbc_sample_metadata = pcbc_mrna
    
    pcbc_mirna = getAnnotations(_MIRNA_QUERY)
    pcbc_mirna.set_index(['sampleName'], inplace=True)

    pcbc_mirna = pd.merge(pcbc_mirna, mirna_manual_meta,
                          left_index=True, right_index=True)

    pcbc_methyl = getAnnotations(_METHYL_QUERY)
    pcbc_methyl['File'] = pcbc_methyl['name'].apply(lambda x: x[:17])
    pcbc_methyl = pcbc_methyl.drop_duplicates()
    pcbc_methyl.set_index('File', inplace=True)

    pcbc_methyl = pd.merge(pcbc_methyl, methyl_manual_meta,
                          left_index=True, right_index=True)
    pcbc_methyl['dataType'] = "methylation"
        
    pcbc_sample_metadata = pd.concat([pcbc_mrna, pcbc_mirna, pcbc_methyl],
                                     ignore_index=True)

    # pcbc_sample_metadata.set_index(['sampleName', 'dataType'], inplace=True)
    keepcols = ["dataType", "UID", "C4_Cell_Line_ID", "Diffname short", "benefactorId", "fileType", "id", "parentId"]
    pcbc_sample_metadata[keepcols].to_csv(args.outfile, sep="\t", index=True)

if __name__ == "__main__":
    main()
