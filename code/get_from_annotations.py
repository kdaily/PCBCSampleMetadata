import synapseclient
import numpy as np
import pandas as pd

import logging
# logging.basicConfig(level=logging.DEBUG)

syn = synapseclient.login()

_PCBC_IDS = (1773109, 2249208)

_PCBC_PRIVATE_ID = 2249208
_PCBC_PUBLIC_ID = 1773109

myqy = "select sampleName,run,lane,dataType,id,parentId,benefactorId,fileType,bamType from entity where fileType=='bam' and dataType=='mRNA'"

def fix_rows(row):
    """Get the first item from column values with more than one.

    """

    for k, v in row.iteritems():
        if type(v) == list:
            row[k] = v[0]

    return row

def main():
    import argparse

    parser = argparse.ArgumentParser()
    args = parser.parse_args()

    pcbcmeta_qry = syn.query(myqy)

    pcbcmeta = map(lambda x: fix_rows(x),
                   pcbcmeta_qry['results'])
    
    pcbcmeta = pd.DataFrame(pcbcmeta)
    print pcbcmeta.head()

if __name__ == "__main__":
    main()
