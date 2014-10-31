import synapseclient
import synapseHelpers
import pandas as pd
import re

def getC4ID(name):
    """Clean up the sample name to get the C4 cell line ID.

    """

    c4id = None

    
    if name.startswith("SC"):
        # Matches a pattern SCNN-NNN
        c4id = re.findall("SC[0-9][0-9]-[0-9][0-9][0-9]", name)
        c4id = c4id[0]
    elif name.startswith("H9"):
        c4id = "H9"
    elif name.startswith("IPS18"):
        c4id = "IPS18"
    else:
        c4id = name

    return c4id

syn = synapseclient.Synapse(skip_checks=True)
syn.login(silent=True)

# possible differentation state suffixes
_DIFF_STATES = ('DE', 'EB', 'MESO', 'ECTO')

# Input from synapse
ARRAYMAP_ID = 'syn2654330'

# Output to synapse
_METHYL_PARENT_ID = "syn2775249"
_CODE_PARENT_ID = "syn2775246"

# Fetch array descriptions from Synapse 
arrayMap = pd.read_csv(syn.get(ARRAYMAP_ID).path, sep='\t')

# Rename columns as they would be from a query
arrayMap.rename(columns={"Sample": "sampleName"}, inplace=True)

# Add columns as they would be from a query
arrayMap["dataType"] = "methylation"
arrayMap["fileType"] = None

# Drop some columns; these should be stored in cell line
# level metadata
arrayMap = arrayMap.drop('Public?', 1)

diffStates = []
cellline = []

for i, s in enumerate(arrayMap.sampleName):

    ## Remove and save the differentiation state
    diffState = re.findall("|".join(_DIFF_STATES), s)
    diffState = diffState[0] if len(diffState) > 0 else ''
    diffStates.append(diffState)

    ## Remove and save the C4 cell line name
    cellid = s.split('_')[0].replace(diffState, '').rstrip('-:')
    cellline.append(cellid)

arrayMap['Differentiation_State'] = diffStates
arrayMap['C4_Cell_Line_ID'] = cellline
arrayMap.set_index('sampleName', inplace=True)

arrayMap.to_csv('Methylation_MetaData.tsv', sep='\t')

# f = syn.store(synapseclient.File('Methylation_MetaData.tsv',
#                                  parentId=_METHYL_PARENT_ID,
#                                  name='methylation metadata file'),
#               used = [ARRAYMAP_ID],
#               executed = synapseHelpers.thisCodeInSynapse(_CODE_PARENT_ID))
