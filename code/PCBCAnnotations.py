import logging
import sys
import csv
import itertools

import pandas as pd
import synapseclient
from synapseclient import table

import synapseHelpers

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)

_PCBC_PRIVATE_ID = 'syn2249208'
_PCBC_PUBLIC_ID = 'syn1773109'

## Generic queries to get various things
_COLS = ["id", "name", "parentId", "benefactorId",
         "UID", "C4_Cell_Line_ID", "Diffname_short",
         "dataType", "fileType"
     ]
_PRIVATE_FILES_QUERY = "select %(cols)s from file where benefactorId=='%(benefactorId)s' AND parentId=='%(parentId)s'" % dict(cols=",".join(_COLS), benefactorId=_PCBC_PRIVATE_ID, parentId='syn2261941')

_PUBLIC_FILES_QUERY = "select %(cols)s from file where benefactorId=='%(benefactorId)s'" % dict(cols=",".join(_COLS), benefactorId=_PCBC_PUBLIC_ID)

def fixRow(row):
    """Get the first item from column values with more than one.
    
    """
    
    for k, v in row.iteritems():
        if type(v) == list:
            row[k] = ",".join(v)
            
    return row

class ExperimentalAnnotations(object):
    def __init__(self, synclient):
        self.syn = synclient

class PCBCAnnotations(ExperimentalAnnotations):
    """PCBC public file anotations.
    
    """
    
    _COLS = ["id", "name", "parentId", "benefactorId",
             "UID", "C4_Cell_Line_ID", "Diffname_short",
             "dataType", "fileType"
    ]
    
    _INPUT_BENEFACTOR_ID = _PCBC_PUBLIC_ID
    
    @staticmethod
    def fixRow(row):
        """Get the first item from column values with more than one.
    
        """
        
        for k, v in row.iteritems():
            if type(v) == list:
                row[k] = ",".join(v)
        
        return row
    
    def get_annotations(self):
        
        _q_params = dict(cols=",".join(self._COLS),
                         id=self._INPUT_BENEFACTOR_ID)
        
        _query = "select %(cols)s from file where benefactorId=='%(id)s'" % _q_params
        
        logger.debug(_query)
        
        # This is a large query, so need to chunk it
        # Make a temporary dict so it looks the same
        qr = self.syn.chunkedQuery(_query)
        annot_qry = dict(results=map(self.fixRow, qr))
        
        self.annots = synapseHelpers.query2df(annot_qry, filterSynapseFields=False)
        
        return self.annots

class PCBCPublicAnnotations(PCBCAnnotations):
    pass
    
class PCBCPrivateAnnotations(PCBCAnnotations):
    _INPUT_BENEFACTOR_ID = _PCBC_PRIVATE_ID

class PCBCAnnotationTableUpdate(object):
    """Make a synapse file that has up-to-date sample metadata annotation information for all files.
    
    """
    
    _FILENAME = "/tmp/pcbcannotation_public_status.tsv"
    _OUTPUT_PARENT_ID = 'syn2805609'
    
    _INPUT_PARENT_IDS = ('syn2653626', # methylation
                         'syn2247164', #miRNA miExpress
                         'syn2247098', #miRNA fastq
                         'syn1773111', # mRNA bam
                         'syn2246521', #mRNA expression
                         'syn1773112' # mRNA fastq
                        )
    
    _OUTPUT_COLS = ["id", "name", "UID", "C4_Cell_Line_ID",
                    "Diffname_short", "dataType", "fileType"]

    _TABLE_NAME = 'PCBC Public Experimental Metadata'

    _TABLE_COLS = [table.Column(name='name', columnType='STRING', maximumSize=100),
                   table.Column(name='id', columnType='ENTITYID'),
                   table.Column(name='UID', columnType='STRING', maximumSize=100),
                   table.Column(name='C4_Cell_Line_ID', columnType='STRING', maximumSize=100),
                   table.Column(name='Diffname_short', columnType='STRING', maximumSize=100),
                   table.Column(name='dataType', columnType='STRING', maximumSize=100),
                   table.Column(name='fileType', columnType='STRING', maximumSize=100)]

    def __init__(self, syn, annotations):
        self.syn = syn
        self.annotations = annotations

    def _get_annotations(self):
        """Convenience function to get annotations and filter based on parent id, and select columns.
        
        """

        self.curr_annots = self.annotations.get_annotations()
        self.curr_annots = self.curr_annots[self.curr_annots["parentId"].isin(self._INPUT_PARENT_IDS)]
        self.curr_annots = self.curr_annots[self._OUTPUT_COLS]
        
    def update_annots_synapse(self, executed=None, dryrun=False):
        self._get_annotations()
        
        self.curr_annots.to_csv(self._FILENAME, sep='\t', index=False)
        myfile = synapseclient.File(self._FILENAME, parentId=self._OUTPUT_PARENT_ID)

        if not dryrun:
            myfile = self.syn.store(myfile, executed=executed)
            logger.info("file ID: %s, filename: %s" % (myfile.id, myfile.name))
        else:
            logger.info("DRYRUN: filename: %s" % (myfile.name))

    def update_annots_table_synapse(self, projectId=_PCBC_PUBLIC_ID, dryrun=False):

        self._get_annotations()
                
        schema = table.Schema(name=self._TABLE_NAME,
                              columns=self._TABLE_COLS,
                              parent=projectId)

        my_table = table.create_table(schema, self.curr_annots)
        logger.debug("schema: %s, table: %s" % (schema, my_table))

        if not dryrun:
            my_table = self.syn.store(my_table)

        return my_table


class PCBCAllAnnotations(PCBCAnnotations):
    """Get experimental and cell line annotations for all PCBC files.
    
    """

    _FILENAME = "/tmp/pcbcannotation_public_status_all.tsv"
    _OUTPUT_PARENT_ID = 'syn2805609'

    _COLS = ["id", "name", "parentId", "benefactorId", "PCBC_Cell_Line_Name",
             "Host_Species", "UID", "C4_Cell_Line_ID", "Diffname_short",
             "dataType", "fileType", "Originating_Lab_ID",
             "Public_Data", "Cell_Type", "Cell_Line_Type",
             "Cell_Type_of_Origin", "Cell_Line_of_Origin", "Tissue_of_Origin",
             "Reprogramming_Vector_Type", "Reprogramming_Gene_Combination",
             "Culture_Conditions", "Small_Molecules",
             "Other_Conditions_During_Reprogramming", "Donor_Life_Stage", "Race",
             "Ethnicity", "Gender", "Disease", "Donor_Phenotype", "Genotype",
             "C4_Karyotype_Result", "Originating_Lab",
             "Other_Significant_Contributors_to_line_generation",
             "Pubmed_ID", "DOI", "Donor_Cell_Link", "High_Confidence_Donor_ID"
         ]

class PCBCAllAnnotationTableUpdate(PCBCAnnotationTableUpdate):
    """Make a synapse file that has up-to-date sample and cell line metadata annotation information for all files.
    """

    _TABLE_NAME = 'PCBC Public All Metadata'

    _OUTPUT_COLS = ["id", "name", "PCBC_Cell_Line_Name",
             "Host_Species", "UID", "C4_Cell_Line_ID", "Diffname_short",
             "dataType", "fileType", "Originating_Lab_ID",
             "Public_Data", "Cell_Type", "Cell_Line_Type",
             "Cell_Type_of_Origin", "Cell_Line_of_Origin", "Tissue_of_Origin",
             "Reprogramming_Vector_Type", "Reprogramming_Gene_Combination",
             "Culture_Conditions", "Small_Molecules",
             "Other_Conditions_During_Reprogramming", "Donor_Life_Stage", "Race",
             "Ethnicity", "Gender", "Disease", "Donor_Phenotype", "Genotype",
             "C4_Karyotype_Result", "Originating_Lab",
             "Other_Significant_Contributors_to_line_generation",
             "Pubmed_ID", "DOI", "Donor_Cell_Link", "High_Confidence_Donor_ID"
         ]

    _TABLE_COLS = [table.Column(name='name', columnType='STRING', maximumSize=100),
                   table.Column(name='id', columnType='ENTITYID'),
                   table.Column(name='UID', columnType='STRING', maximumSize=100),
                   table.Column(name='C4_Cell_Line_ID', columnType='STRING',
                                maximumSize=100),
                   table.Column(name='Diffname_short', columnType='STRING',
                                maximumSize=100),
                   table.Column(name='dataType', columnType='STRING', maximumSize=100),
                   table.Column(name='fileType', columnType='STRING', maximumSize=100),
                   table.Column(name="PCBC_Cell_Line_Name", columnType='STRING',
                                maximumSize=100),
                   table.Column(name="Host_Species", columnType='STRING',
                                maximumSize=100),
                   table.Column(name="Originating_Lab_ID", columnType='STRING',
                                maximumSize=100),
                   table.Column(name="Public_Data", columnType='STRING',
                                maximumSize=100),
                   table.Column(name="Cell_Type", columnType='STRING',
                                maximumSize=100),
                   table.Column(name="Cell_Line_Type", columnType='STRING',
                                maximumSize=100),
                   table.Column(name="Cell_Type_of_Origin", columnType='STRING',
                                maximumSize=100),
                   table.Column(name="Cell_Line_of_Origin", columnType='STRING',
                                maximumSize=100),
                   table.Column(name="Tissue_of_Origin", columnType='STRING',
                                maximumSize=500),
                   table.Column(name="Reprogramming_Vector_Type", columnType='STRING',
                                maximumSize=500),
                   table.Column(name="Reprogramming_Gene_Combination",
                                columnType='STRING', maximumSize=500),
                   table.Column(name="Culture_Conditions", columnType='STRING',
                                maximumSize=500),
                   table.Column(name="Small_Molecules", columnType='STRING',
                                maximumSize=500),
                   table.Column(name="Other_Conditions_During_Reprogramming",
                                columnType='STRING', maximumSize=500),
                   table.Column(name="Donor_Life_Stage", columnType='STRING',
                                maximumSize=100),
                   table.Column(name="Race", columnType='STRING', maximumSize=100),
                   table.Column(name="Ethnicity", columnType='STRING',
                                maximumSize=100),
                   table.Column(name="Gender", columnType='STRING', maximumSize=100),
                   table.Column(name="Disease", columnType='STRING', maximumSize=500),
                   table.Column(name="Donor_Phenotype", columnType='STRING',
                                maximumSize=500),
                   table.Column(name="Genotype", columnType='STRING', maximumSize=500),
                   table.Column(name="C4_Karyotype_Result", columnType='STRING',
                                maximumSize=500),
                   table.Column(name="Originating_Lab", columnType='STRING',
                                maximumSize=500),
                   table.Column(name="Other_Significant_Contributors_to_line_generation",
                                columnType='STRING', maximumSize=500),
                   table.Column(name="Pubmed_ID", columnType='STRING',
                                maximumSize=100),
                   table.Column(name="DOI", columnType='STRING', maximumSize=100),
                   table.Column(name="Donor_Cell_Link", columnType='STRING',
                                maximumSize=500),
                   table.Column(name="High_Confidence_Donor_ID", columnType='STRING',
                                maximumSize=100)
    ]
