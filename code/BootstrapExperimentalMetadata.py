"""Experimental metadata from Synapse files and annotations.

Currently specific for PCBC to assist with bootstrapping the annotations,
and enforce that particular 
Kenneth Daily, 2014

"""

import pandas as pd

from synapseclient import Project, File, Folder, Schema
from synapseclient import table


import logging
logging.basicConfig(level=logging.DEBUG)

class ExperimentalMetadata(object):
    def __init__(self, synclient):
        self.syn = synclient

class PCBCExperimentalMetadata(ExperimentalMetadata):
    _LIMIT_IDS = (1773109, 2249208)

    _PCBC_PRIVATE_ID = 2249208
    _PCBC_PUBLIC_ID = 1773109
    
    # possible differentation state suffixes
    _DIFF_STATES = ('DE', 'EB', 'MESO', 'ECTO')

    def __init__(self, synclient):
        ExperimentalMetadata.__init__(self, synclient)
        self.meta = None

    def __call__(self):

        self.meta = self.getAnnotations()
        
        return self.meta
    
    @staticmethod
    def fixRow(row):
        """Get the first item from column values with more than one.
    
        """
    
        for k, v in row.iteritems():
            if type(v) == list:
                row[k] = v[0]
    
        return row

    @staticmethod
    def cleanColumnNames(x):
        if x.startswith("entity"):
            cleancol = x.split("entity.")[1]
        else:
            cleancol = x

        return cleancol

    def _getAnnotations(self, query):
        pcbcmeta_qry = self.syn.query(query)
    
        pcbcmeta = map(lambda x: self.fixRow(x),
                    pcbcmeta_qry['results'])
        
        pcbcmeta = pd.DataFrame(pcbcmeta)
    
        # Remove entity prefix from column names
        pcbcmeta.rename(columns=self.cleanColumnNames,
                        inplace=True)

        # Only those that are in the PCBC project
        pcbcmeta = pcbcmeta[pcbcmeta['benefactorId'].isin(self._LIMIT_IDS)]

        return pcbcmeta

    def getAnnotations(self):
        manual_meta = pd.read_csv(self.syn.get(self._MANUAL_META_ID).path, sep='\t')
        manual_meta['UID'] = manual_meta[self._MANUAL_UID]
        manual_meta.set_index(self._MANUAL_IDX, inplace=True)
        
        manual_meta.rename(columns=self._COL_RENAME_DICT,
                                inplace=True)
        
        annot_meta = self._getAnnotations(self._QUERY)
        annot_meta.set_index(self._ANNOT_IDX, inplace=True)
        
        self.meta = pd.merge(annot_meta, manual_meta,
                             left_index=True, right_index=True)
        
        return self.meta

    def to_csv(self, of):

        if self.meta is None:
            self()

        keepcols = ["name", "id", "benefactorId", "UID", "C4_Cell_Line_ID",
                    "Diffname short", "dataType", "fileType"]
        self.meta[keepcols].to_csv(of, sep="\t", index=False)

class MRNAMetadata(PCBCExperimentalMetadata):
    _QUERY = "select sampleName,dataType,id,benefactorId,fileType,name from entity where fileType=='bam' AND dataType=='mRNA' AND bamType=='mapped'"
    _MANUAL_META_ID = "syn2278178"

    # Set this as the UID column
    _MANUAL_UID = 'Decorated Name'

    # For indexing the table
    _MANUAL_IDX = "Decorated Name"

    _ANNOT_IDX = 'sampleName'
    _COL_RENAME_DICT = {"CellLine": "C4_Cell_Line_ID"}

    def __call__(self):

        self.meta = self.getAnnotations()
        
        return self.meta

class MIRNAMetadata(PCBCExperimentalMetadata):
    _QUERY = "select sampleName,dataType,id,benefactorId,fileType,name from entity where fileType=='fastq' AND dataType=='miRNA'"
    _MANUAL_META_ID = "syn2278179"

    _MANUAL_UID = 'sample'
    _MANUAL_IDX = "sample"

    _ANNOT_IDX = 'sampleName'
    _COL_RENAME_DICT = {"C4 Cell Line ID": "C4_Cell_Line_ID"}
    

class MethylMetadata(PCBCExperimentalMetadata):
    _QUERY = "select sampleName,dataType,id,benefactorId,fileType,name from entity where parentId=='syn2653626' AND id!='syn2654330'"
    
    _METHYL_FILE_ID = "syn2677043"

    _MANUAL_UID = 'Sample'
    _MANUAL_IDX = "File"

    _ANNOT_IDX = 'File'
    _COL_RENAME_DICT = {"C4 Cell Line ID": "C4_Cell_Line_ID"}

    # Overloads getAnnotations()
    def getAnnotations(self):
        manual_meta = pd.read_csv(self.syn.get(self._METHYL_FILE_ID).path, sep='\t')
        manual_meta['UID'] = manual_meta[self._MANUAL_UID]

        manual_meta.rename(columns=self._COL_RENAME_DICT,
                           inplace=True)

        manual_meta.set_index(self._MANUAL_IDX, inplace=True)
        
        annot_meta = self._getAnnotations(self._QUERY)
        annot_meta['File'] = annot_meta['name'].apply(lambda x: x[:17])
        # annot_meta = annot_meta.drop_duplicates()
        
        annot_meta.set_index(self._ANNOT_IDX, inplace=True)
    
        annot_meta = pd.merge(annot_meta, manual_meta,
                            left_index=True, right_index=True)

        # Add some annotations that aren't there now'
        annot_meta['dataType'] = "methylation"
        annot_meta['fileType'] = "idat"
    
        return annot_meta

def to_table(syn, meta, projectId):

    keepcols = ["name", "id", "UID", "C4_Cell_Line_ID",
                "Diffname short", "dataType", "fileType"]
    fortbl = meta[keepcols]
    
    cols = [table.Column(name='name', columnType='STRING', maximumSize=100),
            table.Column(name='id', columnType='ENTITYID'),
            # table.Column(name='benefactorId', columnType='STRING', maximumSize=100),
            table.Column(name='UID', columnType='STRING', maximumSize=100),
            table.Column(name='C4_Cell_Line_ID', columnType='STRING', maximumSize=100),
            table.Column(name='Diffname short', columnType='STRING', maximumSize=100),
            table.Column(name='dataType', columnType='STRING', maximumSize=20),
            table.Column(name='fileType', columnType='STRING', maximumSize=20)]
        
    schema = table.Schema(name='PCBC Experimental Metadata', columns=cols, parent=projectId)

    my_table = table.create_table(schema, meta)
    logging.debug("schema: %s, table: %s" % (schema, my_table))
    my_table = syn.store(my_table)
    return my_table
    
class AllExperimentalMetadata(PCBCExperimentalMetadata):

    _LIMIT_IDS = (1773109, 2249208)

    _PCBC_PRIVATE_ID = 2249208
    _PCBC_PUBLIC_ID = 1773109
    
    _MRNA_QUERY = "select sampleName,dataType,id,benefactorId,fileType,name from entity where fileType=='bam' AND dataType=='mRNA' AND bamType=='mapped'"
    
    _MIRNA_QUERY = "select sampleName,dataType,id,benefactorId,fileType,name from entity where fileType=='fastq' AND dataType=='miRNA'"
    
    # One file should not be included - it's the array description file: syn2654330
    _METHYL_QUERY = "select sampleName,dataType,id,benefactorId,fileType,name from entity where parentId=='syn2653626' AND id!='syn2654330'"
    
    _METHYL_FILE_ID = "syn2677043"
    _MRNA_MANUAL_META_ID = "syn2278178"
    _MIRNA_MANUAL_META_ID = "syn2278179"
    
    # possible differentation state suffixes
    _DIFF_STATES = ('DE', 'EB', 'MESO', 'ECTO')

    def __init__(self, synclient):
        ExperimentalMetadata.__init__(self, synclient)
        self.pcbcmeta = None
        
    def __call__(self):

        pcbc_mrna = self.getMRNAAnnotations()
        pcbc_mirna = self.getMIRNAAnnotations()
        pcbc_methyl = self.getMethylAnnotations()
        
        pcbc_sample_metadata = pd.concat([pcbc_mrna, pcbc_mirna, pcbc_methyl],
                                         ignore_index=True)

        self.pcbcmeta = pcbc_sample_metadata
        return self.pcbcmeta


    def getMethylAnnotations(self):
        # Fetch manually curated methylation data
        methyl_manual_meta = pd.read_csv(self.syn.get(self._METHYL_FILE_ID).path, sep='\t')
        methyl_manual_meta['UID'] = methyl_manual_meta['Sample']
        methyl_manual_meta.rename(columns={"C4 Cell Line ID": "C4_Cell_Line_ID"},
                                inplace=True)
        methyl_manual_meta.set_index(['File'], inplace=True)
    
        pcbc_methyl = self._getAnnotations(self._METHYL_QUERY)
        pcbc_methyl['File'] = pcbc_methyl['name'].apply(lambda x: x[:17])
        pcbc_methyl = pcbc_methyl.drop_duplicates()
        pcbc_methyl.set_index('File', inplace=True)
    
        pcbc_methyl = pd.merge(pcbc_methyl, methyl_manual_meta,
                            left_index=True, right_index=True)
        pcbc_methyl['dataType'] = "methylation"
        pcbc_methyl['fileType'] = "idat"
    
        return pcbc_methyl
    
    def getMRNAAnnotations(self):
        # Fetch manually curated metadata
        mrna_manual_meta = pd.read_csv(self.syn.get(self._MRNA_MANUAL_META_ID).path, sep='\t')
        mrna_manual_meta['UID'] = mrna_manual_meta['Decorated Name']
        mrna_manual_meta.set_index(['Decorated Name'], inplace=True)
        mrna_manual_meta.rename(columns={"CellLine": "C4_Cell_Line_ID"},
                                inplace=True)
    
        pcbc_mrna = self._getAnnotations(self._MRNA_QUERY)
        pcbc_mrna.set_index(['sampleName'], inplace=True)
    
        pcbc_mrna = pd.merge(pcbc_mrna, mrna_manual_meta,
                            left_index=True, right_index=True)
    
        return pcbc_mrna
    
    def getMIRNAAnnotations(self):
        mirna_manual_meta = pd.read_csv(self.syn.get(self._MIRNA_MANUAL_META_ID).path, sep='\t')
        mirna_manual_meta['UID'] = mirna_manual_meta['sample']
        mirna_manual_meta.set_index('sample', inplace=True)
        mirna_manual_meta.rename(columns={"C4 Cell Line ID": "C4_Cell_Line_ID"},
                                inplace=True)
        
        pcbc_mirna = self._getAnnotations(self._MIRNA_QUERY)
        pcbc_mirna.set_index(['sampleName'], inplace=True)
        
        pcbc_mirna = pd.merge(pcbc_mirna, mirna_manual_meta,
                            left_index=True, right_index=True)
        
        return pcbc_mirna
