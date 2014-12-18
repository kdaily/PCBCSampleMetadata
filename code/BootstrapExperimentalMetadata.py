"""Experimental metadata from Synapse files and annotations.

Currently specific for PCBC to assist with bootstrapping the Synapse file annotations.

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
    """Defines access to Synapse file annotations and manually created file metadata.

    The annotations come from a query defined by a class attribute (_QUERY).

    The metadata comes from a Synapse file defined by its ID (_MANUAL_META_ID).

    Only those that have a benefactor id in _LIMIT_IDS are kept (since query cannot have an OR statement)
    
    These are merged to create the minimum required metadata for PCBC files.
    
    """
    
    _LIMIT_IDS = (1773109, 2249208)
    
    _PCBC_PRIVATE_ID = 2249208
    _PCBC_PUBLIC_ID = 1773109
    
    # possible differentation state suffixes
    _DIFF_STATES = ('DE', 'EB', 'MESO', 'ECTO')
    
    def __init__(self, synclient):
        super(PCBCExperimentalMetadata, self).__init__(synclient)
        self.meta = None
    
    def __call__(self):
        
        self.meta = self.makeMetadata()
        
        return self.meta
    
    @staticmethod
    def fixRow(row):
        """Get the first item from column values with more than one.
        
        """
        
        for k, v in row.iteritems():
            if type(v) == list:
                row[k] = ",".join(v)
        
        return row
    
    @staticmethod
    def cleanColumnNames(x):
        if x.startswith("entity"):
            cleancol = x.split("entity.")[1]
        else:
            cleancol = x
        
        return cleancol
    
    def getAnnotations(self, query):
        """Run a query to get the annotations, clean up the rows.

        Only returns things that have a benefactorId in _LIMIT_IDS.
        
        Returns: pandas.DataFrame
        
        """
        
        # This is a large query, so need to chunk it
        # Make a temporary dict so it looks the same
        qr = self.syn.chunkedQuery(query)
        pcbcmeta = map(self.fixRow, qr)
        
        pcbcmeta = pd.DataFrame(pcbcmeta)
        
        # Remove entity prefix from column names
        pcbcmeta.rename(columns=self.cleanColumnNames,
                        inplace=True)
        
        # Only those that are in the PCBC project
        pcbcmeta = pcbcmeta[pcbcmeta['benefactorId'].isin(self._LIMIT_IDS)]
        
        return pcbcmeta

    def getManualMetadata(self):
        """Get the manually defined metadata from a synapse file defined by _MANUAL_META_ID

        Unique ID (UID) defined by _MANUAL_UID.
        Index for merging defined by _MANUAL_IDX.
        
        """
        
        manual_meta = pd.read_csv(self.syn.get(self._MANUAL_META_ID).path, sep='\t')
        manual_meta['UID'] = manual_meta[self._MANUAL_UID]
        manual_meta.set_index(self._MANUAL_IDX, inplace=True)
        
        manual_meta.rename(columns=self._COL_RENAME_DICT,
                                inplace=True)
        
        return manual_meta

    def getAnnotationMetadata(self):
        """Get existing annotation metadata from results of _QUERY. 

        Index set by _ANNOT_IDX.
        
        """
        
        annot_meta = self.getAnnotations(self._QUERY)
        annot_meta.set_index(self._ANNOT_IDX, inplace=True)
        
        return annot_meta
        
    def makeMetadata(self):
        """Get the manually defined metadata file and the annotations and join them.
        
        Merges based on _MANUAL_IDX and _ANNOT_IDX.
        
        Returns: pandas.DataFrame
        Sets: instance-level meta attribute.
        
        """
        
        manual_meta = self.getManualMetadata()
        annot_meta = self.getAnnotationMetadata()
        
        self.meta = pd.merge(annot_meta, manual_meta,
                             left_index=True, right_index=True)
        
        # Rename column names, Synapse tables don't play well with spaces
        self.meta.rename(columns=lambda x: x.replace(" ", "_"), inplace=True)
        
        return self.meta

    def to_csv(self, of):
        """Write merged annotation and manual metadata to a csv.
        
        If it hasn't been created, create it first.
        
        """
        
        if self.meta is None:
            self()
        
        keepcols = ["name", "id", "UID", "C4_Cell_Line_ID",
                    "Diffname_short", "dataType", "fileType"]
        
        self.meta[keepcols].to_csv(of, sep="\t", index=False)

class MRNAMetadata(PCBCExperimentalMetadata):
    """Default mRNA metadata using mapped bam files.
    
    """

    _QUERY = "select sampleName,dataType,id,benefactorId,fileType,name from entity where fileType=='bam' AND dataType=='mRNA' AND bamType=='mapped'"

    _MANUAL_META_ID = "syn2278178"

    # Set this as the UID column
    _MANUAL_UID = 'Decorated Name'
    
    # For indexing the table
    _MANUAL_IDX = "Decorated Name"
    
    _ANNOT_IDX = 'sampleName'
    _COL_RENAME_DICT = {"CellLine": "C4_Cell_Line_ID"}

class MRNABamMetadata(MRNAMetadata):
    """Default mRNA metadata using mapped bam files.
    
    """
    
    pass

class MRNABedMetadata(MRNAMetadata):
    
    _QUERY = "select sampleName,dataType,id,benefactorId,fileType,name from entity where fileType=='bed' AND dataType=='mRNA'"

class MRNAFastqMetadata(MRNAMetadata):
    
    _QUERY = "select sampleName,dataType,id,benefactorId,fileType,name from entity where fileType=='fastq' AND dataType=='mRNA'"

class MRNAFpkmMetadata(MRNAMetadata):
    
    _QUERY = "select sampleName,dataType,id,benefactorId,fileType,name from entity where fileType=='fpkm' AND dataType=='mRNA'"

class MRNAHTSeqCountMetadata(MRNAMetadata):
    
    _QUERY = "select sampleName,dataType,id,benefactorId,fileType,name from entity where fileType=='count' AND dataType=='mRNA'"

class MIRNAMetadata(PCBCExperimentalMetadata):
    """Default miRNA metadata.
    
    """
    
    _QUERY = "select sampleName,dataType,id,benefactorId,fileType,name from entity where fileType=='fastq' AND dataType=='miRNA'"
    
    _MANUAL_META_ID = "syn2278179"
    
    _MANUAL_UID = 'sample'
    _MANUAL_IDX = "sample"
    
    _ANNOT_IDX = 'sampleName'
    _COL_RENAME_DICT = {"C4 Cell Line ID": "C4_Cell_Line_ID"}

class MIRNAExprMetadata(PCBCExperimentalMetadata):
    _QUERY = "select sampleName,dataType,id,benefactorId,fileType,name from entity where fileType=='expr' AND dataType=='miRNA'"
    
    _MANUAL_META_ID = "syn2278179"
    
    _MANUAL_UID = 'sample'
    _MANUAL_IDX = "sample"
    
    _ANNOT_IDX = 'sampleName'
    _COL_RENAME_DICT = {"C4 Cell Line ID": "C4_Cell_Line_ID"}

class MIRNAFastqMetadata(PCBCExperimentalMetadata):
    
    _QUERY = "select sampleName,dataType,id,benefactorId,fileType,name from entity where fileType=='fastq' AND dataType=='miRNA'"
    
    _MANUAL_META_ID = "syn2278179"
    
    _MANUAL_UID = 'sample'
    _MANUAL_IDX = "sample"
    
    _ANNOT_IDX = 'sampleName'
    _COL_RENAME_DICT = {"C4 Cell Line ID": "C4_Cell_Line_ID"}

class MethylMetadata(PCBCExperimentalMetadata):
    """Class for methylation data.
    
    Overloads the makeMetadata() to add dataType and fileType annotations.
    
    """
    
    _QUERY = "select sampleName,dataType,id,benefactorId,fileType,name from entity where parentId=='syn2653626' AND id!='syn2654330'"
    
    _MANUAL_META_ID = "syn2677043"
    
    _MANUAL_UID = 'Sample'
    _MANUAL_IDX = "File"
    
    _ANNOT_IDX = 'File'
    _COL_RENAME_DICT = {"C4 Cell Line ID": "C4_Cell_Line_ID"}
    
    def getAnnotationMetadata(self):
        """Get existing annotation metadata from results of _QUERY. 
        
        Index set by _ANNOT_IDX.
        
        Overloaded for MethylMetadata to get an index to merge (File) from the name.
        
        """
        
        annot_meta = self.getAnnotations(self._QUERY)
        annot_meta['File'] = annot_meta['name'].apply(lambda x: x[:17])
        annot_meta.set_index(self._ANNOT_IDX, inplace=True)
        
        return annot_meta

def to_table(syn, meta, projectId):
    """Save to a synapse table.
    
    """

    keepcols = ["name", "id", "UID", "C4_Cell_Line_ID",
                "Diffname_short", "dataType", "fileType"]
    fortbl = meta[keepcols]
    
    cols = [table.Column(name='name', columnType='STRING', maximumSize=100),
            table.Column(name='id', columnType='ENTITYID'),
            # table.Column(name='benefactorId', columnType='STRING', maximumSize=100),
            table.Column(name='UID', columnType='STRING', maximumSize=100),
            table.Column(name='C4_Cell_Line_ID', columnType='STRING', maximumSize=100),
            table.Column(name='Diffname_short', columnType='STRING', maximumSize=100),
            table.Column(name='dataType', columnType='STRING', maximumSize=20),
            table.Column(name='fileType', columnType='STRING', maximumSize=20)]
        
    schema = table.Schema(name='PCBC Experimental Metadata', columns=cols, parent=projectId)

    my_table = table.create_table(schema, meta)
    logging.debug("schema: %s, table: %s" % (schema, my_table))
    my_table = syn.store(my_table)
    return my_table
    
