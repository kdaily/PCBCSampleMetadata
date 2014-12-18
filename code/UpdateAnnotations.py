import logging
import sys

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)

class ExperimentalAnnotations(object):
    def __init__(self, synclient):
        self.syn = synclient

# Which column is the UID?
_UID_LOOKUP = dict(mRNA='sampleName', miRNA='sampleName')

class UpdatePCBCAnnotations(ExperimentalAnnotations):
    """Update the PCBC annotations.

    """
    
    # These are the annots that we'll be adding
    _ANNOTS_TO_ADD = ("UID", "C4_Cell_Line_ID", "Diffname_short", "dataType", "fileType")

    def __init__(self, synclient, bootstrapped_data):
        super(UpdatePCBCAnnotations, self).__init__(synclient)
        self.bootstrapped_data = bootstrapped_data

    @staticmethod
    def _test(synobj, record, objannot, recannot):
        """Compare two values from the record and an existing synapse object.
        
        """

        if type(synobj[objannot]) == list:
            assert synobj[objannot][0] == record[recannot], \
                "Mismatch between %s and %s (synobj:%s != record:%s)" % (objannot,
                                                                        recannot,
                                                                        synobj[objannot],
                                                                        record[recannot])
        else:
            assert synobj[objannot] == record[recannot], \
                "Mismatch between %s and %s (synobj:%s != record:%s)" % (objannot,
                                                                         recannot,
                                                                         synobj[objannot],
                                                                         record[recannot])

    def test_annotations(self, synobj, record):
        """Tests that should pass before adding annotations.
        """
        
        self._test(synobj, record, "name", "name")
        
        try:
            self._test(synobj, record,
                       _UID_LOOKUP[synobj.dataType[0]], 'UID')
        except KeyError:
            pass

    def _update_annotation(self, record, overwrite=False, dryrun=False,
                           ignoreerrors=False, verbose=False):
        """Update annotations on a synapse object.

        If overwrite is True and an annotation already exists, this will replace it.

        If ignoreerrors is True, will continue.
        
        """

        # Get the object
        synobj = self.syn.get(record['id'], downloadFile=False)

        # Test to make sure this the right thing to update
        self.test_annotations(synobj, record)

        # Do the updates
        for annot_to_add in self._ANNOTS_TO_ADD:
            
            try:
                # Check if it already has the annotation
                if synobj[annot_to_add]:
                    if overwrite:
                        synobj[annot_to_add] = record[annot_to_add]
                    else:
                        if ignoreerrors and verbose:
                            logger.info("Annotation '%s' exists. Continuing." % (annot_to_add,))
                        else:
                            raise KeyError, \
                                "Annotation '%s' already exists and overwrite is False." % (annot_to_add)
            except KeyError:
                synobj[annot_to_add] = record[annot_to_add]

        if not dryrun:
            if verbose:
                logger.info("Saving annotations for %s..." % (synobj.id,))
            
            synobj = self.syn.store(synobj, forceVersion=False)
        else:
            if verbose:
                logger.info("Dry run, not saving annotations for %s..." % (synobj.id,))

    
    def update_annotations(self, overwrite=False, dryrun=False, verbose=False):
        """Run the update.

        """


        myfun = lambda x: self._update_annotation(x,
                                                  overwrite=overwrite,
                                                  dryrun=dryrun,
                                                  verbose=verbose)

        # No output
        res = map(myfun, self.bootstrapped_data)
        
        
class UpdatePCBCCellLineAnnotations(UpdatePCBCAnnotations):
    """Update the PCBC annotations with cell line information.

    """
    
    # These are the annots that we'll be adding
    _ANNOTS_TO_ADD = ("Originating_Lab_ID", "Public_Data", "Cell_Type", "Cell_Line_Type",
                      "Cell_Type_of_Origin", "Cell_Line_of_Origin", "Tissue_of_Origin",
                      "Reprogramming_Vector_Type", "Reprogramming_Gene_Combination",
                      "Culture_Conditions", "Small_Molecules",
                      "Other_Conditions_During_Reprogramming", "Donor_Life_Stage", "Race",
                      "Ethnicity", "Gender", "Disease", "Donor_Phenotype", "Genotype",
                      "C4_Karyotype_Result", "Originating_Lab",
                      "Other_Significant_Contributors_to_line_generation",
                      "Pubmed_ID", "DOI", "Donor_Cell_Link", "High_Confidence_Donor_ID")
