import synapseclient
import numpy as np
import pandas as pd

import logging
logging.basicConfig(level=logging.DEBUG)

myqy = "select sampleName,run,lane,dataType,id,parentId,benefactorId from entity where fileType=='bam' and dataType=='mRNA'"
