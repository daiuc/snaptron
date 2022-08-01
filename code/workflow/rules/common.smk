import pandas as pa 
import os

# common variables - potentially move 
SampleLookupTable = pd.read_csv(config['RailID_Tissue_Lookup'], 
    sep='\t').set_index(keys=['rail_id', 'tissue', 'ds'], drop=False)

GTExTissues = SampleLookupTable.query('ds == "GTEx"').tissue.unique()
TcgaTissues = SampleLookupTable.query('ds == "TCGA"').tissue.unique()



