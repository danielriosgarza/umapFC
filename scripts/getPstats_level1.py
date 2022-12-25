# -*- coding: utf-8 -*-
"""
Created on Sun Feb 13 12:39:14 2022

@author: u0139894
"""
import pandas as pd
import numpy as np
import os

#folder with the composition tables
dataFolder =  'G:\\My Drive\\FC\\LabelledFiles\\labelled_bhbtri'

#file with results
Results =  'G:\\My Drive\\FC\\results\\predicted_bhbtri_level1.csv'


predictions = os.listdir(dataFolder)


res = {}

for i,v in enumerate(predictions):
    if v.split('.')[-1] == 'csv':
        t = v.split('.')[0].split('_')
        
        if (len(t) == 4) and (t[0][0] == 't') and (t[1][0]=='d'):
            df = pd.read_csv(os.path.join(dataFolder, v), sep=',')
            res[i] = dict(time = int(t[0][1::]), dilution = int(t[1][1::]), condition = t[2], group=t[3], blanks = sum(df['state']=='bkgr'), live=sum(df['state']=='live'), inactive = sum(df['state']=='dead'), groupLabel = t[3][0:-1])
        
        

df = pd.DataFrame.from_dict(res).T

df.to_csv(Results)
            
        
            
        