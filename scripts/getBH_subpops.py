# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 04:21:37 2022

@author: u0139894
"""

import pandas as pd
import numpy as np
import os
from sklearn.mixture import GaussianMixture

#folder with the composition tables
dataFolder1 =  'G:\\My Drive\\FC\\LabelledFiles\\labelled_bhbt'
dataFolder2 =  'G:\\My Drive\\FC\\LabelledFiles\\labelled_bhri'
dataFolder3 = 'G:\\My Drive\\FC\\LabelledFiles\\labelled_bhbtri'

#file with results
Results =  'G:\\My Drive\\FC\\results\\predicted_bh_groups_bhbtri.csv'


predictions = [dataFolder1 + '\\' + i for i in  os.listdir(dataFolder1)] + [dataFolder2 + '\\' + i for i in  os.listdir(dataFolder2)]   + [dataFolder3 + '\\' + i for i in  os.listdir(dataFolder3)]


liveEvents = np.array([])

for i,v in enumerate(predictions):
    
    if v.split('.')[-1] == 'csv':
        t = v.split('\\')[-1].split('.')[0].split('_')
        
        if (len(t) == 4) and (t[0][0] == 't') and (t[1][0]=='d'):
            
            if ((t[3]=='bhA') or (t[3]=='bhB') or (t[3]=='bhC')) and (int(t[0][1::])>28):
                df = pd.read_csv(v, sep=',')
                live = np.array(df['FITC-A'][df['state']=='live'])
                liveEvents = np.append(liveEvents, live)

t_live = np.arcsinh(liveEvents/150)
gm = GaussianMixture(n_components=2, random_state=0).fit(t_live.reshape(-1,1))

pred = {}

for i,v in enumerate(predictions):
    if 'labelled_bhri\\' in v:
        grp = 'bhri'
    elif 'labelled_bhbt\\' in v:
        grp = 'bhbt'
    else:
        grp = 'bhbtri'

    if v.split('.')[-1] == 'csv':
        t = v.split('\\')[-1].split('.')[0].split('_')
        
        if (len(t) == 4) and (t[0][0] == 't') and (t[1][0]=='d'):
            
            
            
            if (t[3]=='bhA') or (t[3]=='bhB') or (t[3]=='bhC'):
                
                df = pd.read_csv(v, sep=',')
                live = np.array(df['FITC-A'][df['state']=='live'])
                live_trans = np.arcsinh(live/150)
                p = gm.predict(live_trans.reshape(-1,1))
                pred[i] = dict(time = int(t[0][1::]), dilution = int(t[1][1::]), condition = t[2], group=t[3], groupLabel =grp, blanks = sum(df['state']=='bkgr'), live=sum(df['state']=='live'), inactive = sum(df['state']=='dead'), groupA = sum(p==0), groupB = sum(p==1))

df = pd.DataFrame.from_dict(pred).T

#df.to_csv(Results)
                
                