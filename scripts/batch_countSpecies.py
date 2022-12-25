# -*- coding: utf-8 -*-
"""
Created on Sun Feb 13 02:44:14 2022

@author: u0139894
"""

from ClassifyCoculture_class import *

import plotly.io as pio
pio.renderers.default='browser'
import plotly.graph_objects as go

#the .csv files from the previous script
dataFolder =  'G:\\My Drive\\FC\\LabelledFiles\\labelled_bhbtri'

#folder for the composition prediction
folderForLabelledData =  'G:\\My Drive\\FC\\PredictedFiles\\bhbtri'

try:
    os.makedirs(folderForLabelledData)

except FileExistsError:
    pass

monocultures = ['riA', 'riB', 'bhA', 'bhB', 'btA', 'btB']

experimentF = os.listdir(dataFolder)

experimentR = os.listdir(folderForLabelledData)

for i,v in enumerate(experimentF):
    
    
    if v.split('.')[-1] == 'csv':
        t = v.split('.')[0].split('_')
        
        if (len(t) == 4) and (t[0][0] == 't') and (t[1][0]=='d'):
        
            print(t)
            
            title = t[0]+ '_' + t[1] + '_' + t[2] + '_' + t[3] + '.csv'
            
            if title in experimentR:
                print('already processed')
            else:
                
                c_monocultures = monocultures[:]
                while t[3] in c_monocultures:
                    c_monocultures.remove(t[3])
                    
                cc = ClassifyCocultures(dataFolder, int(t[0][1::]), t[2], t[3], *c_monocultures, removeCont=True)
                #cc.makePlot()
                cc.writeTable(folderForLabelledData)
            

