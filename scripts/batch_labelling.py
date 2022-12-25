# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 12:54:58 2022

@author: u0139894
"""
import gc
from LabelExperiment_class import *

import plotly.io as pio
pio.renderers.default='browser'
import plotly.graph_objects as go


dataFolder = 'G:\\My Drive\\FC\\rawFCfiles\\BHBTRI_batch_wc'

folderForLabelledData = 'G:\\My Drive\\FC\\LabelledFiles\\labelled_bhbtri'

try:
    os.makedirs(folderForLabelledData)

except FileExistsError:
    pass


fileList = os.listdir(dataFolder)
included = []
excluded = []

for i in fileList:
    
    if i.split('.')[-1] =='fcs':
        
        p = i.split('-')[1].split('_')
        if (len(p)==4) and (p[0][0] == 't') and (p[1][0] == 'd'):
            
            if ('blank' not in p[-1]):
                included.append((int(p[0][1::]), p[2], p[3])) #time, condition, group
        else:
            excluded.append(i)

print('exluded files because of name not matching the expected format:','\n', excluded, '\n', 'a total of ', len(included), ' experiments will be processed\n')




cover = 0.5


for i, experiment in enumerate(included):
    

    print('processing: ', experiment, cover)
        
    labeledExp = LabelExperiment(folder = dataFolder, group = experiment[2], time = experiment[0], condition = experiment[1], cover=cover, sgT = 10**4.5, piT = 10**3.1, 
                             channels = ['FSC-A',
                                         'SSC-H',
                                         'SSC-A',
                                         'FITC-H',
                                         'FITC-A',
                                         'PerCP-H',
                                         'PerCP-A',
                                         'APC-H',
                                         'APC-A',
                                         'APC-A750-H',
                                         'APC-A750-A',
                                         'VSSC-H',
                                         'VSSC-A',
                                         'KO525-H',
                                         'KO525-A',
                                         'PE-H',
                                         'PE-A',
                                         'ECD-H',
                                         'ECD-A',
                                         'PI-H',
                                         'PI-A',
                                         'FSC-Width'])
        
        
                
    #labeledExp.makePlot()
        
    labeledExp.writeTable(folderForLabelledData)
        
    gc.collect()