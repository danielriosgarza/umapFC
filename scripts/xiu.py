# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 11:55:19 2023

@author: danie
"""

from LabelExperiment_class import *



filePath = 'C:\\Users\\danie\\OneDrive\\Documentos\\xiu\\'


channels= ['FSC',
           'SSC',
           'FL1',
           'FL2',
           'FL3',
            'FL4',
            'FL5',
            'FL6',]


blankP = 100
cellP = 100

le = LabelExperiment(filePath, 
                     group='riverA', 
                     time =0,
                     condition = 'river',
                     blanks = ['20230621125452245-t0_d10_river_blank.fcs'],
                     sgT=None,
                     piT=None,
                     channels= channels,
                     sybr_channel='FL1',
                     cover = 1.5,
                     fracBlanks = 1/blankP,
                     fracCells = 1/cellP,
                     n_jobs = 8
                     
                     ) 

le.makePlot()

labels = le.state[le.nBlankFilePoints::]
cellN = sum(labels=='cell')*cellP