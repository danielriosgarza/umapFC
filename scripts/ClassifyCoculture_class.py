# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 15:07:10 2022

@author: u0139894
"""
import plotly.io as pio
pio.renderers.default='browser'
import plotly.graph_objs as go

import os
import numpy as np
import pandas as pd
import fcsparser
import os
import umap

from LabelExperiment_class import LabelExperiment
from sklearn.neighbors import KNeighborsClassifier
from sklearn.semi_supervised import LabelSpreading


class ClassifyCocultures:
    def __init__(self, folder, time, condition, group, *monocultures,
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
                             'FSC-Width'],
                 colors = {'bh': '#083e12', 
                           'bt' : '#f887ff',
                           'ri' : '#ff6e27'
                           }, removeCont = False):
        
        self.folder = folder
        self.time = time
        self.condition = condition
        self.group = group
        self.monocultures = monocultures
        self.channels = channels
        self.colors = colors
        self.contamination = None
        self.species = list(set([i[0:2] for i in self.monocultures]))
        
        self.__getFiles()
        
        
        self.dilution = int(self.experimentFile.split('_')[1][1::])
        self.title = 't' + str(self.time) + '_d' + str(self.dilution) + '_' + self.condition + '_' + self.group
        
        self.__getTables()
        
        
        
        self.nMonoculturePoints = sum([self.speciesCounts[i] for i in self.species])
        self.nExperimentPoints = len(self.dataTable) - self.nMonoculturePoints
        
        if removeCont:
            self.removeContaminants()
        
        self.sData = np.arcsinh(self.dataTable/150)
        self.swData = LabelExperiment.whiten(self.sData)
        self.embbeding = umap.UMAP(n_components=3, n_neighbors=25, min_dist=0.1, metric = 'euclidean', n_jobs=1).fit(self.swData, self.umapLabels).embedding_
        
        self.prediction = self.nnClassifier(self.embbeding[0:self.nMonoculturePoints], self.umapLabels[0:self.nMonoculturePoints], self.embbeding)
    
    def __getFiles(self):
        files = os.listdir(self.folder)
        
        monocultureFiles = {}
        
        
        for i in files:
            if i.split('.')[-1] == 'csv':
                p = i.split('.')[0].split('_')
                if (p[0][0] == 't') and (p[1][0] == 'd'):
                    if int(p[0][1::])==self.time:
                        if p[-1]== self.group:
                            self.experimentFile = i
                        if p[-1] in self.monocultures:
                            monocultureFiles[p[-1]] = i
        self.monocultureFiles = monocultureFiles
        
    
    def __extractData(self, pdDf):
        return np.array(pdDf[self.channels][pdDf['state']=='live'])
        
                    
    
    def __getTables(self):
        self.experimentTable = pd.read_csv(os.path.join(self.folder, self.experimentFile), sep=',')
        
        self.monocultureTables = {i:pd.read_csv(os.path.join(self.folder, self.monocultureFiles[i]), sep=',') for i in self.monocultures}
        
        self.experimentArray = self.__extractData(self.experimentTable)
        self.monocultureArrays = {i:self.__extractData(self.monocultureTables[i][0:min(len(self.monocultureTables[i]), 10000)]) for i in self.monocultureTables}
        
        data = []
        umapLabels = []
        self.speciesCounts = {i:0 for i in self.species}
        self.speciesCounts[self.title] = len(self.experimentArray)
        
        
        
        for i, species in enumerate(self.species):
            for monoC in self.monocultureArrays:
                if monoC[0:2] == species:
                    count = len(self.monocultureArrays[monoC])
                    print(species, count)
                    self.speciesCounts[species] += count
                    t = [i]*count
                    umapLabels +=t
                    data.append(self.monocultureArrays[monoC])
        
        data.append(self.experimentArray)
        
        self.dataTable = np.concatenate(data)
        
        umapLabels += [-1]*self.speciesCounts[self.title]
        self.umapLabels = np.array(umapLabels)
    
    def makeTrainingSpace(self, size=1000, plot=False):
        selected = np.random.choice(np.arange(self.nMonoculturePoints), size= min(size, self.nMonoculturePoints), replace = False)
        tX = self.dataTable.copy()
        tX = np.arcsinh(tX/150)
        tX = LabelExperiment.whiten(tX)
        
        tX = tX[selected]
        
        tY = self.umapLabels[selected]
        space = umap.UMAP(n_components=3, n_neighbors=25, min_dist=0.1, metric = 'euclidean', n_jobs=1).fit(tX, tY)
        
        if plot:
            x = space.embedding_.T[0]
            y = space.embedding_.T[1]
            z = space.embedding_.T[2]
            
            data = []
            
            for i,v in enumerate(self.species):
                    
                data.append(go.Scatter3d(x=x[tY==i], y=y[tY==i], z=z[tY==i],mode='markers', name=v, marker=dict(size=1.0, color=self.colors[v])))
            fig = go.Figure(data=data)
        
            fig.update_layout(title_text='UmapTrainingSpace')
            fig.update_layout(legend= {'itemsizing': 'constant'})
            fig.show()
        return space, tY
        
    
    
    def mapContaminants(self):
        space, tY = self.makeTrainingSpace()
        
        truth = self.umapLabels[0:self.nMonoculturePoints]
        testData = self.dataTable[0:self.nMonoculturePoints]
        testData = np.arcsinh(testData/150)
        testData = LabelExperiment.whiten(testData)
        
        
        emb = space.transform(testData)
        labels = self.nnClassifier(space.embedding_, tY, emb)
        self.contamination  = (labels==truth)
        
        self.contamination = np.concatenate((self.contamination, np.array([True]*self.nExperimentPoints)))
        return self.contamination
        
    def removeContaminants(self):
        contaminants = self.mapContaminants()
        
        self.dataTable = self.dataTable[contaminants]
        
        self.umapLabels = self.umapLabels[contaminants]
        
        self.nMonoculturePoints -=sum(contaminants==False)
    
    
    def writeTable(self, filePath):
        
        t = self.experimentTable.copy()
        
        monocultureD = {i:v for i,v in enumerate(self.species)}
        
        p = np.array([monocultureD[i] for i in self.prediction[self.nMonoculturePoints::]])
        
        t['population'] = ''
        
        t['population'][t['state']=='live'] = p.copy()
        
        t.to_csv(os.path.join(filePath, self.title + '.csv'))
        
        
        
        
        
        
        
        
        
        
        
    def makePlot(self):
        x = self.embbeding.T[0]
        y = self.embbeding.T[1]
        z = self.embbeding.T[2]
        
        xA = x[self.nMonoculturePoints::]
        yA = y[self.nMonoculturePoints::]
        zA = z[self.nMonoculturePoints::]
        
        xB = x[0:self.nMonoculturePoints]
        yB = y[0:self.nMonoculturePoints]
        zB = z[0:self.nMonoculturePoints]
        
        lA = self.prediction[self.nMonoculturePoints::]
        data = []
        
        for i,v in enumerate(self.species):
                
            data.append(go.Scatter3d(x=x[self.umapLabels==i], y=y[self.umapLabels==i], z=z[self.umapLabels==i],mode='markers', name=v, marker=dict(size=1.0, color=self.colors[v])))
            
            
            
            
            data.append(go.Scatter3d(x=xA[lA==i], y=yA[lA==i], z=zA[lA==i],mode='markers', name='predicted '+ v, marker=dict(size=2.0, color=self.colors[v], opacity=.5)))
        
        
        
        fig = go.Figure(data=data)
    
        fig.update_layout(title_text='UmapGating' + ' ' + self.title)
        fig.update_layout(legend= {'borderwidth':1, 'font' : dict(size= 8), 'itemwidth' : 30, 'itemsizing' : 'constant', 'xanchor' : 'left', 'orientation': 'h'})
        fig.show()
    
    @staticmethod
    def nnClassifier(train_x, train_y, predict_x):
        #nbrs = KNeighborsClassifier(n_neighbors=50, algorithm='auto', weights='distance', metric='mahalanobis', n_jobs=8).fit(train_x, train_y)
        nbrs = LabelSpreading(kernel='knn', tol=1e-7).fit(train_x, train_y)
        
        return nbrs.predict(predict_x)
        
        
        