# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 10:19:06 2022

@author: DGarza
"""
import plotly.io as pio
pio.renderers.default='browser'
import plotly.graph_objs as go

import numpy as np
import pandas as pd
import fcsparser
import os
import umap

class LabelExperiment:
    def __init__(self, 
                 folder : str, 
                 group : str, 
                 time : int, 
                 condition : str ='wc', 
                 cover : float = 1.5, 
                 sgT : float = 10**4.5, 
                 piT : float = 10**3.1,
                 channels : list = ['FSC-A',
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
                 n_jobs : int = 8, 
                 blanks : list = None,
                 sybr_channel : str = 'FITC-A',
                 pi_channel : str = 'PI-H',
                 fracCells = 1.0,
                 fracBlanks = 1.0):
        self.folder = folder
        
        self.group = group
        
        self.time = time
        self.condition = condition
        self.cover = cover
        self.sgT = sgT
        self.piT = piT
        self.channels = channels
        self.n_jobs = n_jobs
        
        if self.sgT is not None:
            self.sybrColumn = self.channels.index(sybr_channel)
            self.piColumn = self.channels.index(pi_channel)
        
        self.blanks = blanks
        
        
        
        self.__getDataTables()
        
        self.blankTable = self.blankTable[0: int(len(self.blankTable)*fracBlanks)]
        
        self.experimentTable = self.experimentTable[0: int(len(self.experimentTable)*fracCells)]
        
        
        self.dilution = int(self.experimentFile.split('-')[1].split('_')[1][1::])
        self.title = 't' + str(self.time) + '_d' + str(self.dilution) + '_' + self.condition + '_' + self.group
        
        self.dataTable = np.concatenate((self.blankTable, self.experimentTable))
        self.nExperimentPoints = len(self.experimentTable)
        self.nBlankFilePoints = len(self.blankTable)
        self.nPoints = len(self.dataTable)
        
        self.sData = np.arcsinh(self.dataTable/150)
        self.swData = self.whiten(self.sData)
        
        self.embbeding = umap.UMAP(n_components=3, n_neighbors=500, min_dist=0.01, metric = 'canberra', n_jobs=self.n_jobs).fit(self.swData).embedding_
        
        self.__getLabels()
        self.__getState()
        
    
    
    
    def __getFCtable(self, filePath):
        
        meta, data = fcsparser.parse(filePath, reformat_meta=True)
    
        chData = data[self.channels]
    
        return chData
    
    def __getDataTables(self):
        allFiles = os.listdir(self.folder)
        
        
        
        if self.blanks is None:
            blanks = []
            
        else:
            blanks = self.blanks[:]
        
        
        for i in allFiles:
            
            if i.split('.')[-1] == 'fcs':
                a = i.split('.')[0].split('-')[1].split('_')
                
                print(a)

                if len(a) == 4: 
                    if (int(a[0][1::]) == self.time) and (a[2] == self.condition):
                        if a[3] == 'blank':
                            if self.blanks is None:
                                blanks.append(i)
                        elif a[3] == self.group:
                            
                            self.experimentFile = i
    
        blankTables = [self.__getFCtable(os.path.join(self.folder, i)) for i in blanks]
        self.experimentTable = self.__getFCtable(os.path.join(self.folder,self.experimentFile))
        
        
        self.blankTable = pd.concat(blankTables)
        
    def __getLabels(self):
        
        labels = np.array(['cell']*self.nPoints)

        for i in range(self.nBlankFilePoints):
            cond = np.sqrt(np.sum((self.embbeding-self.embbeding[i])**2, axis=1)) < self.cover * (self.countPointsWithinRadius(self.embbeding[0:self.nBlankFilePoints], self.embbeding[i], self.cover)/self.nBlankFilePoints)
    
            labels[cond] = 'blank'
        
        self.labels = labels
    
    def __getState(self):
        state = np.array(['live']*self.nPoints)
        for i,v in enumerate(self.dataTable):
            
            if self.labels[i] == 'cell':
                state[i] = 'cell'
                
                
                if self.sgT is not None:
                    
                    
                    
                    dead = (v[self.piColumn] > self.piT)
                    
                    debr = (v[self.sybrColumn] <= self.sgT) & (v[self.piColumn] <= self.piT)
                    if dead:
                        state[i] = 'dead'
                    elif debr:
                        state[i] = 'debr'
                
            else:
                state[i] = 'bkgr'
        
        self.state = state
    
    

    def makePlot(self):
        x = self.embbeding.T[0]
        y = self.embbeding.T[1]
        z = self.embbeding.T[2]
        
        data = []
        
        
        if self.sgT is not None:

        

            data.append(go.Scatter3d(x=x[self.state=='live'], y=y[self.state=='live'], z=z[self.state=='live'],mode='markers', name='live cells', marker=dict(size=1.0, color='#39FF14')))
            
            data.append(go.Scatter3d(x=x[self.state=='dead'], y=y[self.state=='dead'], z=z[self.state=='dead'],mode='markers', name='inactive cells<br>(pi positive)', marker=dict(size=1.0, color='#FF10F0')))
            
            data.append(go.Scatter3d(x=x[self.state=='bkgr'], y=y[self.state=='bkgr'], z=z[self.state=='bkgr'],mode='markers', name='background events', marker=dict(size=1.0, color='#1c1c1c')))
            
            data.append(go.Scatter3d(x=x[0:self.nBlankFilePoints], y=y[0:self.nBlankFilePoints], z=z[0:self.nBlankFilePoints],mode='markers', name='blank file events', marker=dict(size=1.0, color='#aaaaaa')))
        
        else:
            
            data.append(go.Scatter3d(x=x[self.state=='cell'], y=y[self.state=='cell'], z=z[self.state=='cell'],mode='markers', name='cells', marker=dict(size=1.0, color='#39FF14')))
            
            
            data.append(go.Scatter3d(x=x[self.state=='bkgr'], y=y[self.state=='bkgr'], z=z[self.state=='bkgr'],mode='markers', name='background events', marker=dict(size=1.0, color='#1c1c1c')))
            
            data.append(go.Scatter3d(x=x[0:self.nBlankFilePoints], y=y[0:self.nBlankFilePoints], z=z[0:self.nBlankFilePoints],mode='markers', name='blank file events', marker=dict(size=1.0, color='#aaaaaa')))
        
        
        
        
        
        fig = go.Figure(data=data)
    
        fig.update_layout(title_text='UmapGating' + ' ' + self.title)
        fig.update_layout(legend= {'borderwidth':1, 'font' : dict(size= 8), 'itemwidth' : 30, 'itemsizing' : 'constant', 'xanchor' : 'left', 'orientation': 'h'})
        fig.update_layout(legend_title_text='<b>populations</b>')
        
        fig.update_layout(scene = dict(
                    xaxis_title='UMAP 1',
                    yaxis_title='UMAP 2',
                    zaxis_title='UMAP 3'),
                    boxgap=0)
        fig.show()
    
    def writeTable(self, filePath):
        data = self.experimentTable.copy()
        data['labels'] = self.labels[self.nBlankFilePoints::]
        data['state'] = self.state[self.nBlankFilePoints::]
        
        data.to_csv(os.path.join(filePath, self.title + '.csv'))
    @staticmethod
    def whiten(X):
    
        X = X.reshape((-1, np.prod(X.shape[1:])))
        X_centered = X - np.mean(X, axis=0)
        Sigma = np.dot(X_centered.T, X_centered) / X_centered.shape[0]

        V_sqrt = np.diag(np.std(X, axis=0))
        P = np.dot(np.dot(np.linalg.inv(V_sqrt), Sigma), np.linalg.inv(V_sqrt))
        G, Theta, _ = np.linalg.svd(P)
        W = np.dot(np.dot(G, np.dot(np.diag(1.0 / np.sqrt(Theta + 1e-5)), G.T)), np.linalg.inv(V_sqrt))

        return np.dot(X_centered, W.T)
    
    @staticmethod
    def countPointsWithinRadius(X, point, radius):
        dist = np.sqrt(np.sum((X-point)**2, axis=1))

        return sum(dist<=radius)
    
        
    
    