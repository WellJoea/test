#!/usr/bin/env python
# -*- encoding: utf-8 -*-

import pandas as pd
import numpy  as np

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype']  = 42
import matplotlib.pyplot as plt
import seaborn as sns

class Visal():
    def __init__(self):
        pass
    
    def query_length(self, _indf, out, X='query_length', Dup='query_name', log=False, title=''):
        if not _indf.empty:
            indef = _indf.copy()
            if Dup:
                indef = indef[[Dup, X]].drop_duplicates(keep='first')

            indef[X] = indef[X].astype(int)
            dp = sns.displot(data=indef, x=X, kde=True, log_scale=log)
            dp.set_xticklabels(rotation=270)

            if title:
                plt.title(title)

            plt.tight_layout()
            plt.savefig( out )
            plt.close()

    def clustmap(self, _indf, out, Trans=False):
        linewidths= 0 if min(_indf.shape) > 60  else 0.01
        hm = sns.clustermap(_indf,
                    method='complete',
                    metric='euclidean',
                    z_score=None,
                    #figsize=figsize,
                    linewidths=linewidths,
                    cmap="viridis_r",
                    cbar_pos=(0.02, 0.83, 0.03, 0.11)
                    #center=0,
                    #fmt='.2f',
                    #square=True, 
                    #cbar=True,
                    #yticklabels=Xa,
                    #xticklabels=Xa,
                    #vmin=-1.1,
                    #max=1.1,
                    )
        hm.savefig(out)
        #hm.fig.subplots_adjust(right=.2, top=.3, bottom=.2)
        plt.close()
