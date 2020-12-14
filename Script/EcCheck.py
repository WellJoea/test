#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd
import numpy  as np
from joblib import Parallel, delayed

from .EcVisual import Visal
from .EcUtilities import Utilities

class CheckBP():
    def __init__(self, arg, log):
        self.arg  = arg
        self.log  = log
        self.mapq = self.arg.minmapQ
        self.overfremin = self.arg.overfremin
        self.overlenmin = self.arg.overlenmin
        self.maxcheck = self.arg.maxchecksofttwoends

    def _getinfo(self, _info):
        self.info = _info
        self.inid = _info.sampleid
        self.inbam = '%s/%s.sorted.bam'%(self.arg.Bam, self.inid)
        self.outdir= '%s/%s'%(self.arg.Cheak, self.inid )
        self.arg.outpre= '%s/%s'%(self.outdir, self.inid)
        os.makedirs(self.outdir, exist_ok=True)
        return self

    def _getbeddb(self):
        self.inbed = '{0}/{1}/{1}.chimeric.bed'.format(self.arg.Fetch, self.inid)
        self.inbed = pd.read_csv( self.inbed, sep='\t', low_memory=False)

        self.inbed[['start', 'end']]  = self.inbed[['start', 'end']].astype(int)
        self.inbed['#chrom']   = self.inbed['#chrom'].astype(str)
        self.inbed['cigarreg'] = self.inbed.cigarreg.map(eval)
        self.inbed['fflag']    = 'DROP'
        self.inbed['raw_order'] = 1

        COLs  = ['#chrom', 'start', 'end',  'SID', 'length', 'forword', 'query_name', 'query_length',
                 'fflag', 'raw_order','query_counts', 'cigarstring',  'cigarreg']
        self.inbed = self.inbed[COLs]
        self.inbed.sort_values(by=['query_name', 'cigarreg', '#chrom', 'start', 'end' ], ascending=[True]*5, inplace=True)
        self.inbed['raw_order'] =  self.inbed.groupby(by=['SID', 'query_name'], sort=False)['raw_order'].apply(np.cumsum)

        self.inBP = pd.read_csv(self.arg.checkbed, sep='\t')
        self.inBP['#chrom'] = self.inBP['#chrom'].astype(str)
        self.inBP['start']  = self.inBP['start'].astype(int) -1
        self.inBP['end']    = self.inBP['end'].astype(int) -1
        self.inBP['lenght'] = self.inBP['end'] - self.inBP['start'] + 1

    def _getkeep(self):
        self.outdir= '%s/%s'%(self.arg.Cheak, 'BPState' )
        self.outpre= '%s/%s'%(self.outdir, 'All.plasmid')
        os.makedirs(self.outdir, exist_ok=True)
        return self

    def Rmerge( self, intervals):
        """
        :param intervals: List[List[int]]
        :return: List[List[int]]
        """
        intervals.sort(key=lambda x: x[0])
        merged = []
        for interval in intervals:
            if not merged or merged[-1][-1] < interval[0]:
                merged.append(interval)
            else:
                merged[-1][-1] = max(merged[-1][-1], interval[-1])
        merged = sum([i[1]-i[0] + 1 for i in merged ])
        return merged

    def BEDfilter1(self, _inbed):
        F = open(_inbed, 'r').readlines()
        F = [ i.strip().split('\t') for i in F if len(i.strip().split('\t')) < 12 ]

        Names = ['#chrom', 'start', 'end', '#chrom_r', 'start_r', 'end_r', 'query_name', 'mapq', 'forword']
        inbed = pd.DataFrame(F, columns=Names)
        inbed[['#chrom', '#chrom_r']] = inbed[['#chrom', '#chrom_r']].astype(str)
        inbed[['start', 'end', 'start_r', 'end_r', 'mapq']] = inbed[['start', 'end', 'start_r', 'end_r', 'mapq']].astype(int)

        checkb = pd.read_csv(self.arg.checkbed, sep='\t')
        checkb['#chrom'] = checkb['#chrom'].astype(str)
        checkb[['start', 'end']] = checkb[['start', 'end']].astype(int)

        inbed = inbed.merge(checkb, on=['#chrom', 'start', 'end'], how='left')

        #####addinfor
        inbed['SID']      = self.inid
        inbed['end_o']    = inbed[['end', 'end_r']].min(1)
        inbed['start_o']  = inbed[['start', 'start_r']].max(1)
        inbed['BPlength'] = inbed['end'] - inbed['start'] + 1
        Names = ['query_name', '#chrom', 'start', 'end', 'SID', 'Links', 'BPlength', 'start_o', 'end_o', 
                 '#chrom_r', 'start_r', 'end_r', 'mapq', 'forword']
        inbed = inbed[Names]
        inbed.sort_values(by=['SID', 'query_name', '#chrom', 'start'], inplace=True)

        #####filter
        inbed1 = inbed[(inbed.mapq >=self.mapq)]
        inbed1 = inbed1.merge(inbed1.groupby(by=['Links', 'query_name', 'start'])\
                                 .apply(lambda x: self.Rmerge(x[['start_o', 'end_o']].values.tolist()))\
                            .to_frame('OVERlen').reset_index(), on=['Links', 'query_name', 'start'], how='left')
        inbed1['OVERfre'] = (inbed1['OVERlen']/inbed1['BPlength']).round(4)
        inbed1 = inbed1[((inbed1.OVERfre >=self.overfremin) | (inbed1.OVERlen >=self.overlenmin))]

        inbed1 = inbed1.merge(inbed1.groupby(by=['Links', 'query_name']).size()\
                                    .to_frame('query_count').reset_index(), on=['Links', 'query_name'], how='left')
        inbed1 = inbed1.merge(inbed1.groupby(by=['Links', 'query_name'])['start'].unique().apply(lambda x:len(x))\
                                    .to_frame('BP_count').reset_index(), on=['Links', 'query_name'], how='left')
        inbed1 = inbed1[((inbed1.query_count >=2) & (inbed1.BP_count >=2))]

        #inbed1 = inbed1.merge(inbed1.groupby(by=['Links', 'query_name'])['forword'].unique().apply(lambda x:len(x))\
        #                            .to_frame('Forwrd_count').reset_index(), on=['Links', 'query_name'], how='left')
        #inbed1 = inbed1[(inbed1.Forwrd_count<=1)]

        inbed2 = pd.concat([inbed, inbed1[Names]], axis=0, sort=False).drop_duplicates(keep=False)

        return inbed1, inbed2

    def BPFetchBed1(self, _inline):
        self._getinfo(_inline)
        self.log.CI('start schecking breakpoin region: ' + self.inid)

        bambed = self.arg.outpre + '.breakpoit.reads.txt'
        if not os.path.exists(bambed):
            Utilities(self.arg, self.log).bambedflter(self.inbam, bambed)
        bedKeep, bedDrop = self.BEDfilter(bambed)
        bedKeep.to_csv(self.arg.outpre + '.breakpoint.Keep.txt', sep='\t', index=False)
        bedDrop.to_csv(self.arg.outpre + '.breakpoint.Drop.txt', sep='\t', index=False)
        self.log.CI('finish schecking breakpoin region: ' + self.inid)

    def BPStat1(self, _info ):
        Parallel( n_jobs=self.arg.njob, verbose=1 )( delayed( self.BPFetchBed1 )(_l) for _n, _l in _info.iterrows())
        self.log.CI('start stating all samples region.')
        self._getkeep()
        BPKEEP = []
        for _n, _l in _info.iterrows():
            EMerge = '{0}/{1}/{1}.breakpoint.Keep.txt'.format(self.arg.Cheak, _l.sampleid)
            if os.path.exists( EMerge ):
                BPKEEP.append( pd.read_csv(EMerge, sep='\t', header=0) )
            else:
                self.log.CW('cannot find the file: '+ EMerge)
        if BPKEEP:
            BPKEEP = pd.concat(BPKEEP, axis=0,sort=False)
            BPKEEP['#chrom'] = BPKEEP['#chrom'].astype(str)
            BPKEEP.to_csv(self.outpre+'.Keep', sep='\t', index=False)
            self.BPKEEP(BPKEEP)
        else:
            self.log.CW('cannot find the valid files.')
        self.log.CI('finish stating all samples region.')

    def BPKEEP1(self, _indf, Lplot=True):
        indf = _indf.groupby(by=['#chrom', 'Links', 'SID', ])['query_name']\
                    .unique().apply(lambda x:len(x)).to_frame('support_ID_num').reset_index()

        pvot = indf.pivot(index='#chrom', columns='SID', values='support_ID_num').fillna(0).astype(int)
        indf = indf.groupby(by=['#chrom', 'Links'])['support_ID_num'].sum().to_frame('support_num').reset_index()
        indf = indf.merge(pvot.reset_index(), on='#chrom').sort_values(by=['support_num', '#chrom'], ascending=[False, True])
        indf.to_csv(self.outpre+'.Keep.matrix', sep='\t', index=False)

        if Lplot:
            Visal().clustmap(pvot, self.outpre+'.Keep.matrix.pdf')
            Visal().clustmap(np.log2(pvot+1), self.outpre+'.Keep.matrix.log2.pdf')

    def BEDfilter(self, inbed):
        #####addinfor
        inbed['cigarreg'] = inbed.cigarreg.map(eval)
        inbed['start_o']  = inbed[['start', 'start_i']].max(1)
        inbed['end_o']    = inbed[['end', 'end_i']].min(1)
        inbed.sort_values(by=['SID', 'query_name', 'raw_order'], inplace=True)

        inbed.rename(columns={'Links_i':'Links'}, inplace=True)
        #####addstander
        GRPBy = ['Links', 'query_name', 'forword', 'end_i']
        inbed = inbed.merge(inbed.groupby(by=GRPBy)\
                                 .apply(lambda x: self.Rmerge(x[['start_o', 'end_o']].values.tolist()))\
                            .to_frame('OVERlen').reset_index(), on=GRPBy, how='left')
        inbed['OVERfre'] = (inbed['OVERlen']/inbed['lenght_i']).round(4)

        GRPBy = inbed.groupby(by=['Links', 'query_name'])
        GROUP = [ GRPBy['end_i'].unique().apply(lambda x:len(x)).to_frame('BP_count'),
                  GRPBy['cigarreg'].first().str[0].to_frame('HeadSoft'),
                  GRPBy['cigarreg'].last().str[1].to_frame('TailSoft') ]
        GROUP = pd.concat(GROUP, axis=1, sort=False).reset_index()
        inbed = inbed.merge(GROUP, on=['Links', 'query_name'], how='left')

        inbed['HeadSoft'] = (inbed['HeadSoft']/inbed['query_length']).round(4)
        inbed['TailSoft'] = (1 - inbed['TailSoft']/inbed['query_length']).round(4)
        
        # add marker
        inbed.loc[((inbed.OVERfre < self.overfremin) & (inbed.OVERlen < self.overlenmin)), 'fflag'] += ';OVERMIN'
        inbed.loc[(inbed.BP_count  < 2), 'fflag'] += ';BPLOW'
        inbed.loc[((inbed.HeadSoft > self.maxcheck) | (inbed.TailSoft > self.maxcheck)),   'fflag'] += ';HEADTAIL'
        inbed.loc[(inbed.fflag=='DROP'), 'fflag'] = 'KEEP'
        return inbed

    def BPFetchBed(self, _inline):
        self._getinfo(_inline)
        self._getbeddb()
        intSect = Utilities(self.arg, self.log)\
                    .bedintersect(self.inbed, self.inBP, s=False, S=False, wa=True, wb=True)
        intSect.to_csv(self.arg.outpre + '.breakpoint.bed.txt', sep='\t', index=False)
        intSect = self.BEDfilter(intSect)
        intSect.to_csv(self.arg.outpre + '.breakpoint.Mark.txt', sep='\t', index=False)
        intSect = intSect[(intSect.fflag=='KEEP')]
        intSect.to_csv(self.arg.outpre + '.breakpoint.Keep.txt', sep='\t', index=False)

    def BPKEEP(self, _indf, Lplot=True):
        indf = _indf.groupby(by=['#chrom', 'Links', 'SID', ])['query_name']\
                    .unique().apply(lambda x:len(x)).to_frame('support_ID_num').reset_index()

        pvot = indf.pivot(index='#chrom', columns='SID', values='support_ID_num').fillna(0).astype(int)
        indf = indf.groupby(by=['#chrom', 'Links'])['support_ID_num'].sum().to_frame('support_num').reset_index()
        indf = indf.merge(pvot.reset_index(), on='#chrom').sort_values(by=['support_num', '#chrom'], ascending=[False, True])
        indf.to_csv(self.outpre+'.Keep.matrix', sep='\t', index=False)

        if Lplot:
            Visal().clustmap(pvot, self.outpre+'.Keep.matrix.pdf')
            Visal().clustmap(np.log2(pvot+1), self.outpre+'.Keep.matrix.log2.pdf')

    def BPStat(self, _info ):
        Parallel( n_jobs=self.arg.njob, verbose=1 )( delayed( self.BPFetchBed )(_l) for _n, _l in _info.iterrows())

        self.log.CI('start stating all samples region.')
        self._getkeep()
        BPKEEP = []
        for _n, _l in _info.iterrows():
            EMerge = '{0}/{1}/{1}.breakpoint.Keep.txt'.format(self.arg.Cheak, _l.sampleid)
            if os.path.exists( EMerge ):
                BPKEEP.append( pd.read_csv(EMerge, sep='\t', header=0) )
            else:
                self.log.CW('cannot find the file: '+ EMerge)
        if BPKEEP:
            BPKEEP = pd.concat(BPKEEP, axis=0,sort=False)
            BPKEEP['#chrom'] = BPKEEP['#chrom'].astype(str)
            BPKEEP.to_csv(self.outpre+'.Keep', sep='\t', index=False)
            self.BPKEEP(BPKEEP)
        else:
            self.log.CW('cannot find the valid files.')
        self.log.CI('finish stating all samples region.')

    def PlotLM(self, _info):
        self._getkeep()
        SampID = _info.sampleid.tolist()
        ecCOL = ['plasmid', 'support_num']  + SampID

        indf = pd.read_csv(self.outpre+'.Keep.matrix', sep='\t')
        indf.rename(columns={'#chrom' : 'plasmid'}, inplace=True)
        indf['plasmid'] = indf['plasmid'].str.upper()

        dPCR = '/data/zhouwei/01Projects/03ecDNA/Nanopore/spikei.info.txt'
        dPCR = pd.read_csv(dPCR, sep='\t')
        dPCR['plasmid'] = dPCR['plasmid'].str.upper()

        ecdf = self.arg.outdir + '/04.EcRegion/All.circle.region.UpMerge_sort'
        ecdf = pd.read_csv(ecdf, sep='\t')
        ecdf.rename(columns={'#chrom' : 'plasmid'}, inplace=True)
        ecdf['plasmid'] = ecdf['plasmid'].str.upper()


        indf = indf.merge(dPCR, on='plasmid', how='right')
        indf.to_csv('./aa.xls', sep='\t', index=False)
        print(indf)
