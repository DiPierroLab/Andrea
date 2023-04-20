#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 13:50:43 2023

@author: andreafalcon
"""

import pandas as pd
import os

def clean_df(df):
    
    for i in df.columns.values:
        
        if((i == 'ChromA') | (i == 'ChromB') | (i == 'OGA') | (i == 'OGB')):
            
            df[i] = df[i].str.replace('A', '')
            df[i] = df[i].str.replace('B', '')
            
    return df

dir = '/Users/andreafalcon/Documents/all_dagchainer_results/'

contenido = os.listdir(dir)

tsv = []

for fichero in contenido:
    if os.path.isfile(os.path.join(dir, fichero)) and fichero.endswith('.aligncoords'):
        tsv.append(fichero)


empty = []

for j in tsv:

    data = []
    
    with open(dir+j, 'r') as f:
        for line in f:
            if not line.startswith('##'):
                fields = line.strip().split('\t')
                data.append(fields)
                
    df = pd.DataFrame(data, columns=['ChromA', 'OGA', 'GeneIndexA', 'GeneIndexAA', 'ChromB', 'OGB', 'GeneIndexB', 'GeneIndexBB', 'Score0', 'Score'])
    
    if(len(df) != 0):
    
        df = df.drop(['GeneIndexAA', 'GeneIndexBB' ,'Score0', 'Score'], axis = 1)
        
        df = df.drop_duplicates()
                
        df_split = df.loc[df['ChromA'].str.contains('B')]        
        
        diff_df = df.merge(df_split, how='outer', indicator=True)
        diff_df = diff_df.loc[diff_df['_merge'] == 'left_only'].drop('_merge', axis=1)
        
        df_split.columns = ['ChromB', 'OGB', 'GeneIndexB' ,'ChromA', 'OGA', 'GeneIndexA']
        df_split = df_split[['ChromA', 'OGA', 'GeneIndexA', 'ChromB', 'OGB', 'GeneIndexB']]  
        
        df_unsplit = pd.concat([diff_df, df_split], axis=0)
        df_unsplit = df_unsplit.reset_index(drop=True)
        
        df_unsplit = clean_df(df_unsplit)
        
        df_unsplit['ChromA'] = df_unsplit.apply(lambda row: str(j.split('-')[0] + '_Ch-' + row['ChromA']), axis=1)
        df_unsplit['ChromB'] = df_unsplit.apply(lambda row: str(j.split('-')[1] + '_Ch-' + row['ChromB']), axis=1)
             
        df_unsplit['GeneA'] = df_unsplit.apply(lambda row: str(j.split('-')[0] + '_G-' + row['GeneIndexA']), axis=1)
        df_unsplit['GeneB'] = df_unsplit.apply(lambda row: str(j.split('-')[1] + '_G-' + row['GeneIndexB']), axis=1)
        
        df_net = df_unsplit[['ChromA','GeneA', 'OGA' ,'ChromB','GeneB', 'OGB']]
        
        name_file = j.split('-')[0] + '-' + j.split('-')[1] + '.csv' 
        
        df_net.to_csv('/Users/andreafalcon/Documents/NetworkFiles-dragchainer/'+name_file, index=False, header=True)
        
    else:
        
        empty.append(j)
        
print(len(empty))
