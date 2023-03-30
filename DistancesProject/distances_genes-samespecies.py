#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 12:47:04 2023

@author: andreafalcon
"""

import pandas as pd
import matplotlib.pyplot as plot
import seaborn as sb

def chrom_sep (df):
    
    hi_c = df.seq_id.unique()
    
    df_list = []
    
    for h in range(len(hi_c)):
        
        df_list.append(df.loc[df['seq_id'] == hi_c[h]])
        df_list[h] = df_list[h].reset_index(drop = True)
        
    return df_list

def distance (df):
    
    df_list = chrom_sep(df)
            
    df_chrom_list = []
    
    for i in range(len(df_list)):
        
        d = {'chrom':[],'g1_id':[], 'g2_id':[], 'OG_g1':[], 'OG_g2':[] ,'discrete_distance':[], 'gap_distance':[]}
        #d = {'chrom':[],'g1_id':[], 'g2_id':[], 'OG_g1g2':[],'discrete_distance':[], 'gap_distance':[]}
    
        for j in range(len(df_list[i])-1):
            for k in range(j+1, len(df_list[i])):
                d['chrom'].append(i+1)
                d['g1_id'].append(df_list[i].gene_id[j])
                d['g2_id'].append(df_list[i].gene_id[k])
                #d['OG_g1g2'].append(sorted([df_list[i].orthogroup_id[j],df_list[i].orthogroup_id[k]]))
                d['OG_g1'].append(df_list[i].orthogroup_id[j])
                d['OG_g2'].append(df_list[i].orthogroup_id[k])
                d['discrete_distance'].append(df_list[i].index[k]-df_list[i].index[j]-1)
                d['gap_distance'].append(df_list[i].end[k]-df_list[i].start[j])
        
        df_chrom_list.append(pd.DataFrame(data=d))
        
    df_chrom = pd.concat(df_chrom_list, ignore_index=True)    
    
    return df_chrom


def graph (df1, df2,variable):
    
    # maximo = max(max(df1[variable]), max(df2[variable]))
    # minimo = min(min(df1[variable]), min(df2[variable]))
    
    # intervalos = range(minimo, maximo + 2)
    
    sb.set_style("darkgrid")

    sb.histplot(df1[variable], color='dodgerblue', kde=True)
    sb.histplot(df2[variable], color='lightgreen', kde=True)

    #configuramos en Matplotlib
    plot.ylabel('Frequency')
    plot.xlabel(variable)
    plot.title('')
    
    plot.legend(labels=['Cheetah', 'Panthera'])

    plot.show()
    
    return

def boxgraph (df1, df2,variable):
    
    data = [df1[variable], df2[variable]]
    
    labels = ['Cheetah','Panthera']
    
    fig = plot.figure(figsize =(10, 7))
    ax = fig.add_subplot(111)
    
    bplot = ax.boxplot(data,
                       vert=True,  # vertical box alignment
                       patch_artist=True,  # fill with color
                       labels=labels)  # will be used to label x-ticks
    
    colors = ['dodgerblue', 'lightgreen']
 
    for patch, color in zip(bplot['boxes'], colors):
        patch.set_facecolor(color)
        
    for flier in bplot['fliers']:
        flier.set(marker ='D',
                  color ='grey',
                  alpha = 0.2)
        
    ax.yaxis.grid(True)
    ax.set_ylabel(variable)
    
    plot.show()
    
    return


df_sp1 = pd.read_csv('/Users/andreafalcon/Downloads/aciJub1.csv')
df_sp2 = pd.read_csv('/Users/andreafalcon/Downloads/PanPar1.csv')

ortho_intersection = list(set(df_sp1.orthogroup_id).intersection(set(df_sp2.orthogroup_id)))

df_s1_ortho = df_sp1[df_sp1.orthogroup_id.isin(ortho_intersection)]
df_s2_ortho = df_sp2[df_sp2.orthogroup_id.isin(ortho_intersection)]


df_s1_ortho = df_s1_ortho.reset_index(drop=False)
df_s2_ortho = df_s2_ortho.reset_index(drop=False)

df_dis_s1 = distance(df_s1_ortho)
df_dis_s2 = distance(df_s2_ortho)

del df_sp1
del df_sp2

del df_s1_ortho
del df_s2_ortho

names = df_dis_s1.columns.values

df_ortho = pd.merge(df_dis_s1, df_dis_s2, on=['OG_g1','OG_g2'], how  = 'inner' )

names_ortho = df_ortho.columns.values

del df_dis_s1
del df_dis_s2

df_ortho_s1 = df_ortho[names_ortho[0:7]]
df_ortho_s2 = df_ortho[names_ortho[[7,8,9,3,4,10,11]]]

df_ortho_s1.columns = names
df_ortho_s2.columns = names

del df_ortho

df_ortho_s1 = df_ortho_s1.drop_duplicates()
df_ortho_s2 = df_ortho_s2.drop_duplicates()

sb.jointplot(df_ortho_s1.discrete_distance, df_ortho_s2.discrete_distance, cmap='hot_r', kind = 'kde', fill = True)
sb.jointplot(df_ortho_s1.gap_distance, df_ortho_s2.gap_distance, cmap='hot_r', kind = 'kde', fill = True)


graph(df_ortho_s1, df_ortho_s2, 'discrete_distance')
graph(df_ortho_s1, df_ortho_s2, 'gap_distance')

boxgraph(df_ortho_s1, df_ortho_s2, 'discrete_distance')
boxgraph(df_ortho_s1, df_ortho_s2, 'gap_distance')

sb.jointplot(df_ortho_s1.discrete_distance, df_ortho_s2.discrete_distance, cmap='hot_r', kind = 'kde', fill = True)
sb.jointplot(df_ortho_s1.gap_distance, df_ortho_s2.gap_distance, cmap='hot_r', kind = 'kde', fill = True)
