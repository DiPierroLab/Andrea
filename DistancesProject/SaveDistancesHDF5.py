#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 13:33:46 2023

@author: andreafalcon
"""


######## LIBRERIAS NECESARIAS PARA CORRER EL PROGRAMA

import pandas as pd
import numpy as np
import datetime
import h5py as hdf

start_time = datetime.datetime.now()

######## FUNCIONES VARIAS PARA EVITAR CODIGO REPETITIVO

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
        
        #d = {'g1g2_name':[], 'OG_g1g2':[],'discrete_distance':[]}
        #d = {'g1g2_name':[], 'OG_g1g2':[], 'gap_distance':[]}
        d = {'g1g2_name':[], 'OG_g1g2':[], 'discrete_distance':[] ,'gap_distance':[]}
    
        for j in range(len(df_list[i])-1):
            for k in range(j+1, len(df_list[i])):
                d['g1g2_name'].append([df_list[i].gene_name[j],df_list[i].gene_name[k]])
                d['OG_g1g2'].append([int(df_list[i].orthogroup_id[j]),int(df_list[i].orthogroup_id[k])])
                d['discrete_distance'].append(df_list[i].index[k]-df_list[i].index[j]-1)
                d['gap_distance'].append(df_list[i].end[k]-df_list[i].start[j])
        
        df_chrom_list.append(pd.DataFrame(data=d))
        
    df_chrom = pd.concat(df_chrom_list, ignore_index=True)    
    
    return df_chrom

######### CARGAR LAS ANOTACIONES DE LA ESPECIE

dir = '/Users/andreafalcon/Documents/GeneAnnotation/'

species = pd.read_csv('/Users/andreafalcon/Dropbox/PostDoc_NEU/DNAZoo/GeneAnnotation_SpeciesOG.csv')

csv = list(species.Specie)

csv.sort()

csv_aux = list(map(lambda x: x.replace('.csv', ''), csv))

df = pd.DataFrame(columns=['g1g2_name', 'OG_g1g2', 'discrete_distance', 'gap_distance'])

for j in range(len(csv)):
    
    df_sp = pd.read_csv(dir + csv[j])
    
    if(len(df_sp) != 0):
    
        df_pairgenes = distance(df_sp[0:100])
    
        df = pd.concat([df, df_pairgenes], axis = 0, ignore_index= True)    
    
    
    
og_list = list(df.OG_g1g2)

tuplas = [tuple(x) for x in  og_list]

unique_og = [list(x) for x in set(tuplas)]

d = {'key':[], 'len_species':[]}
    
for value in unique_og:

    filtro = df['OG_g1g2'].apply(lambda x: sorted(x)).isin([sorted(x) for x in [value]])
    df_filter = df[filtro]
    
    df_filter = df_filter.drop('OG_g1g2', axis = 1)
    
    d['len_species'].append(len(df_filter))
    
    key = 'OG' + str(value[0]) + '-' + 'OG' + str(value[1]) 
    
    d['key'].append(key)
    df_info = pd.DataFrame(data=d)

    df_filter.to_hdf('/Users/andreafalcon/Documents/DistancesAmongGenes/AllDistances.h5', key=key)
    df_info.to_csv('/Users/andreafalcon/Documents/DistancesAmongGenes/AllDistances_Info.csv', encoding='utf-8', index=False, header=True)


        
end_time = datetime.datetime.now()  # Detener el contador de tiempo

elapsed_time = end_time - start_time  # Calcular el tiempo transcurrido

days, seconds = elapsed_time.days, elapsed_time.seconds
hours = seconds // 3600
minutes = (seconds % 3600) // 60
seconds = seconds % 60

print(f"Tiempo de ejecucion: {days} días, {hours} horas, {minutes} minutos, {seconds} segundos")