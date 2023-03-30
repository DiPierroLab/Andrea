#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 15:52:54 2023

@author: andreafalcon
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 14:47:45 2023

@author: andreafalcon
"""

######## LIBRERIAS NECESARIAS PARA CORRER EL PROGRAMA

import pandas as pd
import statistics as stat
import numpy as np
import datetime
import os

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

######### CARGAR LAS ANOTACIONES DE LA ESPECIE DE REFERENCIA
######### ESPECIE A LA QUE SE LE VAN A RASTREAR LOS PARES DE GENES EN LAS ESPECIES RESTANTES

dir = '/home/j.falcn/GeneAnnotation/'

species = pd.read_csv('/home/j.falcn/Info/GeneAnnotation_SpeciesOG.csv')

csv = list(species.Specie)

csv.sort()

for j in [0]:
    
    csv_aux = list(map(lambda x: x.replace('.csv', ''), csv))

    df_sp = pd.read_csv(dir + csv[j])
    
    df_pairgenes = distance(df_sp[0:1750])
       
    num_col = ['genes_name'] + [csv_aux.pop(j)] + ['orthopair'] +  csv_aux
    
    df_pairtrack_dd = pd.DataFrame(columns = num_col)
    df_pairtrack_gd = pd.DataFrame(columns = num_col)
    
    vector_dd = np.empty(3, dtype=object)
    vector_gd = np.empty(3, dtype=object)
    
    
    for i in range(len(df_pairgenes)):
        
        vector_dd[0] = df_pairgenes.g1g2_name[i]
        vector_dd[1] = df_pairgenes.discrete_distance[i]
        vector_dd[2] = df_pairgenes.OG_g1g2[i]
        
        vector_gd[0] = df_pairgenes.g1g2_name[i]
        vector_gd[1] = df_pairgenes.gap_distance[i]
        vector_dd[2] = df_pairgenes.OG_g1g2[i]
        
        
        vec_temp_dd = np.empty(len(num_col)-3, dtype = object)
        vec_temp_gd = np.empty(len(num_col)-3, dtype = object)
        
        h = 0
        
        ####### CARGAR LAS ANOTACIONES DE LAS ESPECIES RESTANTES Y EMPEZAR A CONSTRUIR EL DATASET 
        
        for k in range(len(csv)):
            
            if(k != j):
        
                df_sp_track = pd.read_csv(dir + csv[k])
    
                
                df_track = df_sp_track[(df_sp_track['orthogroup_id']== df_pairgenes.OG_g1g2[i][0]) | (df_sp_track['orthogroup_id']== df_pairgenes.OG_g1g2[i][1])]
                
                if (len(df_track) != 0):
                    
                    df_track_distance = distance(df_track)      
                    
                    func = lambda x: set(x) == set(df_pairgenes.OG_g1g2[i])
                    
                    df_track_distance = df_track_distance[df_track_distance['OG_g1g2'].apply(func)]
                    
                    if(len(df_track_distance)!= 0):
                    
                        med_discrete_distance = stat.median(df_track_distance.discrete_distance)                 
                        med_gap_distance = stat.median(df_track_distance.gap_distance)
                        
                        vec_temp_dd[h] = med_discrete_distance
                        vec_temp_gd[h] = med_gap_distance
                        
                h = h+1
                
        
        vec_concat_dd = np.concatenate((vector_dd,vec_temp_dd))
        vec_concat_gd = np.concatenate((vector_gd,vec_temp_gd))
                
                               
        df_pairtrack_dd.loc[len(df_pairtrack_dd)] = vec_concat_dd
        df_pairtrack_gd.loc[len(df_pairtrack_gd)] = vec_concat_gd
        
        
    df_pairtrack_dd.to_hdf('/home/j.falcn/TrackPairwaiseGenes/'+csv[j].replace('.csv','') +'.h5', key='discrete_distance')
    df_pairtrack_gd.to_hdf('/home/j.falcn/TrackPairwaiseGenes/'+csv[j].replace('.csv','') +'.h5', key='gap_distance')


end_time = datetime.datetime.now()  # Detener el contador de tiempo

elapsed_time = end_time - start_time  # Calcular el tiempo transcurrido

days, seconds = elapsed_time.days, elapsed_time.seconds
hours = seconds // 3600
minutes = (seconds % 3600) // 60
seconds = seconds % 60

print(f"Tiempo de ejecucion: {days} días, {hours} horas, {minutes} minutos, {seconds} segundos")
