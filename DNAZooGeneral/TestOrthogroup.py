#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 15:20:51 2023

@author: andreafalcon
"""

import pandas as pd
import os

dir = '/Users/andreafalcon/Documents/GeneAnnotation/'

contenido = os.listdir(dir)

csv = []

for fichero in contenido:
    if os.path.isfile(os.path.join(dir, fichero)) and fichero.endswith('.csv'):
        csv.append(fichero)
        
csv.sort()
        
og_empty = []

d = {'Specie':[]}
        
for j in range(len(csv)):

    df_sp = pd.read_csv(dir + csv[j])
    
    a = df_sp.loc[df_sp.isna().any(axis=1)]["orthogroup_id"]
    
    if(len(a) != 0 ):
        
        og_empty.append(csv[j])
        
    else:
        
        d['Specie'].append(csv[j])
        
species = pd.DataFrame(d)

species.to_csv('/Users/andreafalcon/Dropbox/PostDoc_NEU/DNAZoo/GeneAnnotation_SpeciesOG.csv', encoding='utf-8', index=False, header=True)
        
