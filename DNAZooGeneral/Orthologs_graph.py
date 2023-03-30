#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 15:23:00 2023

@author: andreafalcon
"""

import pandas as pd
import matplotlib.pyplot as plot
import datetime

start_time = datetime.datetime.now()


df_sp1 = pd.read_csv('/Users/andreafalcon/Dropbox/PostDoc_NEU/DNAZoo/GeneAnnotation/Acinonyx_jubatus.csv')
df_sp2 = pd.read_csv('/Users/andreafalcon/Dropbox/PostDoc_NEU/DNAZoo/GeneAnnotation/Acinonyx_jubatus.csv')

ortho_intersection = list(set(df_sp1.orthogroup_id).intersection(set(df_sp2.orthogroup_id)))

df_s1_ortho = df_sp1[df_sp1.orthogroup_id.isin(ortho_intersection)]
df_s2_ortho = df_sp2[df_sp2.orthogroup_id.isin(ortho_intersection)]


df_s1_ortho = df_s1_ortho.reset_index(drop=False)
df_s2_ortho = df_s2_ortho.reset_index(drop=False)

df = pd.merge(df_s1_ortho, df_s2_ortho, on = 'orthogroup_id', how = 'inner')

plot.scatter(df.index_x, df.index_y, c = 'blue', s = 0.00001)


end_time = datetime.datetime.now()  # Detener el contador de tiempo

elapsed_time = end_time - start_time  # Calcular el tiempo transcurrido

days, seconds = elapsed_time.days, elapsed_time.seconds
hours = seconds // 3600
minutes = (seconds % 3600) // 60
seconds = seconds % 60

print(f"Tiempo de ejecucion: {days} días, {hours} horas, {minutes} minutos, {seconds} segundos")
