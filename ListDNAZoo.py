#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 14:10:54 2023

@author: andreafalcon
"""

import pandas as pd
import pickle
import requests
import json
import urllib
import datetime

start_time = datetime.datetime.now()

species = pd.read_csv('/Users/andreafalcon/Dropbox/PostDoc_NEU/DNAZoo/DNAZoo_Species.txt')

with open('/Users/andreafalcon/Dropbox/PostDoc_NEU/DNAZoo/dnazoo_species_ASM.pkl', 'rb') as f:
    data = pickle.load(f)
    
for i in range(len(data)):
    data[i] = data[i].replace(' ', '_')
    
data1 = []

for i in range(len(species)): 
    data1.append(species.Specie[i].replace('/',''))
    
datafin = list(set(data).intersection(data1))

datafin.append('Balaenoptera_physalus')
datafin.append('Canis_lupus_dingo')
datafin.append('Mustela_putorius_furo')

for i in range(len(datafin)):
    datafin[i] = datafin[i] + '/'
    
datafin[52] = 'Mesocricetus_auratus__MesAur1.0/'
datafin[61] = 'Chinchilla_lanigera/'  

main_url = 'https://dnazoo.s3.wasabisys.com/'

cont_file = 0
cont_json = 0

for j in range(len(datafin)):
    
    url_specie =  main_url + datafin[j]
    
    file_url = url_specie + 'README.json'
    
    response_json = requests.head(file_url)
    
    if(response_json.status_code == 200):
        
        cont_json = cont_json + 1
    
        with urllib.request.urlopen(file_url) as url:
            data = json.load(url)
            
        if (datafin[j] == 'Lutra_lutra/'):
            
            name = 'mLutLut1_HiC'
            
        else:
        
            name = data['chromlengthAssembly']['name']
            
            
        print(datafin[j], name)
    
        gff3_name = name + '.fasta_v2.functional.gff3.gz'
        
        response = requests.head(url_specie + gff3_name)
        
        if (response.status_code == 200):
            
            cont_file = cont_file + 1
            
        else:
            print(datafin[j], j)
            
    else:
        print(datafin[j], j)
        
datafin.sort()
                  

df = pd.DataFrame(datafin, columns = ['Specie'])

df.to_csv('/Users/andreafalcon/Dropbox/PostDoc_NEU/DNAZoo/DNAZoo_SpeciesReduced.csv', encoding='utf-8', index=False, header=True)

end_time = datetime.datetime.now()  # Detener el contador de tiempo

elapsed_time = end_time - start_time  # Calcular el tiempo transcurrido

days, seconds = elapsed_time.days, elapsed_time.seconds
hours = seconds // 3600
minutes = (seconds % 3600) // 60
seconds = seconds % 60

print(f"Total time: {days} days, {hours} hours, {minutes} minutes, {seconds} seconds")
