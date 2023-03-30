#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 12:25:34 2023

@author: andreafalcon
"""

# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.

Andrea Falcon

"""


import pandas as pd
import urllib.request
import json
import requests
import datetime

start_time = datetime.datetime.now()
            

species = pd.read_csv('/Users/andreafalcon/Dropbox/PostDoc_NEU/DNAZoo/DNAZoo_SpeciesReduced.csv')

main_url = 'https://dnazoo.s3.wasabisys.com/'


for j in range(len(species)):
    
    url_specie =  main_url + species.Specie[j]
    
    file_url = url_specie + 'README.json'
    
    response_json = requests.head(file_url)
    
    if(response_json.status_code == 200):
    
        with urllib.request.urlopen(file_url) as url:
            data = json.load(url)
        
        if (species.Specie[j] == 'Lutra_lutra/'):
            
            name = 'mLutLut1_HiC'
            
        else:
        
            name = data['chromlengthAssembly']['name']
    
        gff3_name = name + '.fasta_v2.functional.gff3.gz'
        
        response = requests.head(url_specie + gff3_name)
        
        n_chrom = data['chromlengthAssembly']['karyotype']
        
        if (response.status_code == 200):
            
            print(species.Specie[j], name, n_chrom)


end_time = datetime.datetime.now()  # Detener el contador de tiempo

elapsed_time = end_time - start_time  # Calcular el tiempo transcurrido

days, seconds = elapsed_time.days, elapsed_time.seconds
hours = seconds // 3600
minutes = (seconds % 3600) // 60
seconds = seconds % 60

print(f"Tiempo de ejecucion: {days} días, {hours} horas, {minutes} minutos, {seconds} segundos")