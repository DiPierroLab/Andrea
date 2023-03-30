#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 12:51:49 2023

@author: andreafalcon
"""

# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.

Andrea Falcon

"""


import gffpandas.gffpandas as gffpd
import pandas as pd
import numpy as np
import urllib.request
import json
import gzip
import shutil
import requests
import os
import datetime

start_time = datetime.datetime.now()

with open('/home/j.falcn/Orthogroups.txt', 'r') as archivo:
    
    orthogroups = {}
    
    for line in archivo:

        part = line.strip().split(':')
        
        if len(part) == 2:

            genes = part[1].strip().split()
            genes = list(map(lambda x: x[:-3], genes))
            orthogroups[part[0]] = genes
            

species = pd.read_csv('/home/j.falcn/DNAZoo_SpeciesReduced.csv')   #there are 103 species in this file

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
        
        if (response.status_code == 200):
            
            print(j, species.Specie[j])
            
            urllib.request.urlretrieve(url_specie + gff3_name, '/home/j.falcn/' + gff3_name)
            
            with gzip.open('/home/j.falcn/' + gff3_name, 'rb') as f_in:
                with open('/home/j.falcn/' + gff3_name.replace('.gz',''), 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
                    
            os.remove('/home/j.falcn/' + gff3_name)
    
            annotation = gffpd.read_gff3('/home/j.falcn/' + gff3_name.replace('.gz',''))
            
            df_annotation = annotation._read_gff3_to_df()
            
            del annotation
            
            names = list(df_annotation.columns.values)
            
            df_gene = df_annotation.loc[df_annotation['type'] == 'gene']
            
            del df_annotation
            
            df_gene = df_gene.reset_index(drop = True)
            
            df_work= df_gene[['seq_id','start','end','strand','attributes']]
            
            del df_gene
            
            df_work['gene_id'] = np.nan
            df_work['gene_name'] = np.nan
            
            
            for i in range(len(df_work)):
                df_work['gene_id'][i] = df_work['attributes'][i].split(';')[0].split('=')[1]
                df_work['gene_name'][i] = df_work['attributes'][i].split(';')[1].split('=')[1]
                
            df_work = df_work[['seq_id','start','end','strand','gene_id','gene_name','attributes']]
            
            unique_scaffold = df_work.seq_id.unique()
            
            n_scaffold = np.empty(len(unique_scaffold), dtype=object)
            
            for i in range(len(unique_scaffold)):
                n_scaffold[i] = (df_work['seq_id'].isin([unique_scaffold[i]])==True).sum()
                
            d = {'seq_id': unique_scaffold, 'size_scaffold': n_scaffold}
            
            df_scaffold = pd.DataFrame(data=d)
            
            df = df_work.merge(df_scaffold, on='seq_id', how='left')
            
            del df_work
            
            df = df.sort_values(by='size_scaffold', ascending=False)
            
            
            if(species.Specie[j] == 'Hippotragus_niger/'):
                
                n_chrom = 30
            
            elif(species.Specie[j] == 'Struthio_camelus/'):
                
                n_chrom * 40
            
            else:
                n_chrom = int(data['chromlengthAssembly']['karyotype'].split('=')[1])/2
            
            array = ['HiC_scaffold_' + str(i) for i in range(1, int(n_chrom+1))]
            
            df_list = []
            
            for i in range(len(array)):
            
                df_list.append(df.loc[df["seq_id"] == array[i]])
                df_list[i] = df_list[i].sort_values(by = 'start', ascending = True)
                df_list[i] = df_list[i].reset_index(drop = True)
                
            df_chrom = pd.concat(df_list, ignore_index=True)
            
            del df
            del df_list
            
            
            df_chrom = df_chrom[['seq_id','start','end','strand','gene_id','gene_name','attributes']]
                        
            df_chrom['orthogroup_id'] = np.nan        
                        
            for k in range(len(df_chrom)):
                gen = df_chrom.gene_name[k]
                
                for clave,valor in orthogroups.items():
                    if gen in valor:
                        if clave == 'OG0000000':
                            df_chrom['orthogroup_id'][k] = '0'
                        else:
                            df_chrom['orthogroup_id'][k] = clave.lstrip('OG0')
                        break
            
            df_chrom = df_chrom[['seq_id','start','end','strand','gene_id','gene_name','orthogroup_id','attributes']]
            
            binomial = data['organism']['binomial']
            
            df_chrom.to_csv('/home/j.falcn/GeneAnnotation/' + binomial.replace(' ','_') + '.csv', encoding='utf-8', index=False, header=True)
            
            os.remove('/home/j.falcn/' + gff3_name.replace('.gz',''))


end_time = datetime.datetime.now()  # Detener el contador de tiempo

elapsed_time = end_time - start_time  # Calcular el tiempo transcurrido

days, seconds = elapsed_time.days, elapsed_time.seconds
hours = seconds // 3600
minutes = (seconds % 3600) // 60
seconds = seconds % 60

print(f"Tiempo de ejecucion: {days} días, {hours} horas, {minutes} minutos, {seconds} segundos")
