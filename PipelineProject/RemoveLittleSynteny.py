#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 14:15:17 2023

@author: andreafalcon
"""

import re

# Abrir el archivo TSV
with open('/Users/andreafalcon/Documents/all_dagchainer_results/Eulemur_flavifrons-Acinonyx_jubatus-dagchainer_input_index.tsv.aligncoords', 'r') as f:
    data = f.read()

# Definir la expresión regular para buscar los bloques de información
pattern = r'(?s)##\s*(.*?)\s*##'

# Buscar todos los bloques de información que coincidan con la expresión regular
matches = re.findall(pattern, data)

# Eliminar los bloques de información que sean de longitud menor a 4
new_data = ''
for match in matches:
    if len(match.split('\n')) > 4:
        #print(match.split('\n'))
        new_data = new_data + '##' + match + '\n'

# Guardar el archivo TSV con los bloques de información eliminados
with open('/Users/andreafalcon/Eulemur_flavifrons-Acinonyx_jubatus.tsv', 'w') as f:
    f.write(new_data)