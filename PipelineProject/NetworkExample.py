#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 31 15:56:22 2023

@author: andreafalcon
"""

import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import igraph as ig
import itertools
import numpy as np

dir = '/Users/andreafalcon/Documents/NetworkFiles-dragchainer/'

df1 = pd.read_csv(dir + 'Bison_bison-Acinonyx_jubatus.csv')
df2 = pd.read_csv(dir + 'Bison_bison-Eulemur_flavifrons.csv')
df3 = pd.read_csv(dir + 'Eulemur_flavifrons-Acinonyx_jubatus.csv')

df = pd.concat([df1, df2, df2], axis = 0, ignore_index= True)

df = df.drop_duplicates()

#df.to_csv('/Users/andreafalcon/Documents/net_example.csv',  encoding='utf-8', index=False, header=True)

G = nx.from_pandas_edgelist(df, 'GeneA', 'GeneB')

for nodo in G.nodes():
    
    G.nodes[nodo]['type'] = nodo.split('_')[0]
    G.nodes[nodo]['degree'] = G.degree(nodo)
    G.nodes[nodo]['name'] = nodo


component_sizes = [len(c) for c in nx.connected_components(G)]
component_unique = list(set(component_sizes))
num_cc = [component_sizes.count(u) for u in component_unique]
frac_cc = [x/len(component_sizes) for x in num_cc]

plt.plot(component_unique, num_cc, marker='o', c='dodgerblue')
plt.ylabel('Number')
plt.xlabel('CC Size')
plt.yscale('log')
plt.title('')
plt.legend(labels=[])
plt.show()


cc = nx.connected_components(G)
cc_filt=[]
for c in cc:
    if(len(c)==5):        
        cc_filt.append(list(c))
        
subgraph = G.subgraph(list(itertools.chain.from_iterable(cc_filt[100:300])))
df_net = nx.to_pandas_edgelist(subgraph, source='source', target='target')

sub_igraph = ig.Graph.from_networkx(subgraph)

n_tipos = len(set(sub_igraph.vs['type']))
colormap = plt.cm.get_cmap('viridis', n_tipos)  
tipos = sorted(set(sub_igraph.vs['type']))  
colores = [colormap(i) for i in np.linspace(0, 1, n_tipos)] 
colores_dict = {t: c for t, c in zip(tipos, colores)} 
colors = [colores_dict[tipo] for tipo in sub_igraph.vs['type']]
sizes = [5*grado for grado in sub_igraph.vs['degree']]

ig.plot(sub_igraph, layout = sub_igraph.layout_fruchterman_reingold(), vertex_size = sizes, vertex_color = colors)

keys = list(set(sub_igraph.vs['degree']))
keys = [str(i) for i in keys]
diccionario = dict.fromkeys(keys)

for key in list(diccionario.keys()):
    deg = [elemento for elemento in sub_igraph.vs['name'] if subgraph.nodes[elemento]['degree']==int(key)]
    add = []
    for type in tipos:
        add.append(len([elemento for elemento in deg if subgraph.nodes[elemento]['type'] == type ]))
    diccionario[key] = add
    
df_degree =  pd.DataFrame.from_dict(diccionario, orient='index', columns=tipos)

n = len(df_degree.index)
x = np.arange(n)
width = 0.25
plt.bar(x - width, df_degree.Acinonyx, width=width, label='Cheeta', color = colores[0])
plt.bar(x, df_degree.Bison, width=width, label='Bison', color = colores[1])
plt.bar(x + width, df_degree.Eulemur, width=width, label='Lemur', color = colores[2])
plt.xticks(x, df_degree.index)
plt.legend(loc='best')
plt.show()

g = ig.Graph.from_networkx(G)
infomap = g.community_infomap()
communities_unique = list(set(infomap.sizes()))
num_com = [infomap.sizes().count(u) for u in communities_unique]
frac_com = [x/len(infomap.sizes()) for x in num_com]

plt.plot(communities_unique, num_com, marker='o', c='lightgreen')
plt.ylabel('Number')
plt.xlabel('Community Size')
plt.yscale('log')
plt.title('')
plt.legend(labels=[])
plt.show()
