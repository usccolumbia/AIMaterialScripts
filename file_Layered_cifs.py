#!/usr/bin/env python
# coding: utf-8

# In[1]:


import warnings
warnings.simplefilter('ignore')


# In[2]:


# chemically connected groups are constructed
class Graph:
 
    # init function to declare class variables
    def __init__(self, V):
        self.V = V
        self.adj = [[] for i in range(V)]
 
    def DFSUtil(self, temp, v, visited):
 
        # Mark the current vertex as visited
        visited[v] = True
 
        # Store the vertex to list
        temp.append(v)
 
        # Repeat for all vertices adjacent
        # to this vertex v
        for i in self.adj[v]:
            if visited[i] == False:
 
                # Update the list
                temp = self.DFSUtil(temp, i, visited)
        return temp
 
    # method to add an undirected edge
    def addEdge(self, v, w):
        self.adj[v].append(w)
        self.adj[w].append(v)
 
    # Method to retrieve connected components
    # in an undirected graph
    def connectedComponents(self):
        visited = []
        cc = []
        for i in range(self.V):
            visited.append(False)
        for v in range(self.V):
            if visited[v] == False:
                temp = []
                cc.append(self.DFSUtil(temp, v, visited))
        return cc


# In[ ]:


import glob
import numpy as np
from pymatgen.io.cif import CifParser
from pymatgen.analysis.molecule_structure_comparator import CovalentRadius
from pymatgen.analysis.local_env import MinimumOKeeffeNN

covalent_radius_dict = CovalentRadius().radius
maybe = []
voronoi = MinimumOKeeffeNN()
start = 0
index = 0
files = glob.glob('icsd_cif/*.cif')

for name in range(start, len(files)):
    name = files[index]
    
    print(index)
    
    try: 
        parser = CifParser(f"{name}")
        structure = parser.get_structures()[0]
        s = structure


        s.make_supercell([3,3,3])

        n = len(s.sites)
        M = np.zeros((n,n))

        for i in range(len(s.sites)):
            site = s.sites[i]

            for neighbor in voronoi.get_nn_info(s, i):
                j = neighbor['site_index']
                other = neighbor['site']

                r1 = covalent_radius_dict[site.specie.symbol]
                r2 = covalent_radius_dict[other.specie.symbol]

                if site == other:
                    continue

                if site.distance(other) < (r1 + r2) * 1.2:
                    M[i,j] = 1
                else:
                    # print("vdw")
                    pass

        print(index,name)

        g = Graph(n)
        for i in range(n):
            for j in range(n): 
                if M[i,j] == 1:
                    g.addEdge(i, j)
        cc = g.connectedComponents()

        i = 0
        while i < len(cc):
            if len(cc[i]) < 3:
                cc.pop(i)
            else:
                i += 1

        shortest_distance = 9223372036854775807
        for r in range(len(cc)):
            for r2 in range(r + 1,len(cc)):
                for c in cc[r]:
                    for c2 in cc[r2]:
                        site = s.sites[c]
                        other = s.sites[c2]
                        if other.distance(site) < shortest_distance:
                            shortest_distance = other.distance(site)

        if(shortest_distance < 10 and shortest_distance > 3):
            maybe.append([name, shortest_distance])
            print("new", maybe)
    
        print()
    except:
        print(name, 'failed')
        print(maybe)
        print()
        
    index += 1


# In[ ]:


print(maybe)


# In[ ]:




