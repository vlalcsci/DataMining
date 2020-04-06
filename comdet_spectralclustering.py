# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 19:14:26 2019

@author: Vishal Lal
"""

import networkx as nx
import sys
import numpy as np




def generate_laplacianmatrix(cluster, G):
    N = len(cluster)
    cluster = sorted(cluster)
    
    A = np.zeros((N, N))
    for i in cluster:
        for j in list(G.adj[i]):
            if j in cluster:
                A[cluster.index(i)][cluster.index(j)] = 1
                A[cluster.index(j)][cluster.index(i)] = 1
    D = np.zeros((N, N))
    for i in range(N):
        D[i][i] = G.degree[cluster[i]]

    L = D - A   
    return L

def modify_graph(pos, neg, G):
    for i in pos:
        for j in list(G.adj[i]):
            if j in neg:
                G.remove_edge(i, j)

    return G

if __name__ == "__main__":        
    inputfile=open(sys.argv[1],"r")
    k=int(sys.argv[2])
    outputfile=sys.argv[3]
    clusters=[]
    node_set=set()
    G=nx.Graph()
    for entry in inputfile.readlines(): 
        entry = entry.split()
        node_set.add(int(entry[0]))
        node_set.add(int(entry[1]))
        G.add_node(int(entry[0]))
        G.add_node(int(entry[1]))
        G.add_edge(int(entry[0]), int(entry[1]))
    clusters.append(sorted(list(node_set)))

    for iter in range(k - 1):
        # pick largest cluster 
        len_maxcluster=float('-inf')
        maxcluster=[]
        for cl in clusters:
            if len(cl)>len_maxcluster:
                maxcluster=cl
                len_maxcluster=len(cl)
           
    
        L = generate_laplacianmatrix(maxcluster,G)
        
        clusters.remove(maxcluster)
        eigen_values,eigen_vectors = np.linalg.eig(L)
        
        # pick the second smallest eigenvalue and corresponding eigenvector
        minimum,second_minimum=float('inf'),float('inf')
        for ev in eigen_values:
            if ev<=minimum:
                minimum, min2 = ev, minimum
            elif ev<min2:
                min2=ev

        eigen_index=list(eigen_values).index(min2)
        second_eigenvector = eigen_vectors[:, eigen_index]
        
        # split the eigenvector at 0 into 2 clusters
        positive=[]
        negative = []
        for i in range(len(second_eigenvector)):
            if second_eigenvector[i] >= 0:
                positive.append(maxcluster[i])
            else:
                negative.append(maxcluster[i])
        clusters.append(sorted(positive))
        clusters.append(sorted(negative))
        
        G = modify_graph(positive, negative, G)
        
        

    result=[]
    for cl in clusters:
        cl = sorted(cl)
        
        result.append(",".join(str(i) for i in cl))
    
    output = open(outputfile, "w")
    for i in result:
        output.write(i + "\n")
    output.close()
