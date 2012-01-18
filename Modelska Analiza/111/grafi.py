import networkx as nx
import sys
import random
import math

def random_edge(N):
    a = random.randint(0,N-1)
    b = a
    while a == b:
        b = random.randint(0,N-1)
    return a,b
    
def add_random_edge(G):
    while True:
        a, b = random_edge(G.number_of_nodes())
        if not G.has_edge(a, b):
            break
    G.add_edge(a,b)
    return a, b
    

def sirjenje(graph, N, p):
    G = graph.copy()
    start_info(G)
    T = 0
    nk = 1
    pol = False
    vsi = False
    while nk < N:
        knowing = [n for n in G if 'i' in G.node[n] ]
        for u in knowing:
            for v in [n for n in G.neighbors_iter(u) if 'i' not in G.node[n] ]:
                if random.random() < p:
                    G.node[v]['i'] = True
                    nk = nk + 1
        T = T + p
        if not pol and nk > N/2:
            pol = True
            yield T
    yield T
    
def start_info(G):
    kn = random.choice(G.nodes())
    G.node[kn]['i'] = True
    
def local_graph(N, E):
    G = nx.Graph()
    G.add_nodes_from(range(N))
    
    for d in range(E):
        G.add_edges_from([ (i,(i+1+d) % N) for i in range(N) ])
    return G
