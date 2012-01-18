from grafi import *

N = int(sys.argv[1])
E = int(sys.argv[2])
p = 0.1

step = int(sys.argv[3]) * E

G = local_graph(N,E)
orig_edges = G.edges()
kn = random.choice(G.nodes())
G.node[kn]['i'] = True

def sirjenje(graph, eta):
    G = graph.copy()
    T = 0
    nk = 1
    with open("g_sirjenje_%d_%d.dat" % (N, int(100*eta)), 'w') as f:
        f.write("%g %g\n" % (T, nk/N))
        while nk < N:
            knowing = [n for n in G if 'i' in G.node[n] ]
            for u in knowing:
                for v in [n for n in G.neighbors_iter(u) if 'i' not in G.node[n] ]:
                    if random.random() < p:
                        G.node[v]['i'] = True
                        nk = nk + 1
            T = T + p
            f.write("%g %g\n" % (T, nk/N))

for i in range(E*N):
    if (i % step == 0):
        sirjenje(G, i*1.0/E/N)
    G.remove_edge(*orig_edges.pop(random.randint(0, len(orig_edges)-1)))
    a, b = add_random_edge(G)
    while not nx.is_connected(G):
        G.remove_edge(a,b)
        a,b = add_random_edge(G)
        
sirjenje(G, 1)