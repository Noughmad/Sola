from grafi import *

N = int(sys.argv[1])
E = int(sys.argv[2])
p = 0.1

step = int(sys.argv[3])
limit = int(sys.argv[4])

G = local_graph(N,E)
orig_edges = [e for e in G.edges()]

with open('g_hitrost_%d.dat' % N, 'w') as f:
    for i in range(limit):
        if (i % step == 0):
            sirjenje(G, i*1.0/E/N, f)
        G.remove_edge(*orig_edges.pop(random.randint(0, len(orig_edges)-1)))
        a, b = add_random_edge(G)
        while not nx.is_connected(G):
            G.remove_edge(a,b)
            a,b = add_random_edge(G)
    sirjenje(G, limit*1.0/E/N, f)
