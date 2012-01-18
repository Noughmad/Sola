from grafi import *

N = int(sys.argv[1])
E = int(sys.argv[2])

step = int(E*N/100)

G = local_graph(N, E)
orig_edges = [e for e in G.edges()]
    
for i in range(E*N):
    if (i % step == 0):
        half, full = sirjenje(G, N, 0.05)
        print ('%g\t%g\t%g\t%g\t%g' % (
            i*1.0/E/N, 
            nx.average_shortest_path_length(G)*math.log(E)/math.log(N), 
            nx.average_clustering(G), 
            half,
            full
            ))
    G.remove_edge(*orig_edges.pop(random.randint(0, len(orig_edges)-1)))
    a, b = add_random_edge(G)
    while not nx.is_connected(G):
        G.remove_edge(a,b)
        a,b = add_random_edge(G)