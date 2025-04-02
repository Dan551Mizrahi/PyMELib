import networkx as nx
from PyMELib.labels import F_sigma, F_omega, F_rho

MAX_CHR = 1114111

def read_hypergraph(path_to_file: str) -> nx.Graph:
    """Reads a hypergraph from a file.
    Outputs the incidence graph, with proper constraints.
    :param path_to_file: str
    :return: nx.Graph"""
    with open(path_to_file, 'r') as f:
        lines = f.readlines()

        hyperedges = []
        vertices = set()
        for line in lines:
            hyperedges.append(tuple([int(x) for x in line.replace(',', ' ').strip().split()]))
            vertices.update(hyperedges[-1])

        max_vertex = max(vertices)

        G = nx.Graph()
        G.add_nodes_from(vertices)

        G.add_node(MAX_CHR)
        G.nodes[MAX_CHR]['options'] = [F_sigma.SI, F_sigma.S1, F_sigma.S0]

        for node in vertices:
            G.add_edge(MAX_CHR, node)
            G.nodes[node]['options'] = [F_sigma.S0, F_sigma.S1, F_omega.W0, F_omega.W1]

        for i, hyperedge in enumerate(hyperedges):
            hyperedge_id = i+max_vertex+1
            G.add_node(hyperedge_id)
            for vertex in hyperedge:
                G.add_edge(vertex, hyperedge_id)
            G.nodes[hyperedge_id]['options'] = [F_rho.R0, F_rho.R1, F_rho.R2, F_omega.W0, F_omega.W1]

    return G