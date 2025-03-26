import networkx as nx
from typing import Iterable
from PyMELib.labels import F_sigma, F_omega, F_rho

def add_constraints_on_graph(G: nx.Graph, include_in_ds: Iterable = [], exclude_from_ds: Iterable = []) -> int:
    """
    This function is a helper function to add constraints on a graph, before the running of the preprocessing phase or the enumeration phase.
    There are two kinds of constraints: include_in_ds and exclude_from_ds.
    After using this function, remember to call the preprocessing phase and the enumeration phase with options_for_labels=True.
    :param G: The graph on which the constraints will be added.
    :param include_in_ds: The nodes that you want to force to be included in the minimal dominating set (or hitting set).
    :param exclude_from_ds: The nodes that you want to force to be excluded from the minimal dominating set (or hitting set).
    :return: 0 if successful, -1 if you asked for invalid constraints.
    """

    for node in include_in_ds:
        if 'options' not in G.nodes[node]:
            G.nodes[node]['options'] = [F_sigma.SI, F_sigma.S1, F_sigma.S0]
        else:
            intersection = set(G.nodes[node]['options']).intersection({F_sigma.SI, F_sigma.S1, F_sigma.S0})
            if len(intersection) == 0:
                print("You asked for invalid constraints. The node", node, "cannot be forced to be in the dominating set.")
                return -1
            else:
                G.nodes[node]['options'] = list(intersection)

    for node in exclude_from_ds:
        if 'options' not in G.nodes[node]:
            G.nodes[node]['options'] = [F_rho.R0, F_rho.R1, F_rho.R2, F_omega.W0, F_omega.W1]
        else:
            intersection = set(G.nodes[node]['options']).intersection({F_rho.R0, F_rho.R1, F_rho.R2, F_omega.W0, F_omega.W1})
            if len(intersection) == 0:
                print("You asked for invalid constraints. The node", node, "cannot be forced to be excluded from the dominating set.")
                return -1
            else:
                G.nodes[node]['options'] = list(intersection)

    for node in G.nodes:
        if 'options' not in G.nodes[node]:
            G.nodes[node]['options'] = [F_sigma.SI, F_sigma.S1, F_sigma.S0, F_omega.W0, F_omega.W1, F_rho.R0, F_rho.R1, F_rho.R2]

    return 0
