import networkx as nx
from typing import Dict

def find_leaves_and_depths(tree: nx.Graph, root) -> Dict[int, int]:
    """Finds the leaves of a tree and their depths from the root.
    :param tree: nx.Graph
    :param root: root node
    :return: dict"""
    leaves = [node for node in tree.nodes if tree.degree(node) == 1]
    return_dict = {}
    for leaf in leaves:
        return_dict[leaf] = nx.shortest_path_length(tree, root, leaf)

    return return_dict