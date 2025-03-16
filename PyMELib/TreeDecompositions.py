import networkx as nx
from typing import Dict, List
from Factors import MemoTable
from itertools import product
from frozendict import frozendict

# This code was created by Dan Mizrahi with the help of GitHub Copilot and Gemini
# The code supports graphs up to 1,114,112 nodes (maximum chr value in python).

def generate_dictionaries_from_sets(data: Dict[str, List]) -> List:
    # Extract keys and sets of values
    keys = list(data.keys())
    values = list(data.values())

    # Generate Cartesian product of the sets of values
    combinations = product(*values)

    # Create a list of dictionaries for each combination
    result = [dict(zip(keys, combination)) for combination in combinations]

    return result

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

# Possible labels for the vertices
SI = 0
S0 = 1
S1 = 2
W0 = 3
W1 = 4
R0 = 5
R1 = 6
R2 = 7
F_sigma = {SI, S0, S1}
F_omega = {W0, W1}
F_rho = {R0, R1, R2}
F = F_sigma.union(F_omega.union(F_rho))

def in_sigma(label: int):
    return label < 3

def in_omega(label: int):
    return (label < 5  and label > 2)

def in_rho(label: int):
    return label > 4

def same_label_group(l0: int, l1: int):
    return ((in_sigma(l0) and in_sigma(l1)) or
            (in_omega(l0) and in_omega(l1)) or
            (in_rho(l0) and in_rho(l1)))

def join_labels(l0: int, l1: int):

    flag = True
    label_original_vertex = None

    sum_of_labels = l0 + l1

    if in_sigma(l0) and in_sigma(l1):
        if sum_of_labels == 0:
            label_original_vertex = SI
        elif sum_of_labels == 1:
            flag = False
        elif sum_of_labels == 2:
            if (l0 == SI and l1 == S1) or (l0 == S1 and l1 == SI):
                label_original_vertex = S1
            else:
                label_original_vertex = S0
        else:
            label_original_vertex = S1
    elif in_omega(l0) and in_omega(l1):
        if sum_of_labels > 7:
            flag = False
        elif sum_of_labels == 7:
            label_original_vertex = W1
        else:
            label_original_vertex = W0
    elif in_rho(l0) and in_rho(l1):
        if sum_of_labels == 10:
            label_original_vertex = R0
        elif sum_of_labels == 11:
            label_original_vertex = R1
        else:
            label_original_vertex = R2
    elif in_rho(l0) and l1 == W0:
        label_original_vertex = l0
    elif l0 == W0  and in_rho(l1):
        label_original_vertex = l1
    else:
        flag = False

    return flag, label_original_vertex


trans_dict = {"SI":SI,
              "S0":S0,
              "S1":S1,
              "W0":W0,
              "W1":W1,
              "R0":R0,
              "R1":R1,
              "R2":R2}


def V_label_range(theta: Dict[str, int], label_lower: int, label_upper: int):
    return_set = set()
    for vertex, label in theta.items():
        if label_lower <= label <= label_upper:
            return_set.add(vertex)
    return return_set


def V_label(label_a: str, theta: Dict[str, int]):
    if label_a == "S":
        return V_label_range(theta, SI, S1)
    elif label_a == "W":
        return V_label_range(theta, W0, W1)
    elif label_a == "R":
        return V_label_range(theta, R0, R2)
    else:
        return V_label_range(theta, trans_dict[label_a], trans_dict[label_a])

def V_label_S_W(theta: Dict[str, int]):
    return_set_S = set()
    return_set_W = set()
    for vertex, label in theta.items():
        if label < W0:
            return_set_S.add(vertex)
        elif label < R0:
            return_set_W.add(vertex)
    return return_set_S, return_set_W


LEAF = 0
INTRODUCE = 1
FORGET = 2
JOIN = 3
JOIN_INTRODUCE = 4
BIG_JOIN_INTRODUCE = 5
JOIN_FORGET = 6
ROOT = 7


class RootedTreeDecomposition(nx.classes.digraph.DiGraph):
    """
    This class provides functionality to generate a rooted tree decomposition of a given graph.
    The decomposition is based on the junction tree of the graph and allows for subsequent operations and analysis.
    """

    def __init__(self, G: nx.classes.graph.Graph, root: tuple = tuple(), root_heuristic="leaf", first_label = 0, *args, **kwargs):
        """
        Initializes the RootedTreeDecomposition object.

        Args:
            G: The input graph (NetworkX Graph object) to be decomposed.
            root: (Optional) (tuple) The root node of the decomposition. If not provided, a root is chosen from the junction tree.
            *args, **kwargs: Additional arguments passed to the parent DiGraph class.
        """

        super().__init__(**kwargs)

        # self.original_graph = nx.convert_node_labels_to_integers(G, label_attribute="original_name", first_label=first_label)
        self.original_graph = G
        for node in self.original_graph.nodes:
            self.original_graph.nodes[node]["original_name"] = node

        # TODO: junction tree heuristic handling
        T = nx.junction_tree(self.original_graph)

        root_flag = root == tuple()

        if root_flag:
            root = next(iter(T.nodes))

        # Root the junction tree (initially)
        self.bfs_tree = nx.bfs_tree(T, root)

        if root_heuristic == "leaf" and not root_flag:
            leaves = [node for node in T.nodes if T.degree(node) == 1]
            if len(leaves) != 0:
                root = leaves[0]
                self.bfs_tree = nx.bfs_tree(T, root)

        # Some manipulation on the nodes of this tree decomposition
        new_nodes = [(i, {"bag": set(vertex)})
                     for i, vertex in enumerate(nx.dfs_postorder_nodes(self.bfs_tree, root))]

        self.new_nodes_dict = {t[0]: t[1]["bag"] for t in new_nodes}

        # Adding the post-manipulation nodes to the tree.
        self.add_nodes_from(new_nodes)

        reversed_dict = {frozenset(v): k for k, v in self.new_nodes_dict.items()}

        # Adding directed edges to the tree.
        for edge in self.bfs_tree.edges:
            self.add_edge(reversed_dict[frozenset(edge[0])], reversed_dict[frozenset(edge[1])])

        self.width = max([len(node[1]) for node in self.nodes(data="bag")]) - 1
        self.root = reversed_dict[frozenset(root)]
        self.original_root = self.root
        self.max_id = max(self.nodes)
        #self.draw_regular()

    def get_original_root(self):
        return self.original_root

    def set_root(self, root: int) -> None:
        """
        Sets the root node of the tree decomposition.

        Args:
            root: (int) The new root node identifier.
        """
        self.root = root

    def get_root(self) -> int:
        """
        Gets the root node of the tree decomposition.

        Returns:
            int: The identifier of the root node.
        """
        return self.root

    def add_node_bag(self, bag_of_node: set) -> int:
        """
        Adds a new node to the tree decomposition with the specified bag of vertices.

        Args:
            bag_of_node: A set containing the vertices (or copies) to be included in the bag.

        Returns:
            int: The identifier (ID) of the newly added node.
        """
        new_node = self.max_id + 1
        self.max_id += 1
        self.add_node(new_node)

        self.nodes[new_node]["bag"] = bag_of_node
        self.new_nodes_dict[new_node] = bag_of_node

        return new_node


class RootedNiceTreeDecomposition(RootedTreeDecomposition):
    """
    This class provides functionality to generate a rooted nice tree decomposition of a given graph.
    """
    def __init__(self, G: nx.classes.graph, root=tuple(), semi_nice=True, *args, **kwargs):
        """
        Initializes the RootedNiceTreeDecomposition object.
        :param G: The input graph (NetworkX Graph object) to be decomposed.
        :param root: The root node of the decomposition.
        """
        super().__init__(G, root, *args, **kwargs)

        new_root = self.add_node_bag(set())
        self.add_edge(new_root, self.get_root())
        self.set_root(new_root)
        if semi_nice:
            self.transform_to_semi_nice_rec(self.get_root())
        else:
            self.transform_to_nice_rec(self.get_root())
        self.nodes[new_root]["type"] = ROOT

        # create a complete order of the vertices (for enumeration process)
        self.Q = []
        self.create_Q(self.get_root())

        #self.draw_nice()

    def create_Q(self, current_node):
        if self.nodes[current_node]["type"] == LEAF:
            return
        if self.nodes[current_node]["type"] == FORGET:
            v = self.nodes[list(self.successors(current_node))[0]]["bag"].difference(
                self.nodes[current_node]["bag"]).pop()
            self.Q.append(v)
        elif self.nodes[current_node]["type"] != ROOT and self.nodes[list(self.predecessors(current_node))[0]][
            "type"] == ROOT:
            self.Q.append(self.nodes[current_node]["bag"][0])
        for child in self.successors(current_node):
            self.create_Q(child)

    def transform_to_nice_rec(self, current_node):
        """
        Recursive function that constructs nice form tree decomposition (instead of the existing tree).
        :param current_node: The current node that we are on TD.
        """
        bag_of_node = self.nodes[current_node]["bag"]
        children = list(self.successors(current_node))
        num_of_children = len(children)

        # Leaf node
        if num_of_children == 0:
            if len(bag_of_node) != 0:
                new_node = self.add_node_bag(set())
                self.add_edge(current_node, new_node)
                self.transform_to_nice_rec(current_node)
            else:
                self.nodes[current_node]["type"] = LEAF

        elif num_of_children == 1:

            child = children[0]
            bag_of_child = self.nodes[child]["bag"]

            diff1 = bag_of_node.difference(bag_of_child)
            diff2 = bag_of_child.difference(bag_of_node)

            # Introduce node
            if len(diff1) > 1 or (len(diff1) == 1 and len(diff2) >= 1):
                # creates a new Introduce node
                new_node_bag = bag_of_node.diffrence({diff1.pop()})
                new_node = self.add_node_bag(new_node_bag)

                self.add_edge(current_node, new_node)
                self.add_edge(new_node, child)
                self.remove_edge(current_node, child)

                self.nodes[current_node]["type"] = INTRODUCE
                self.transform_to_nice_rec(new_node)

            elif len(diff1) == 1 and len(diff2) == 0:
                self.nodes[current_node]["type"] = INTRODUCE
                self.transform_to_nice_rec(child)

            # Forget node
            elif len(diff2) > 1:
                # creates a Forget node
                new_node_bag = bag_of_node.union({diff2.pop()})
                new_node = self.add_node_bag(new_node_bag)

                self.add_edge(current_node, new_node)
                self.add_edge(new_node, child)
                self.remove_edge(current_node, child)

                self.nodes[current_node]["type"] = FORGET
                self.transform_to_nice_rec(new_node)

            elif len(diff1) == 0 and len(diff2) == 1:
                self.nodes[current_node]["type"] = FORGET
                self.transform_to_nice_rec(child)

            else:
                # print("Warning: same two bags one after another. (not as in join node)")
                parent = next(iter(self.predecessors(current_node)))
                self.add_edge(parent, child)
                self.remove_edge(current_node, child)
                self.remove_node(current_node)
                self.transform_to_nice_rec(child)

        # multiple children
        else:

            # remove redundancy inside join nodes by creating introduce nodes
            vertices_in_children = set()
            for child in children:
                vertices_in_children = vertices_in_children.union(self.nodes[child]["bag"])

            redundant_vertices = list(bag_of_node.difference(vertices_in_children))
            essential_vertices = list(bag_of_node.difference(redundant_vertices))

            # create introduce nodes for the redundant vertices if needed [by recursion]
            if len(redundant_vertices) > 0:
                # create the new join node
                new_node_bag = set(essential_vertices)
                new_node = self.add_node_bag(new_node_bag)
                for child in children:
                    self.add_edge(new_node, child)
                    self.remove_edge(current_node, child)
                self.add_edge(current_node, new_node)
                self.transform_to_nice_rec(current_node)
            else:
                # Join node
                self.nodes[current_node]["type"] = JOIN
                child_1 = children[0]

                new_node_1 = self.add_node_bag(self.nodes[current_node]["bag"])

                self.add_edge(current_node, new_node_1)
                self.add_edge(new_node_1, child_1)
                self.remove_edge(current_node, child_1)

                self.transform_to_nice_rec(new_node_1)

                new_node_2 = self.add_node_bag(self.nodes[current_node]["bag"])
                self.add_edge(current_node, new_node_2)

                for child in children[1:]:
                    self.add_edge(new_node_2, child)
                    self.remove_edge(current_node, child)

                self.transform_to_nice_rec(new_node_2)

    def transform_to_semi_nice_rec(self, current_node):
        """
        Recursive function that constructs nice form tree decomposition (instead of the existing tree).
        :param current_node: The current node that we are on TD.
        """
        bag_of_node = self.nodes[current_node]["bag"]
        children = list(self.successors(current_node))
        num_of_children = len(children)

        # Leaf node
        if num_of_children == 0:
            if len(bag_of_node) != 0:
                new_node = self.add_node_bag(set())
                self.add_edge(current_node, new_node)
                self.transform_to_semi_nice_rec(current_node)
            else:
                self.nodes[current_node]["type"] = LEAF

        elif num_of_children == 1:

            child = children[0]
            bag_of_child = self.nodes[child]["bag"]

            diff1 = bag_of_node.difference(bag_of_child)
            diff2 = bag_of_child.difference(bag_of_node)

            # Introduce node
            if len(diff1) > 1 or (len(diff1) == 1 and len(diff2) >= 1):
                # creates a new Introduce node
                new_node_bag = bag_of_node.difference({diff1.pop()})
                new_node = self.add_node_bag(new_node_bag)

                self.add_edge(current_node, new_node)
                self.add_edge(new_node, child)
                self.remove_edge(current_node, child)

                self.nodes[current_node]["type"] = INTRODUCE
                self.transform_to_semi_nice_rec(new_node)

            elif len(diff1) == 1 and len(diff2) == 0:
                self.nodes[current_node]["type"] = INTRODUCE
                self.transform_to_semi_nice_rec(child)

            # Forget node
            elif len(diff2) > 1:
                # creates a Forget node
                new_node_bag = bag_of_node.union({diff2.pop()})
                new_node = self.add_node_bag(new_node_bag)

                self.add_edge(current_node, new_node)
                self.add_edge(new_node, child)
                self.remove_edge(current_node, child)

                self.nodes[current_node]["type"] = FORGET
                self.transform_to_semi_nice_rec(new_node)

            elif len(diff1) == 0 and len(diff2) == 1:
                self.nodes[current_node]["type"] = FORGET
                self.transform_to_semi_nice_rec(child)

            else:
                #print("Warning: same two bags one after another. (not as in join node)")
                parent = next(iter(self.predecessors(current_node)))
                self.add_edge(parent, child)
                self.remove_edge(current_node, child)
                self.remove_node(current_node)
                self.transform_to_semi_nice_rec(child)

        # multiple children
        else:

            # remove redundancy inside join nodes by creating introduce nodes
            vertices_in_children = set()
            for child in children:
                vertices_in_children = vertices_in_children.union(self.nodes[child]["bag"])

            redundant_vertices = list(bag_of_node.difference(vertices_in_children))
            essential_vertices = list(bag_of_node.difference(redundant_vertices))

            # create introduce nodes for the redundant vertices if needed [by recursion]
            if len(redundant_vertices) > 0:
                # create the new join node
                new_node_bag = set(essential_vertices)
                new_node = self.add_node_bag(new_node_bag)
                for child in children:
                    self.add_edge(new_node, child)
                    self.remove_edge(current_node, child)
                self.add_edge(current_node, new_node)
                self.transform_to_semi_nice_rec(current_node)
            else:
                # Join node
                self.nodes[current_node]["type"] = JOIN
                if len(children) > 2:
                    # we want to make our tree binary
                    new_node = self.add_node_bag(self.nodes[current_node]["bag"])
                    self.add_edge(current_node, new_node)
                    self.add_edge(new_node, children[0])
                    self.remove_edge(current_node, children[0])
                    self.add_edge(new_node, children[1])
                    self.remove_edge(current_node, children[1])
                    self.transform_to_semi_nice_rec(current_node)
                else:
                    for child in children:
                        new_node = self.add_node_bag(self.nodes[current_node]["bag"])
                        self.add_edge(current_node, new_node)
                        self.add_edge(new_node, child)
                        self.remove_edge(current_node, child)
                        self.transform_to_semi_nice_rec(new_node)


class RootedDisjointBranchNiceTreeDecomposition(RootedNiceTreeDecomposition):

    def __init__(self, G: nx.classes.graph, root: tuple = tuple(), semi_dntd = True, *args, **kwargs):
        super().__init__(G, root, semi_nice = semi_dntd, *args, **kwargs)


        self.first_appear = {vertex: None for vertex in self.original_graph.nodes}
        if semi_dntd:
            self.semi_ntd_to_semi_dntd(self.get_root(), debug_flag=False)
        else:
            self.ntd_to_dntd(self.get_root())
        self.all_vertices = {v for node in self.nodes for v in self.nodes[node]["bag"]}
        self.local_neighbors(self.get_root())
        self.create_factors()
        self.Q = []
        self.create_Q(self.get_root())
        self.first_appear_update(self.get_root())
        self.trans = {vertex: None for vertex in self.original_graph.nodes}
        # self.draw_nice()

    def get_number_of_join_nodes(self):
        return len([node for node in self.nodes if self.nodes[node]["type"] == JOIN])

    def first_appear_update(self, current_node):
        if self.nodes[current_node]["type"] == LEAF:
            return
        for vertex in self.nodes[current_node]["bag"]:
            if vertex not in self.first_appear.keys() or self.first_appear[vertex] is None:
                self.first_appear[vertex] = current_node
        for child in self.successors(current_node):
            self.first_appear_update(child)

    def create_Q(self, current_node):
        if self.nodes[current_node]["type"] == LEAF:
            return
        if self.nodes[current_node]["type"] != ROOT and self.nodes[list(self.predecessors(current_node))[0]][
            "type"] == ROOT:
            self.Q.append(list(self.nodes[current_node]["bag"])[0])
        if self.nodes[current_node]["type"] == FORGET or self.nodes[current_node]["type"] == JOIN_FORGET:
            v = self.nodes[list(self.successors(current_node))[0]]["bag"].difference(
                self.nodes[current_node]["bag"]).pop()
            self.Q.append(v)

        for child in self.successors(current_node):
            self.create_Q(child)

    def local_neighbors(self, current_node):

        # This function is tentative and should be changed in appropriate way to the conjunctions of factors
        self.nodes[current_node]["local_neighbors"] = dict()

        if self.nodes[current_node]["type"] == LEAF:
            return

        else:
            children = list(self.successors(current_node))
            for child in children:
                self.local_neighbors(child)

            if self.nodes[current_node]["type"] == INTRODUCE:
                child_bag = self.nodes[children[0]]["bag"]
                child_bag = {v[0] for v in child_bag}
                v = self.nodes[current_node]["bag"].difference(self.nodes[children[0]]["bag"]).pop()
                for vertex in self.nodes[current_node]["bag"]:
                    if vertex == v:
                        self.nodes[current_node]["local_neighbors"][vertex] = \
                            {chr(n) for n in self.original_graph.neighbors(ord(v[0]))}.intersection(child_bag)
                    elif ord(v[0]) in self.original_graph.neighbors(ord(vertex[0])):
                        self.nodes[current_node]["local_neighbors"][vertex] = \
                            self.nodes[children[0]]["local_neighbors"][vertex].union({v[0]})
                    else:
                        self.nodes[current_node]["local_neighbors"][vertex] = \
                            self.nodes[children[0]]["local_neighbors"][vertex]
            elif self.nodes[current_node]["type"] == JOIN_INTRODUCE or self.nodes[current_node][
                "type"] == BIG_JOIN_INTRODUCE:
                v = self.nodes[current_node]["bag"].difference(self.nodes[children[0]]["bag"]).pop()
                for vertex in self.nodes[current_node]["bag"]:
                    if vertex == v:
                        set_of_neighbors = set()
                        for v1 in self.nodes[current_node]["bag"]:
                            if v1 != v and v1[:-1] == v:
                                set_of_neighbors = set_of_neighbors.union(
                                    self.nodes[children[0]]["local_neighbors"][v1])
                        self.nodes[current_node]["local_neighbors"][vertex] = set_of_neighbors
                    else:
                        self.nodes[current_node]["local_neighbors"][vertex] = \
                            self.nodes[children[0]]["local_neighbors"][vertex]
            elif self.nodes[current_node]["type"] == FORGET or self.nodes[current_node]["type"] == ROOT or \
                    self.nodes[current_node]["type"] == JOIN_FORGET:
                for vertex in self.nodes[current_node]["bag"]:
                    self.nodes[current_node]["local_neighbors"][vertex] = \
                        self.nodes[children[0]]["local_neighbors"][vertex]

            else:  # Join node
                for vertex in self.nodes[current_node]["bag"]:
                    if vertex in self.nodes[children[0]]["bag"]:
                        self.nodes[current_node]["local_neighbors"][vertex] = \
                            self.nodes[children[0]]["local_neighbors"][vertex]
                    else:
                        self.nodes[current_node]["local_neighbors"][vertex] = \
                            self.nodes[children[1]]["local_neighbors"][vertex]

    def ntd_to_dntd(self, current_node, debug_flag=False):
        """
        Recursive function that transforms the tree disjoint branch nice form tree decomposition
        (after it is already nice form).
        :param current_node: The current node that we are on TD.
        :param debug_flag: If True, prints the current node and its information.
        :return: None
        """

        bag_of_node = self.nodes[current_node]["bag"]

        if debug_flag:
            print("current node id:" + str(current_node))
            print("current bag:" + str(self.nodes[current_node]["bag"]))
            try:
                print("father:" + str(list(self.predecessors(current_node))[0]))
            except IndexError:
                print("father: None")
            print("children:" + str(list(self.successors(current_node))))
            print("current type:" + self.nodes[current_node]["type"])
            print("-" * 30 + "\n")

        if self.nodes[current_node]["type"] == LEAF:
            return

        children = list(self.successors(current_node))
        if self.nodes[current_node]["type"] == ROOT:
            self.nodes[current_node]["br"] = ""
            return self.ntd_to_dntd(children[0])

        parent_node = next(iter(self.predecessors(current_node)))
        if self.nodes[parent_node]["type"] == JOIN:
            if self.nodes[current_node]["leftCh"]:
                self.nodes[current_node]["br"] = self.nodes[parent_node]["br"] + "0"
            else:
                self.nodes[current_node]["br"] = self.nodes[parent_node]["br"] + "1"
        else:
            self.nodes[current_node]["br"] = self.nodes[parent_node]["br"]

        new_bag = set()
        for vertex in bag_of_node:
            if self.first_appear[vertex] is None:
                self.first_appear[vertex] = current_node

            new_bag += {chr(vertex) + self.nodes[current_node]["br"]}

        self.nodes[current_node]["bag"] = new_bag
        self.new_nodes_dict[current_node] = new_bag

        if debug_flag:
            print("updated current bag:" + str(self.nodes[current_node]["bag"]))
            print("-" * 30 + "\n")
        if self.nodes[current_node]["type"] == JOIN:
            self.nodes[children[0]]["leftCh"] = True
            self.nodes[children[1]]["leftCh"] = False
            self.ntd_to_dntd(children[0])
            self.ntd_to_dntd(children[1])

            new_join_node_bag = self.nodes[children[0]]["bag"] + self.nodes[children[1]]["bag"]
            new_join_node = self.add_node_bag(new_join_node_bag)
            self.add_edge(current_node, new_join_node)
            self.remove_edge(current_node, children[0])
            self.remove_edge(current_node, children[1])
            self.nodes[current_node]["type"] = FORGET
            self.add_edge(new_join_node, children[0])
            self.add_edge(new_join_node, children[1])
            self.nodes[new_join_node]["type"] = JOIN
            self.nodes[new_join_node]["br"] = self.nodes[current_node]["br"]

            current_forget_node = current_node
            for vertex in new_join_node_bag:
                new_forget_node_bag = self.nodes[current_forget_node]["bag"].union({vertex})
                new_forget_node = self.add_node_bag(new_forget_node_bag)
                self.add_edge(current_forget_node, new_forget_node)
                self.remove_edge(current_forget_node, new_join_node)
                self.nodes[current_forget_node]["type"] = JOIN_FORGET
                self.add_edge(new_forget_node, new_join_node)
                self.nodes[new_forget_node]["br"] = self.nodes[current_forget_node]["br"]
                current_forget_node = new_forget_node

            current_introduce_node = current_forget_node
            self.nodes[current_introduce_node]["type"] = BIG_JOIN_INTRODUCE
            for vertex in list(self.nodes[current_node]["bag"])[1:]:
                new_introduce_node_bag = self.nodes[current_introduce_node]["bag"].diffrnce({vertex})
                new_introduce_node = self.add_node_bag(new_introduce_node_bag)
                self.add_edge(current_introduce_node, new_introduce_node)
                self.remove_edge(current_introduce_node, new_join_node)
                self.add_edge(new_introduce_node, new_join_node)
                self.nodes[new_introduce_node]["br"] = self.nodes[current_introduce_node]["br"]
                current_introduce_node = new_introduce_node
                self.nodes[current_introduce_node]["type"] = JOIN_INTRODUCE
        else:
            self.ntd_to_dntd(children[0])

    def semi_ntd_to_semi_dntd(self, current_node, debug_flag=False):
        """
        Recursive function that transforms the tree disjoint branch nice form tree decomposition
        (after it is already nice form).
        :param current_node: The current node that we are on TD.
        :param debug_flag: If True, prints the current node and its information.
        :return: None
        """

        bag_of_node = self.nodes[current_node]["bag"]

        if debug_flag:
            print("current node id:" + str(current_node))
            print("current bag:" + str(self.nodes[current_node]["bag"]))
            try:
                print("father:" + str(list(self.predecessors(current_node))[0]))
            except IndexError:
                print("father: None")
            print("children:" + str(list(self.successors(current_node))))
            print("current type:" + self.nodes[current_node]["type"])
            print("-" * 30 + "\n")

        if self.nodes[current_node]["type"] == LEAF:
            return

        children = list(self.successors(current_node))
        if self.nodes[current_node]["type"] == ROOT:
            self.nodes[current_node]["br"] = ""
            return self.semi_ntd_to_semi_dntd(children[0], debug_flag=debug_flag)

        parent_node = next(iter(self.predecessors(current_node)))
        if self.nodes[parent_node]["type"] == JOIN:
            if self.nodes[current_node]["leftCh"]:
                self.nodes[current_node]["br"] = self.nodes[parent_node]["br"] + "0"
            else:
                self.nodes[current_node]["br"] = self.nodes[parent_node]["br"] + "1"
        else:
            self.nodes[current_node]["br"] = self.nodes[parent_node]["br"]

        new_bag = set()
        for vertex in bag_of_node:
            if self.first_appear[vertex] is None:
                self.first_appear[vertex] = current_node

            new_bag = new_bag.union({chr(vertex) + self.nodes[current_node]["br"]})

        self.nodes[current_node]["bag"] = new_bag
        self.new_nodes_dict[current_node] = new_bag

        if debug_flag:
            print("updated current bag:" + str(self.nodes[current_node]["bag"]))
            print("-" * 30 + "\n")
        if self.nodes[current_node]["type"] == JOIN:
            self.nodes[children[0]]["leftCh"] = True
            self.nodes[children[1]]["leftCh"] = False
            self.semi_ntd_to_semi_dntd(children[0], debug_flag=debug_flag)
            self.semi_ntd_to_semi_dntd(children[1], debug_flag=debug_flag)

            new_join_node_bag = self.nodes[children[0]]["bag"] | self.nodes[children[1]]["bag"]
            new_join_node = self.add_node_bag(new_join_node_bag)
            self.add_edge(current_node, new_join_node)
            self.remove_edge(current_node, children[0])
            self.remove_edge(current_node, children[1])
            self.nodes[current_node]["type"] = FORGET
            self.add_edge(new_join_node, children[0])
            self.add_edge(new_join_node, children[1])
            self.nodes[new_join_node]["type"] = JOIN
            self.nodes[new_join_node]["br"] = self.nodes[current_node]["br"]

            current_forget_node = current_node
            for vertex in sorted(new_join_node_bag):
                new_forget_node_bag = self.nodes[current_forget_node]["bag"].union({vertex})
                new_forget_node = self.add_node_bag(new_forget_node_bag)
                self.add_edge(current_forget_node, new_forget_node)
                self.remove_edge(current_forget_node, new_join_node)
                self.nodes[current_forget_node]["type"] = JOIN_FORGET
                self.add_edge(new_forget_node, new_join_node)
                self.nodes[new_forget_node]["br"] = self.nodes[current_forget_node]["br"]
                current_forget_node = new_forget_node

            current_introduce_node = current_forget_node
            self.nodes[current_introduce_node]["type"] = BIG_JOIN_INTRODUCE
            for vertex in sorted(self.nodes[current_node]["bag"])[1:]:
                new_introduce_node_bag = self.nodes[current_introduce_node]["bag"].difference({vertex})
                new_introduce_node = self.add_node_bag(new_introduce_node_bag)
                self.add_edge(current_introduce_node, new_introduce_node)
                self.remove_edge(current_introduce_node, new_join_node)
                self.add_edge(new_introduce_node, new_join_node)
                self.nodes[new_introduce_node]["br"] = self.nodes[current_introduce_node]["br"]
                current_introduce_node = new_introduce_node
                self.nodes[current_introduce_node]["type"] = JOIN_INTRODUCE
        else:
            self.semi_ntd_to_semi_dntd(children[0], debug_flag=debug_flag)

    def create_factors(self):

        for node in self.nodes:
            self.nodes[node]["order_for_factor"] = tuple()
            if self.nodes[node]["type"] in [LEAF, JOIN, ROOT, BIG_JOIN_INTRODUCE, JOIN_INTRODUCE, INTRODUCE,
                                            JOIN_FORGET, FORGET]:
                self.nodes[node]["factor"] = MemoTable()
            else:
                print("Error: unknown node type")
                exit(-1)

    def question_factor(self, assignment: dict, node_id: int):
        # build the right tuple by the factor order
        order = self.nodes[node_id]["order_for_factor"]
        print("order: " + str(order))
        right_tuple = tuple()
        for vertex in order:
            right_tuple += (assignment[vertex],)
        try:
            return self.nodes[node_id]["factor"].get_value(right_tuple)
        except KeyError:
            return "Key Error"

    def calculate_factors(self, current_node, with_options = False):
        """
        This is a dynamic programming algorithm that calculates the factors of the TD.
        In order to implement the algorithm for enumeration of dominating sets.
        This what we call in our paper: pre-processing phase.

        :param current_node: The current node that we are on TD.
        :return: None
        """

        def calculate_k_range(vertex: str, label_lower: int, label_upper: int, assignment: dict):

            #TODO: why not local neighbors?
            neighbors = {chr(n) for n in self.original_graph.neighbors(ord(vertex[0]))}
            v_label = {w[0] for w, label in assignment.items() if label_lower <= label <= label_upper}

            return len(neighbors.intersection(v_label))

        def calculate_k(vertex: str, label_a: str, assignment: dict):

            if label_a == "S":
                return calculate_k_range(vertex, SI, S1, assignment)
            elif label_a == "W":
                return calculate_k_range(vertex, W0, W1, assignment)
            elif label_a == "R":
                return calculate_k_range(vertex, R0, R2, assignment)
            else:
                return calculate_k_range(vertex, trans_dict[label_a], trans_dict[label_a], assignment)

        type_of_node = self.nodes[current_node]["type"]
        children = list(self.successors(current_node))

        if type_of_node != LEAF:
            for child in children:
                self.calculate_factors(child)
            child_node = children[0]

        if type_of_node == LEAF:
            self.nodes[current_node]["factor"].set_value(frozendict(), 1)
        elif type_of_node == BIG_JOIN_INTRODUCE:

            set_of_real_vertices = {v[0] + self.nodes[current_node]["br"] for v in self.nodes[current_node]["bag"]}
            dict_of_copies = {v:[] for v in set_of_real_vertices}

            for v in set_of_real_vertices:
                if (v + "0") in self.nodes[current_node]["bag"]:
                    dict_of_copies[v].append(v + "0")
                if (v + "1") in self.nodes[current_node]["bag"]:
                    dict_of_copies[v].append(v + "1")

            # the new vertex that was introduced
            v = self.nodes[current_node]["bag"].difference(self.nodes[child_node]["bag"]).pop()

            for key in self.nodes[child_node]["factor"].get_all_keys():

                phi = dict(key)
                phi[v] = "N"
                flag = True

                for original_vertex in set_of_real_vertices:
                    if flag is False:
                        break
                    if len(dict_of_copies[original_vertex]) == 1:
                        first_copy = dict_of_copies[original_vertex][0]
                        label = phi[first_copy]
                        if original_vertex in phi.keys(): # TODO: why isn't it always in phi.keys()?
                            phi[original_vertex] = label
                        continue

                    first_copy = dict_of_copies[original_vertex][0]
                    second_copy = dict_of_copies[original_vertex][1]

                    label_0 = phi[first_copy]
                    label_1 = phi[second_copy]

                    flag, label_original_vertex = join_labels(label_0, label_1)
                    if not flag:
                        continue
                    phi[original_vertex] = label_original_vertex
                if with_options and not (phi[v] in self.original_graph.nodes[ord(v[0])]["options"]):
                    flag = False
                if flag:
                    self.nodes[current_node]["factor"].set_value(frozendict(phi), 1)

        elif type_of_node == JOIN_INTRODUCE:

            # the new vertex that was introduced
            v = self.nodes[current_node]["bag"].difference(self.nodes[child_node]["bag"]).pop()

            for key in self.nodes[child_node]["factor"].get_all_keys():
                old_value = self.nodes[child_node]["factor"].get_value(key)
                if old_value == 0:
                    continue
                for label in F:
                    if with_options and not (label in self.original_graph.nodes[ord(v[0])]["options"]):
                        continue
                    phi = dict(key)
                    phi[v] = label
                    self.nodes[current_node]["factor"].set_value(frozendict(phi), 1)

        elif type_of_node == INTRODUCE:

            # the new vertex that was introduced
            v = self.nodes[current_node]["bag"].difference(self.nodes[child_node]["bag"]).pop()

            for key in self.nodes[child_node]["factor"].get_all_keys():
                if self.nodes[child_node]["type"] != LEAF and key == frozendict(): # TODO: why is it even possible?
                    print(f"Error: empty key at node {self.nodes[child_node]['bag']}")
                    exit(-1)
                old_value = self.nodes[child_node]["factor"].get_value(key)
                if old_value == 0:
                    continue
                phi = dict(key)
                for label in [SI, R0, W0, S0]:
                    if with_options and not (label in self.original_graph.nodes[ord(v[0])]["options"]):
                        continue
                    phi[v] = label
                    if label == SI:
                        k_v_s = calculate_k(v, "S", key)
                        if k_v_s == 0:
                            self.nodes[current_node]["factor"].set_value(frozendict(phi), 1)
                        else:
                            self.nodes[current_node]["factor"].set_value(frozendict(phi), 0)
                    elif label == S0:
                        k_v_si = calculate_k(v, "SI", key)
                        if k_v_si == 0:
                            self.nodes[current_node]["factor"].set_value(frozendict(phi), 1)
                        else:
                            self.nodes[current_node]["factor"].set_value(frozendict(phi), 0)
                    else:
                        self.nodes[current_node]["factor"].set_value(frozendict(phi), 1)

        elif type_of_node == JOIN_FORGET:

            v = self.nodes[child_node]["bag"].difference(self.nodes[current_node]["bag"]).pop()

            for key in self.nodes[child_node]["factor"].get_all_keys():
                # Just copy the factor
                old_value = self.nodes[child_node]["factor"].get_value(key)
                if old_value == 0:
                    continue
                new_key = dict(key)
                del new_key[v]
                new_key = frozendict(new_key)
                self.nodes[current_node]["factor"].set_value(new_key, 1)

        elif type_of_node == FORGET or type_of_node == ROOT:

            v = self.nodes[child_node]["bag"].difference(self.nodes[current_node]["bag"]).pop()

            # All possible keys pass
            for key in self.nodes[child_node]["factor"].get_all_keys():
                old_value = self.nodes[child_node]["factor"].get_value(key)
                if old_value == 0:
                    continue

                v_label = key[v]
                k_v_si = calculate_k(v, "SI", key)
                k_v_s = calculate_k(v, "S", key)
                k_v_w = calculate_k(v, "W", key)
                k_v_w1 = calculate_k(v, "W1", key)
                k_v_w0 = calculate_k(v, "W0", key)

                if v_label == SI:
                    if not (k_v_s == 0 and k_v_w == 0):
                        self.nodes[child_node]["factor"].set_value(key, 0)
                        continue
                elif in_rho(v_label):
                    j = v_label - R0
                    if not (j + k_v_s >= 2):
                        self.nodes[child_node]["factor"].set_value(key, 0)
                        continue
                elif v_label == W1:
                    if not (k_v_s == 0):
                        self.nodes[child_node]["factor"].set_value(key, 0)
                        continue
                elif v_label == W0:
                    if not (k_v_si == 0 and k_v_s == 1):
                        self.nodes[child_node]["factor"].set_value(key, 0)
                        continue
                elif v_label == S1:
                    if not (k_v_si == 0 and k_v_w1 == 0):
                        self.nodes[child_node]["factor"].set_value(key, 0)
                        continue
                elif v_label == S0:
                    if not (k_v_si == 0 and k_v_w1 == 0 and k_v_w0 >= 1):
                        self.nodes[child_node]["factor"].set_value(key, 0)
                        continue

                options_for_new_labels = {z: set() for z in self.nodes[current_node]["bag"]}

                flag = True
                for w in key.keys():
                    if w==v:
                        continue
                    if v_label == SI:
                        if not (ord(w[0]) in self.original_graph.neighbors(ord(v[0]))):
                            options_for_new_labels[w].add(key[w])
                        elif in_rho(key[w]):
                            for j in range(3):
                                if max(0, j-1) == key[w]-R0:
                                    options_for_new_labels[w].add(trans_dict["R" + str(j)])
                                    # TODO: can't there be break here?
                        else:
                            flag = False
                            break
                    elif v_label == S0:
                        if not (ord(w[0]) in self.original_graph.neighbors(ord(v[0]))) or (S0 <= key[w] <= S1):
                            options_for_new_labels[w].add(key[w])
                        elif in_rho(key[w]):
                            for j in range(3):
                                if max(0, j - 1) == key[w] - R0:
                                    options_for_new_labels[w].add(trans_dict["R" + str(j)])
                        elif key[w] == W0:
                            options_for_new_labels[w].add(W1)
                        else:
                            flag = False
                            break
                    elif v_label == S1:
                        if not (ord(w[0]) in self.original_graph.neighbors(ord(v[0]))) or (S0 <= key[w] <= S1):
                            options_for_new_labels[w].add(key[w])
                        elif in_rho(key[w]):
                            for j in range(3):
                                if max(0, j - 1) == key[w] - R0:
                                    options_for_new_labels[w].add(trans_dict["R" + str(j)])
                        elif key[w] == W0:
                            options_for_new_labels[w].add(W1)
                        else:
                            flag = False
                            break
                    elif v_label == W1:
                        if not (ord(w[0]) in self.original_graph.neighbors(ord(v[0]))) or not (in_sigma(key[w])):
                            options_for_new_labels[w].add(key[w])
                        else:
                            flag = False
                            break
                    elif v_label == W0:
                        if (not (ord(w[0]) in self.original_graph.neighbors(ord(v[0]))) or not (in_sigma(key[w]))):
                            options_for_new_labels[w].add(key[w])
                        elif S0 <= key[w] <= S1:
                            options_for_new_labels[w].add(S1)
                        else:
                            flag = False
                            break
                    else:
                        options_for_new_labels[w].add(key[w])
                if flag:
                    for z in options_for_new_labels.keys():
                        if len(options_for_new_labels[z]) == 0:
                            flag = False
                            break
                        else:
                            options_for_new_labels[z] = list(options_for_new_labels[z])
                if flag:
                    for opt in generate_dictionaries_from_sets(options_for_new_labels):
                        self.nodes[current_node]["factor"].set_value(frozendict(opt), 1)


        else:
            child_node_1 = children[0]
            child_node_2 = children[1]

            # taking a simple conjunction of the two factors
            # Here is a kind of natural join algorithm on the two factors, where we check if eah node has the same label
            # in both factors, and if so on all the nodes in that key, we add the combine key it to the new factor.


            for key_1 in self.nodes[child_node_1]["factor"].get_all_keys():
                if self.nodes[child_node_1]["factor"].get_value(key_1) == 0:
                    continue
                for key_2 in self.nodes[child_node_2]["factor"].get_all_keys():
                    if self.nodes[child_node_2]["factor"].get_value(key_2) == 0:
                        continue
                    flag = True
                    new_key = {z: 0 for z in self.nodes[current_node]["bag"]}
                    for w in key_1.keys():
                        label1 = key_1[w]
                        w2 = w[:-1] + str((int(w[-1]) + 1) % 2)
                        if w2 in key_2.keys():
                            label2 = key_2[w2]
                            if ((same_label_group(label1, label2)) or (in_omega(label1) and in_rho(label2)) or
                                    (in_rho(label1) and in_omega(label2))):
                                new_key[w] = label1
                                new_key[w2] = label2
                            else:
                                flag = False
                                break
                        else:
                            new_key[w] = label1
                    for w in key_2.keys():
                        if new_key[w] == 0:
                            new_key[w] = key_2[w]
                    if flag:
                        self.nodes[current_node]["factor"].set_value(frozendict(new_key), 1)

        if type_of_node == ROOT:

            return

    def calculate_factors_iterative(self, with_options = False):
        """
        This is a for loop version of calculate_factors, using a stack.
        """
        def calculate_k_range(vertex: str, label_lower: int, label_upper: int, assignment: dict):

            #TODO: why not local neighbors?
            neighbors = {chr(n) for n in self.original_graph.neighbors(ord(vertex[0]))}
            v_label = {w[0] for w, label in assignment.items() if label_lower <= label <= label_upper}

            return len(neighbors.intersection(v_label))

        def calculate_k(vertex: str, label_a: str, assignment: dict):

            if label_a == "S":
                return calculate_k_range(vertex, SI, S1, assignment)
            elif label_a == "W":
                return calculate_k_range(vertex, W0, W1, assignment)
            elif label_a == "R":
                return calculate_k_range(vertex, R0, R2, assignment)
            else:
                return calculate_k_range(vertex, trans_dict[label_a], trans_dict[label_a], assignment)

        stack = [self.get_root()]
        while len(stack) != 0:

            current_node = stack[-1]

            type_of_node = self.nodes[current_node]["type"]
            children = list(self.successors(current_node))

            # Check if children have been processed
            children_processed = True
            if type_of_node != LEAF:
                for child in children:
                    if "processed" not in self.nodes[child] or not self.nodes[child]["processed"]:
                        children_processed = False
                        stack.append(child)
                child_node = children[0]

            if not children_processed:
                continue  # Continue processing children

            # All children processed (or leaf node), process current node
            stack.pop()  # Remove the current node from the stack

            if type_of_node == LEAF:
                self.nodes[current_node]["factor"].set_value(frozendict(), 1)
            elif type_of_node == BIG_JOIN_INTRODUCE:

                set_of_real_vertices = {v[0] + self.nodes[current_node]["br"] for v in self.nodes[current_node]["bag"]}
                dict_of_copies = {v: [] for v in set_of_real_vertices}

                for v in set_of_real_vertices:
                    if (v + "0") in self.nodes[current_node]["bag"]:
                        dict_of_copies[v].append(v + "0")
                    if (v + "1") in self.nodes[current_node]["bag"]:
                        dict_of_copies[v].append(v + "1")

                # the new vertex that was introduced
                v = self.nodes[current_node]["bag"].difference(self.nodes[child_node]["bag"]).pop()

                for key in self.nodes[child_node]["factor"].get_all_keys():

                    phi = dict(key)
                    phi[v] = "N"
                    flag = True

                    for original_vertex in set_of_real_vertices:
                        if flag is False:
                            break
                        if len(dict_of_copies[original_vertex]) == 1:
                            first_copy = dict_of_copies[original_vertex][0]
                            label = phi[first_copy]
                            if original_vertex in phi.keys():  # TODO: why isn't it always in phi.keys()?
                                phi[original_vertex] = label
                            continue

                        first_copy = dict_of_copies[original_vertex][0]
                        second_copy = dict_of_copies[original_vertex][1]

                        label_0 = phi[first_copy]
                        label_1 = phi[second_copy]

                        flag, label_original_vertex = join_labels(label_0, label_1)
                        if not flag:
                            continue
                        phi[original_vertex] = label_original_vertex
                    if with_options and not (phi[v] in self.original_graph.nodes[ord(v[0])]["options"]):
                        flag = False
                    if flag:
                        self.nodes[current_node]["factor"].set_value(frozendict(phi), 1)

            elif type_of_node == JOIN_INTRODUCE:

                # the new vertex that was introduced
                v = self.nodes[current_node]["bag"].difference(self.nodes[child_node]["bag"]).pop()

                for key in self.nodes[child_node]["factor"].get_all_keys():
                    old_value = self.nodes[child_node]["factor"].get_value(key)
                    if old_value == 0:
                        continue
                    for label in F:
                        if with_options and not (label in self.original_graph.nodes[ord(v[0])]["options"]):
                            continue
                        phi = dict(key)
                        phi[v] = label
                        self.nodes[current_node]["factor"].set_value(frozendict(phi), 1)

            elif type_of_node == INTRODUCE:

                # the new vertex that was introduced
                v = self.nodes[current_node]["bag"].difference(self.nodes[child_node]["bag"]).pop()

                for key in self.nodes[child_node]["factor"].get_all_keys():
                    if self.nodes[child_node]["type"] != LEAF and key == frozendict():  # TODO: why is it even possible?
                        print(f"Error: empty key at node {self.nodes[child_node]['bag']}")
                        exit(-1)
                    old_value = self.nodes[child_node]["factor"].get_value(key)
                    if old_value == 0:
                        continue
                    phi = dict(key)
                    for label in [SI, R0, W0, S0]:
                        if with_options and not (label in self.original_graph.nodes[ord(v[0])]["options"]):
                            continue
                        phi[v] = label
                        if label == SI:
                            k_v_s = calculate_k(v, "S", key)
                            if k_v_s == 0:
                                self.nodes[current_node]["factor"].set_value(frozendict(phi), 1)
                            else:
                                self.nodes[current_node]["factor"].set_value(frozendict(phi), 0)
                        elif label == S0:
                            k_v_si = calculate_k(v, "SI", key)
                            if k_v_si == 0:
                                self.nodes[current_node]["factor"].set_value(frozendict(phi), 1)
                            else:
                                self.nodes[current_node]["factor"].set_value(frozendict(phi), 0)
                        else:
                            self.nodes[current_node]["factor"].set_value(frozendict(phi), 1)

            elif type_of_node == JOIN_FORGET:

                v = self.nodes[child_node]["bag"].difference(self.nodes[current_node]["bag"]).pop()

                for key in self.nodes[child_node]["factor"].get_all_keys():
                    # Just copy the factor
                    old_value = self.nodes[child_node]["factor"].get_value(key)
                    if old_value == 0:
                        continue
                    new_key = dict(key)
                    del new_key[v]
                    new_key = frozendict(new_key)
                    self.nodes[current_node]["factor"].set_value(new_key, 1)

            elif type_of_node == FORGET or type_of_node == ROOT:

                v = self.nodes[child_node]["bag"].difference(self.nodes[current_node]["bag"]).pop()

                # All possible keys pass
                for key in self.nodes[child_node]["factor"].get_all_keys():
                    old_value = self.nodes[child_node]["factor"].get_value(key)
                    if old_value == 0:
                        continue

                    v_label = key[v]
                    k_v_si = calculate_k(v, "SI", key)
                    k_v_s = calculate_k(v, "S", key)
                    k_v_w = calculate_k(v, "W", key)
                    k_v_w1 = calculate_k(v, "W1", key)
                    k_v_w0 = calculate_k(v, "W0", key)

                    if v_label == SI:
                        if not (k_v_s == 0 and k_v_w == 0):
                            self.nodes[child_node]["factor"].set_value(key, 0)
                            continue
                    elif in_rho(v_label):
                        j = v_label - R0
                        if not (j + k_v_s >= 2):
                            self.nodes[child_node]["factor"].set_value(key, 0)
                            continue
                    elif v_label == W1:
                        if not (k_v_s == 0):
                            self.nodes[child_node]["factor"].set_value(key, 0)
                            continue
                    elif v_label == W0:
                        if not (k_v_si == 0 and k_v_s == 1):
                            self.nodes[child_node]["factor"].set_value(key, 0)
                            continue
                    elif v_label == S1:
                        if not (k_v_si == 0 and k_v_w1 == 0):
                            self.nodes[child_node]["factor"].set_value(key, 0)
                            continue
                    elif v_label == S0:
                        if not (k_v_si == 0 and k_v_w1 == 0 and k_v_w0 >= 1):
                            self.nodes[child_node]["factor"].set_value(key, 0)
                            continue

                    options_for_new_labels = {z: set() for z in self.nodes[current_node]["bag"]}

                    flag = True
                    for w in key.keys():
                        if w == v:
                            continue
                        if v_label == SI:
                            if not (ord(w[0]) in self.original_graph.neighbors(ord(v[0]))):
                                options_for_new_labels[w].add(key[w])
                            elif in_rho(key[w]):
                                for j in range(3):
                                    if max(0, j - 1) == key[w] - R0:
                                        options_for_new_labels[w].add(trans_dict["R" + str(j)])
                                        # TODO: can't there be break here?
                            else:
                                flag = False
                                break
                        elif v_label == S0:
                            if not (ord(w[0]) in self.original_graph.neighbors(ord(v[0]))) or (S0 <= key[w] <= S1):
                                options_for_new_labels[w].add(key[w])
                            elif in_rho(key[w]):
                                for j in range(3):
                                    if max(0, j - 1) == key[w] - R0:
                                        options_for_new_labels[w].add(trans_dict["R" + str(j)])
                            elif key[w] == W0:
                                options_for_new_labels[w].add(W1)
                            else:
                                flag = False
                                break
                        elif v_label == S1:
                            if not (ord(w[0]) in self.original_graph.neighbors(ord(v[0]))) or (S0 <= key[w] <= S1):
                                options_for_new_labels[w].add(key[w])
                            elif in_rho(key[w]):
                                for j in range(3):
                                    if max(0, j - 1) == key[w] - R0:
                                        options_for_new_labels[w].add(trans_dict["R" + str(j)])
                            elif key[w] == W0:
                                options_for_new_labels[w].add(W1)
                            else:
                                flag = False
                                break
                        elif v_label == W1:
                            if not (ord(w[0]) in self.original_graph.neighbors(ord(v[0]))) or not (in_sigma(key[w])):
                                options_for_new_labels[w].add(key[w])
                            else:
                                flag = False
                                break
                        elif v_label == W0:
                            if (not (ord(w[0]) in self.original_graph.neighbors(ord(v[0]))) or not (in_sigma(key[w]))):
                                options_for_new_labels[w].add(key[w])
                            elif S0 <= key[w] <= S1:
                                options_for_new_labels[w].add(S1)
                            else:
                                flag = False
                                break
                        else:
                            options_for_new_labels[w].add(key[w])
                    if flag:
                        for z in options_for_new_labels.keys():
                            if len(options_for_new_labels[z]) == 0:
                                flag = False
                                break
                            else:
                                options_for_new_labels[z] = list(options_for_new_labels[z])
                    if flag:
                        for opt in generate_dictionaries_from_sets(options_for_new_labels):
                            self.nodes[current_node]["factor"].set_value(frozendict(opt), 1)


            elif type_of_node == JOIN:
                child_node_1 = children[0]
                child_node_2 = children[1]

                # taking a simple conjunction of the two factors
                # Here is a kind of natural join algorithm on the two factors, where we check if eah node has the same label
                # in both factors, and if so on all the nodes in that key, we add the combine key it to the new factor.

                for key_1 in self.nodes[child_node_1]["factor"].get_all_keys():
                    if self.nodes[child_node_1]["factor"].get_value(key_1) == 0:
                        continue
                    for key_2 in self.nodes[child_node_2]["factor"].get_all_keys():
                        if self.nodes[child_node_2]["factor"].get_value(key_2) == 0:
                            continue
                        flag = True
                        new_key = {z: 0 for z in self.nodes[current_node]["bag"]}
                        for w in key_1.keys():
                            label1 = key_1[w]
                            w2 = w[:-1] + str((int(w[-1]) + 1) % 2)
                            if w2 in key_2.keys():
                                label2 = key_2[w2]
                                if ((same_label_group(label1, label2)) or (in_omega(label1) and in_rho(label2)) or
                                        (in_rho(label1) and in_omega(label2))):
                                    new_key[w] = label1
                                    new_key[w2] = label2
                                else:
                                    flag = False
                                    break
                            else:
                                new_key[w] = label1
                        for w in key_2.keys():
                            if new_key[w] == 0:
                                new_key[w] = key_2[w]
                        if flag:
                            self.nodes[current_node]["factor"].set_value(frozendict(new_key), 1)

            # Mark the current node as processed
            self.nodes[current_node]["processed"] = True


    def EnumDS(self, theta: Dict[str, int], i=0, debug_flag = False) -> None:
        """
        This algorithm means to enumerate all the dominating sets of the graph.
        :param theta: An extendable labeling.
        :param i: The index of the vertex in the graph (in Q).
        :return:
        """
        if i == len(self.all_vertices):
            yield frozenset({self.original_graph.nodes[ord(x[0])]["original_name"] for x in V_label("S", theta)})
            return
        V_label_S, V_label_W = V_label_S_W(theta)
        for c in F:
            if debug_flag:
                print("Current theta: " + str(theta))
                print("Current vertex: " + str(self.Q[i]))
                print("Current node: " + str(self.nodes[self.first_appear[self.Q[i]]]["bag"]))
                print("Current br: " + str(self.nodes[self.first_appear[self.Q[i]]]["br"]))
                print("Optional label: " + str(c))
            counter = 0
            for v in self.nodes[self.first_appear[self.Q[i]]]["bag"]:
                if v[0] == self.Q[i][0]:
                    counter += 1
            if counter == 1:
                new_theta = self.IncrementLabeling(theta, i, c, V_label_S, V_label_W)
            elif counter == 2:
                new_theta = self.IncrementLabeling2(theta, i, c)
            elif counter == 3:
                original_copy = self.Q[i][0] + self.nodes[self.first_appear[self.Q[i]]]["br"]
                original_c = theta[original_copy]
                first_copy = self.Q[i][0] + self.nodes[self.first_appear[self.Q[i]]]["br"] + "0"
                first_c = theta[first_copy]
                if in_rho(original_c):
                    if in_rho(c) and in_rho(first_c) and original_c - R0 == c - R0 + first_c - R0:
                        if first_c == R1 and c == R1:
                            new_theta = self.IncrementLabeling2(theta, i, c)
                        elif first_c == R0 and c == R0:
                            new_theta = self.IncrementLabeling2(theta, i, c)
                        elif first_c == R0 and c == R1:
                            new_theta = self.IncrementLabeling2(theta, i, c)
                        else:
                            if debug_flag:
                                print("Not Valid Labeling")
                                print("-" * 20)
                            continue
                    elif original_c == R1 and first_c == R1 and c == W0:
                        new_theta = self.IncrementLabeling2(theta, i, c)
                    elif original_c == R2 and first_c == W0 and c == R2:
                        new_theta = self.IncrementLabeling2(theta, i, c)
                    elif original_c == R2 and first_c == R2 and c == W0:
                        new_theta = self.IncrementLabeling2(theta, i, c)
                    else:
                        if debug_flag:
                            print("Not Valid Labeling")
                            print("-" * 20)
                        continue
                elif in_sigma(original_c) and in_sigma(first_c) and in_sigma(c):
                    if original_c == SI and first_c == SI and c == SI:
                        new_theta = self.IncrementLabeling2(theta, i, c)
                    elif original_c == S0 and first_c == S0 and c == S0:
                        new_theta = self.IncrementLabeling2(theta, i, c)
                    elif original_c == S1 and first_c == S1 and c == S0:
                        new_theta = self.IncrementLabeling2(theta, i, c)
                    elif original_c == S1 and first_c == S0 and c == S1:
                        new_theta = self.IncrementLabeling2(theta, i, c)
                    elif original_c == S1 and first_c == S1 and c == S1:
                        new_theta = self.IncrementLabeling2(theta, i, c)
                    else:
                        if debug_flag:
                            print("Not Valid Labeling")
                            print("-" * 20)
                        continue
                elif in_omega(original_c) and in_omega(first_c) and in_omega(c):
                    if original_c == W0 and first_c == W0 and c == W0:
                        new_theta = self.IncrementLabeling2(theta, i, c)
                    elif original_c == W1 and first_c == W0 and c == W1:
                        new_theta = self.IncrementLabeling2(theta, i, c)
                    elif original_c == W1 and first_c == W1 and c == W0:
                        new_theta = self.IncrementLabeling2(theta, i, c)
                    else:
                        if debug_flag:
                            print("Not Valid Labeling")
                            print("-" * 20)
                        continue
                else:
                    if debug_flag:
                        print("Not Valid Labeling")
                        print("-" * 20)
                    continue
            else:
                print("Error - First Appear isn't good")
                return -1
            if debug_flag:
                print("IncrementLabeling: " + str(new_theta))
                print("-" * 20)
            if new_theta is None or not new_theta:
                continue
            for option in new_theta:
                if debug_flag:
                    print("Option: " + str(option))
                    print("IsExtendable: " + str(self.IsExtendable(option, i)))
                    print("-" * 20)
                if self.IsExtendable(option, i):
                    yield from self.EnumDS(option, i + 1, debug_flag=debug_flag)


    def EnumHS(self, theta: Dict[str, int], i=0, debug_flag = False) -> None:
        """
        This algorithm means to enumerate all the minimal hitting sets of a hypergraph (gets it reduction).
        :param theta: An extendable labeling.
        :param i: The index of the vertex in the graph (in Q).
        :return:
        """
        if i == len(self.all_vertices):
            yield frozenset({self.original_graph.nodes[ord(x[0])]["original_name"] for x in V_label("S", theta)})
            return
        options_for_label = self.original_graph.nodes[ord(self.Q[i][0])]["options"]
        V_label_S, V_label_W = V_label_S_W(theta)
        for c in options_for_label:
            if debug_flag:
                print("Current theta: " + str(theta))
                print("Current vertex: " + str(ord(self.Q[i][0])))
                print("Current node: " + str(self.nodes[self.first_appear[self.Q[i]]]["bag"]))
                print("Current br: " + str(self.nodes[self.first_appear[self.Q[i]]]["br"]))
                print("Optional label: " + str(c))
            counter = 0
            for v in self.nodes[self.first_appear[self.Q[i]]]["bag"]:
                if v[0] == self.Q[i][0]:
                    counter += 1
            if counter == 1:
                new_theta = self.IncrementLabeling(theta, i, c, V_label_S, V_label_W)
            elif counter == 2:
                new_theta = self.IncrementLabeling2(theta, i, c)
            elif counter == 3:
                original_copy = self.Q[i][0] + self.nodes[self.first_appear[self.Q[i]]]["br"]
                original_c = theta[original_copy]
                first_copy = self.Q[i][0] + self.nodes[self.first_appear[self.Q[i]]]["br"] + "0"
                first_c = theta[first_copy]
                if in_rho(original_c):
                    if in_rho(c) and in_rho(first_c) and original_c - R0 == c - R0 + first_c - R0:
                        if first_c == R1 and c == R1:
                            new_theta = self.IncrementLabeling2(theta, i, c)
                        elif first_c == R0 and c == R0:
                            new_theta = self.IncrementLabeling2(theta, i, c)
                        elif first_c == R0 and c == R1:
                            new_theta = self.IncrementLabeling2(theta, i, c)
                        else:
                            if debug_flag:
                                print("Not Valid Labeling")
                                print("-" * 20)
                            continue
                    elif original_c == R1 and first_c == R1 and c == W0:
                        new_theta = self.IncrementLabeling2(theta, i, c)
                    elif original_c == R2 and first_c == W0 and c == R2:
                        new_theta = self.IncrementLabeling2(theta, i, c)
                    elif original_c == R2 and first_c == R2 and c == W0:
                        new_theta = self.IncrementLabeling2(theta, i, c)
                    else:
                        if debug_flag:
                            print("Not Valid Labeling")
                            print("-" * 20)
                        continue
                elif in_sigma(original_c) and in_sigma(first_c) and in_sigma(c):
                    if original_c == SI and first_c == SI and c == SI:
                        new_theta = self.IncrementLabeling2(theta, i, c)
                    elif original_c == S0 and first_c == S0 and c == S0:
                        new_theta = self.IncrementLabeling2(theta, i, c)
                    elif original_c == S1 and first_c == S1 and c == S0:
                        new_theta = self.IncrementLabeling2(theta, i, c)
                    elif original_c == S1 and first_c == S0 and c == S1:
                        new_theta = self.IncrementLabeling2(theta, i, c)
                    elif original_c == S1 and first_c == S1 and c == S1:
                        new_theta = self.IncrementLabeling2(theta, i, c)
                    else:
                        if debug_flag:
                            print("Not Valid Labeling")
                            print("-" * 20)
                        continue
                elif in_omega(original_c) and in_omega(first_c) and in_omega(c):
                    if original_c == W0 and first_c == W0 and c == W0:
                        new_theta = self.IncrementLabeling2(theta, i, c)
                    elif original_c == W1 and first_c == W0 and c == W1:
                        new_theta = self.IncrementLabeling2(theta, i, c)
                    elif original_c == W1 and first_c == W1 and c == W0:
                        new_theta = self.IncrementLabeling2(theta, i, c)
                    else:
                        if debug_flag:
                            print("Not Valid Labeling")
                            print("-" * 20)
                        continue
                else:
                    if debug_flag:
                        print("Not Valid Labeling")
                        print("-" * 20)
                    continue
            else:
                print("Error - First Appear isn't good")
                return -1
            if debug_flag:
                print("IncrementLabeling: " + str(new_theta))
                print("-" * 20)
            if new_theta is None or not new_theta:
                continue
            for option in new_theta:
                if debug_flag:
                    print("Option: " + str(option))
                    print("IsExtendable: " + str(self.IsExtendable(option, i)))
                    print("-" * 20)
                if self.IsExtendable(option, i):
                    yield from self.EnumHS(option, i + 1, debug_flag=debug_flag)

    def EnumHS_iterative(self, debug_flag = False) -> None:
        """
        This is a for loop version of calculate_factors, using a stack.
        """
        stack = [(dict(), 0)]

        while stack:

            theta, i = stack.pop()

            if i == len(self.all_vertices):
                yield frozenset({self.original_graph.nodes[ord(x[0])]["original_name"] for x in V_label("S", theta)})
                continue

            options_for_label = self.original_graph.nodes[ord(self.Q[i][0])]["options"]
            V_label_S, V_label_W = V_label_S_W(theta)

            for c in options_for_label:
                if debug_flag:
                    print("Current theta: " + str(theta))
                    print("Current vertex: " + str(ord(self.Q[i][0])))
                    print("Current node: " + str(self.nodes[self.first_appear[self.Q[i]]]["bag"]))
                    print("Current br: " + str(self.nodes[self.first_appear[self.Q[i]]]["br"]))
                    print("Optional label: " + str(c))
                counter = 0
                for v in self.nodes[self.first_appear[self.Q[i]]]["bag"]:
                    if v[0] == self.Q[i][0]:
                        counter += 1
                if counter == 1:
                    new_theta = self.IncrementLabeling(theta, i, c, V_label_S, V_label_W)
                elif counter == 2:
                    new_theta = self.IncrementLabeling2(theta, i, c)
                elif counter == 3:
                    original_copy = self.Q[i][0] + self.nodes[self.first_appear[self.Q[i]]]["br"]
                    original_c = theta[original_copy]
                    first_copy = self.Q[i][0] + self.nodes[self.first_appear[self.Q[i]]]["br"] + "0"
                    first_c = theta[first_copy]
                    if in_rho(original_c):
                        if in_rho(c) and in_rho(first_c) and original_c - R0 == c - R0 + first_c - R0:
                            if first_c == R1 and c == R1:
                                new_theta = self.IncrementLabeling2(theta, i, c)
                            elif first_c == R0 and c == R0:
                                new_theta = self.IncrementLabeling2(theta, i, c)
                            elif first_c == R0 and c == R1:
                                new_theta = self.IncrementLabeling2(theta, i, c)
                            else:
                                if debug_flag:
                                    print("Not Valid Labeling")
                                    print("-" * 20)
                                continue
                        elif original_c == R1 and first_c == R1 and c == W0:
                            new_theta = self.IncrementLabeling2(theta, i, c)
                        elif original_c == R2 and first_c == W0 and c == R2:
                            new_theta = self.IncrementLabeling2(theta, i, c)
                        elif original_c == R2 and first_c == R2 and c == W0:
                            new_theta = self.IncrementLabeling2(theta, i, c)
                        else:
                            if debug_flag:
                                print("Not Valid Labeling")
                                print("-" * 20)
                            continue
                    elif in_sigma(original_c) and in_sigma(first_c) and in_sigma(c):
                        if original_c == SI and first_c == SI and c == SI:
                            new_theta = self.IncrementLabeling2(theta, i, c)
                        elif original_c == S0 and first_c == S0 and c == S0:
                            new_theta = self.IncrementLabeling2(theta, i, c)
                        elif original_c == S1 and first_c == S1 and c == S0:
                            new_theta = self.IncrementLabeling2(theta, i, c)
                        elif original_c == S1 and first_c == S0 and c == S1:
                            new_theta = self.IncrementLabeling2(theta, i, c)
                        elif original_c == S1 and first_c == S1 and c == S1:
                            new_theta = self.IncrementLabeling2(theta, i, c)
                        else:
                            if debug_flag:
                                print("Not Valid Labeling")
                                print("-" * 20)
                            continue
                    elif in_omega(original_c) and in_omega(first_c) and in_omega(c):
                        if original_c == W0 and first_c == W0 and c == W0:
                            new_theta = self.IncrementLabeling2(theta, i, c)
                        elif original_c == W1 and first_c == W0 and c == W1:
                            new_theta = self.IncrementLabeling2(theta, i, c)
                        elif original_c == W1 and first_c == W1 and c == W0:
                            new_theta = self.IncrementLabeling2(theta, i, c)
                        else:
                            if debug_flag:
                                print("Not Valid Labeling")
                                print("-" * 20)
                            continue
                    else:
                        if debug_flag:
                            print("Not Valid Labeling")
                            print("-" * 20)
                        continue
                else:
                    print("Error - First Appear isn't good")
                    return -1
                if debug_flag:
                    print("IncrementLabeling: " + str(new_theta))
                    print("-" * 20)
                if new_theta is None or not new_theta:
                    continue
                for option in new_theta:
                    if debug_flag:
                        print("Option: " + str(option))
                        print("IsExtendable: " + str(self.IsExtendable(option, i)))
                        print("-" * 20)
                    if self.IsExtendable(option, i):
                        stack.append((option, i + 1))

    def IncrementLabeling2(self, theta: Dict[str, int], i, c: int):
        """
        Procedure IncrementLabeling receives as input a labeling which we assume to be extendable (see EnumDS),
        and a label. It generates a new assignment and updates the labels of vertices based on the given label, so that
        the new assignment is legal. [Taken from paper]
        :param theta: Previous labeling.
        :param i: The index of the vertex in the graph (in Q).
        :param c: The label to be added to the vertex.
        :return:
        """
        new_theta = dict(theta)
        new_theta[self.Q[i]] = c
        if i == 0:
            return [new_theta]
        return [new_theta]


    def IncrementLabeling(self, theta: Dict[str, int], i, c: int, V_label_S, V_label_W):
        """
        Procedure IncrementLabeling receives as input a labeling which we assume to be extendable (see EnumDS),
        and a label. It generates a new assignment and updates the labels of vertices based on the given label, so that
        the new assignment is legal. [Taken from paper]
        :param theta: Previous labeling.
        :param i: The index of the vertex in the graph (in Q).
        :param c: The label to be added to the vertex.
        :return:
        """
        new_theta = dict(theta)
        new_theta[self.Q[i]] = c
        if i == 0:
            return [new_theta]
        current_vertex = self.Q[i]

        K_i = self.nodes[self.first_appear[current_vertex]]["local_neighbors"][current_vertex].intersection(
            {w[0] for w in self.Q[:i]})
        K_i = {x+self.nodes[self.first_appear[current_vertex]]["br"] for x in K_i}
        N_i = K_i.intersection(V_label_S)
        W_i = K_i.intersection(V_label_W)

        flag_of_two = False

        if in_sigma(c):
            for v in K_i:
                if in_rho(theta[v]):
                    new_theta[v] = max(R0, theta[v] - 1)

        if c == SI and (len(N_i) != 0 or len(W_i) != 0):
            return False
        if S0 <= c <= S1:
            if len([w for w in K_i if theta[w] in {SI, W0}]) != 0 or \
                    (c == S0 and len(W_i) == 0):
                return False
            else:
                for w in W_i:
                    if theta[w] == W1:
                        new_theta[w] = W0
        if in_omega(c):
            if len([w for w in N_i if theta[w] == SI]) != 0 or \
                    len(N_i) >= 2 or \
                    (len(N_i) == 0 and c == W0) or \
                    (len(N_i) != 0 and c == W1):
                return False
            elif c == W0:
                v = N_i.pop()
                if theta[v] == S0:
                    return False
                flag_of_two = v
        if in_rho(c) and max(0, 2 - len(N_i)) != c - R0:
             return False
        if flag_of_two:
            new_theta[flag_of_two] = S0
            new_theta2 = dict(new_theta)
            new_theta2[flag_of_two] = S1
            return [new_theta, new_theta2]
        return [new_theta]

    def IsExtendable(self, theta, i):
        """
        This function uses the pre-processing phase of the TD to check if the labeling is extendable.
        :param theta: A labeling or False.
        :param i: The index of the vertex in the graph (in Q).
        :return: True if the labeling is extendable, False otherwise.
        """
        if not theta:
            return False
        first_bag = self.first_appear[self.Q[i]]
        bag = self.nodes[first_bag]["bag"]
        frozen_theta = frozendict({key: theta[key] for key in bag})
        if frozen_theta in self.nodes[first_bag]["factor"].get_all_keys():
            return self.nodes[first_bag]["factor"].get_value(frozen_theta) == 1
        else:
            return False


if __name__ == '__main__':

    paper_graph = nx.Graph()
    x = 0
    paper_graph.add_nodes_from([i for i in range(x, x+7)])
    paper_graph.add_edges_from([(x, x+1),
                                (x+1, x+2),
                                (x+2, x+3),
                                (x+3, x+1),
                                (x+3, x+4),
                                (x+4, x+5),
                                (x+4, x+6)])

    rooted_dntd = RootedDisjointBranchNiceTreeDecomposition(paper_graph)
    rooted_dntd.calculate_factors_iterative()
    G1 = rooted_dntd.EnumDS(dict(), debug_flag=False)

    for s in rooted_dntd.EnumDS(dict(), debug_flag=False):
        print("DS - ", s)