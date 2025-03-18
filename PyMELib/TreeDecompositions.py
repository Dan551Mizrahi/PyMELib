import networkx as nx
from enum import IntEnum
import plotly.graph_objects as go
import EoN


# The code is based on the paper "Enumeration of minimal hitting sets parameterized by treewidth" by Batya Kenig and Dan Shlomo Mizrahi
# This code was created by Dan Mizrahi with the help of GitHub Copilot and Gemini.
# The code supports graphs up to 1,114,112 nodes (maximum chr value in python).

class NodeType(IntEnum):
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

    def get_original_root(self):
        """
        Gets the original root node of the tree decomposition.
        """
        return self.original_root

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

    def draw(self, save_path: str = None) -> None:
        """
        Draws the rooted tree decomposition using a hierarchical layout. This visualization includes the bags of the nodes.
        :param save_path: (Optional) The path to save the visualization as an HTML file (don't write the .html ending).
        """
        pos = EoN.hierarchy_pos(self, root=self.root)
        nx.draw(self, pos, with_labels=True)

        # Extract node information
        node_x = []
        node_y = []
        node_text = []
        hover_text = []
        for node, (x, y) in pos.items():
            node_x.append(x)
            node_y.append(y)
            node_text.append(self.nodes[node]['bag'])
            hover_text.append(f"ID: {node}<br>")

        # Create Plotly trace for nodes
        node_trace = go.Scatter(
            x=node_x, y=node_y,
            mode='markers+text',
            hoverinfo='text',
            text=node_text,
            hovertext=hover_text,
            textposition='bottom center',
            marker=dict(showscale=False, symbol='circle-dot', size=20 )
        )

        # Create Plotly trace for edges
        edge_x = []
        edge_y = []
        for edge in self.edges():
            x0, y0 = pos[edge[0]]
            x1, y1 = pos[edge[1]]
            edge_x.extend([x0, x1, None])
            edge_y.extend([y0, y1, None])

        edge_trace = go.Scatter(
            x=edge_x, y=edge_y,
            line=dict(width=0.5, color='#888'),
            hoverinfo='none',
            mode='lines')

        # Create the figure
        fig = go.Figure(data=[edge_trace, node_trace],
                        layout=go.Layout(
                            showlegend=False,
                            hovermode='closest',
                            margin=dict(b=0, l=0, r=0, t=0)
                        ))
        if save_path:
            fig.to_html(save_path + ".html")
        else:
            fig.show()


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
        self.root = new_root
        if semi_nice:
            self.transform_to_semi_nice_rec(self.get_root())
        else:
            self.transform_to_nice_rec(self.get_root())

        self.nodes[new_root]["type"] = NodeType.ROOT

        # create a complete order of the vertices (for enumeration process)
        self.Q = []
        self.create_Q(self.get_root())

    def create_Q(self, current_node):
        if self.nodes[current_node]["type"] == NodeType.LEAF:
            return
        if self.nodes[current_node]["type"] == FORGET:
            v = self.nodes[list(self.successors(current_node))[0]]["bag"].difference(
                self.nodes[current_node]["bag"]).pop()
            self.Q.append(v)
        elif self.nodes[current_node]["type"] != NodeType.ROOT and self.nodes[list(self.predecessors(current_node))[0]][
            "type"] == NodeType.ROOT:
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
                self.nodes[current_node]["type"] = NodeType.LEAF

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
                self.nodes[current_node]["type"] = NodeType.LEAF

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

'''
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
        if self.nodes[current_node]["type"] == NodeType.LEAF:
            return
        for vertex in self.nodes[current_node]["bag"]:
            if vertex not in self.first_appear.keys() or self.first_appear[vertex] is None:
                self.first_appear[vertex] = current_node
        for child in self.successors(current_node):
            self.first_appear_update(child)

    def create_Q(self, current_node):
        if self.nodes[current_node]["type"] == NodeType.LEAF:
            return
        if self.nodes[current_node]["type"] != NodeType.ROOT and self.nodes[list(self.predecessors(current_node))[0]][
            "type"] == NodeType.ROOT:
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

        if self.nodes[current_node]["type"] == NodeType.LEAF:
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
            elif self.nodes[current_node]["type"] == FORGET or self.nodes[current_node]["type"] == NodeType.ROOT or \
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

        if self.nodes[current_node]["type"] == NodeType.LEAF:
            return

        children = list(self.successors(current_node))
        if self.nodes[current_node]["type"] == NodeType.ROOT:
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

        if self.nodes[current_node]["type"] == NodeType.LEAF:
            return

        children = list(self.successors(current_node))
        if self.nodes[current_node]["type"] == NodeType.ROOT:
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
'''

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

    td = RootedTreeDecomposition(paper_graph)
    td.draw()
    # rooted_dntd = RootedDisjointBranchNiceTreeDecomposition(paper_graph)
    # rooted_dntd.calculate_factors_for_mds_enum_iterative()
    # G1 = rooted_dntd.EnumMDS(dict(), debug_flag=False)
    #
    # for s in rooted_dntd.EnumMDS(dict(), debug_flag=False):
    #     print("DS - ", s)