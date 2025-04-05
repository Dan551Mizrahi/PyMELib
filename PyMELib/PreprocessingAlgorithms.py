from PyMELib.TreeDecompositions import RootedDisjointBranchNiceTreeDecomposition, NodeType
from PyMELib.Factors import MemoTable
from frozendict import frozendict
from PyMELib.utils.labels_utils import *
from PyMELib.utils.comb_utils import generate_dictionaries_from_sets

def create_factors(td: RootedDisjointBranchNiceTreeDecomposition) -> None:
    """
    Create empty factors for the nodes of the TD before the preprocessing phase.
    :param td: RootedDisjointBranchNiceTreeDecomposition
    :return: None
    """
    for node in td.nodes:
        td.nodes[node]["factor"] = MemoTable()

def calculate_factors_for_mds_enum(td: RootedDisjointBranchNiceTreeDecomposition, current_node: int, options_for_labels=False) -> None:
    """
    This is a dynamic programming algorithm that calculates the factors of the TD.
    In order to run the algorithm for enumeration of dominating sets.
    This what we call in our paper: preprocessing phase.

    :param td: The disjoint branch nice tree decomposition we are working on.
    :param current_node: The current node that we are on TD.
    :param options_for_labels: If True, we consider the options for labels of the vertices.
    :return: None
    """

    def calculate_k_range(vertex: str, label_lower: int, label_upper: int, assignment: dict):

        # TODO: why not local neighbors?
        neighbors = {chr(n) for n in td.original_graph.neighbors(ord(vertex[0]))}
        v_label_set = {w1[0] for w1, l1 in assignment.items() if label_lower <= l1 <= label_upper}

        return len(neighbors.intersection(v_label_set))

    def calculate_k(vertex: str, label_a: str, assignment: dict):

        if label_a == "S":
            return calculate_k_range(vertex, F_sigma.SI, F_sigma.S1, assignment)
        elif label_a == "W":
            return calculate_k_range(vertex, F_omega.W0, F_omega.W1, assignment)
        elif label_a == "R":
            return calculate_k_range(vertex, F_rho.R0, F_rho.R2, assignment)
        else:
            return calculate_k_range(vertex, trans_dict[label_a], trans_dict[label_a], assignment)

    type_of_node = td.nodes[current_node]["type"]
    children = list(td.successors(current_node))

    if type_of_node != NodeType.LEAF:
        for child in children:
            calculate_factors_for_mds_enum(td, child, options_for_labels=options_for_labels)
        child_node = children[0]

    if type_of_node == NodeType.LEAF:
        td.nodes[current_node]["factor"].set_value(frozendict(), 1)
    elif type_of_node == NodeType.BIG_JOIN_INTRODUCE:

        set_of_real_vertices = {v[0] + td.nodes[current_node]["br"] for v in td.nodes[current_node]["bag"]}
        dict_of_copies = {v: [] for v in set_of_real_vertices}

        for v in set_of_real_vertices:
            if (v + "0") in td.nodes[current_node]["bag"]:
                dict_of_copies[v].append(v + "0")
            if (v + "1") in td.nodes[current_node]["bag"]:
                dict_of_copies[v].append(v + "1")

        # the new vertex that was introduced
        v = td.nodes[current_node]["bag"].difference(td.nodes[child_node]["bag"]).pop()

        for key in td.nodes[child_node]["factor"].get_all_keys():

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
            if options_for_labels and not (phi[v] in td.original_graph.nodes[ord(v[0])]["options"]):
                flag = False
            if flag:
                td.nodes[current_node]["factor"].set_value(frozendict(phi), 1)

    elif type_of_node == NodeType.JOIN_INTRODUCE:

        # the new vertex that was introduced
        v = td.nodes[current_node]["bag"].difference(td.nodes[child_node]["bag"]).pop()

        for key in td.nodes[child_node]["factor"].get_all_keys():
            old_value = td.nodes[child_node]["factor"].get_value(key)
            if old_value == 0:
                continue
            for label in F:
                if options_for_labels and not (label in td.original_graph.nodes[ord(v[0])]["options"]):
                    continue
                phi = dict(key)
                phi[v] = label
                td.nodes[current_node]["factor"].set_value(frozendict(phi), 1)

    elif type_of_node == NodeType.INTRODUCE:

        # the new vertex that was introduced
        v = td.nodes[current_node]["bag"].difference(td.nodes[child_node]["bag"]).pop()

        for key in td.nodes[child_node]["factor"].get_all_keys():
            if td.nodes[child_node]["type"] != NodeType.LEAF and key == frozendict():  # TODO: why is it even possible?
                print(f"Error: empty key at node {td.nodes[child_node]['bag']}")
                exit(-1)
            old_value = td.nodes[child_node]["factor"].get_value(key)
            if old_value == 0:
                continue
            phi = dict(key)
            for label in [F_sigma.SI, F_rho.R0, F_omega.W0, F_sigma.S0]:
                if options_for_labels and not (label in td.original_graph.nodes[ord(v[0])]["options"]):
                    continue
                phi[v] = label
                if label == F_sigma.SI:
                    k_v_s = calculate_k(v, "S", key)
                    if k_v_s == 0:
                        td.nodes[current_node]["factor"].set_value(frozendict(phi), 1)
                    else:
                        td.nodes[current_node]["factor"].set_value(frozendict(phi), 0)
                elif label == F_sigma.S0:
                    k_v_si = calculate_k(v, "SI", key)
                    if k_v_si == 0:
                        td.nodes[current_node]["factor"].set_value(frozendict(phi), 1)
                    else:
                        td.nodes[current_node]["factor"].set_value(frozendict(phi), 0)
                else:
                    td.nodes[current_node]["factor"].set_value(frozendict(phi), 1)

    elif type_of_node == NodeType.JOIN_FORGET:

        v = td.nodes[child_node]["bag"].difference(td.nodes[current_node]["bag"]).pop()

        for key in td.nodes[child_node]["factor"].get_all_keys():
            # Just copy the factor
            old_value = td.nodes[child_node]["factor"].get_value(key)
            if old_value == 0:
                continue
            new_key = dict(key)
            del new_key[v]
            new_key = frozendict(new_key)
            td.nodes[current_node]["factor"].set_value(new_key, 1)

    elif type_of_node == NodeType.FORGET or type_of_node == NodeType.ROOT:

        v = td.nodes[child_node]["bag"].difference(td.nodes[current_node]["bag"]).pop()

        # All possible keys pass
        for key in td.nodes[child_node]["factor"].get_all_keys():
            old_value = td.nodes[child_node]["factor"].get_value(key)
            if old_value == 0:
                continue

            v_label = key[v]
            k_v_si = calculate_k(v, "SI", key)
            k_v_s = calculate_k(v, "S", key)
            k_v_w = calculate_k(v, "W", key)
            k_v_w1 = calculate_k(v, "W1", key)
            k_v_w0 = calculate_k(v, "W0", key)

            if v_label == F_sigma.SI:
                if not (k_v_s == 0 and k_v_w == 0):
                    td.nodes[child_node]["factor"].set_value(key, 0)
                    continue
            elif v_label.in_rho:
                j = v_label - F_rho.R0 # TODO: check id still holds
                if not (j + k_v_s >= 2):
                    td.nodes[child_node]["factor"].set_value(key, 0)
                    continue
            elif v_label == F_omega.W1:
                if not (k_v_s == 0):
                    td.nodes[child_node]["factor"].set_value(key, 0)
                    continue
            elif v_label == F_omega.W0:
                if not (k_v_si == 0 and k_v_s == 1):
                    td.nodes[child_node]["factor"].set_value(key, 0)
                    continue
            elif v_label == F_sigma.S1:
                if not (k_v_si == 0 and k_v_w1 == 0):
                    td.nodes[child_node]["factor"].set_value(key, 0)
                    continue
            elif v_label == F_sigma.S0:
                if not (k_v_si == 0 and k_v_w1 == 0 and k_v_w0 >= 1):
                    td.nodes[child_node]["factor"].set_value(key, 0)
                    continue

            options_for_new_labels = {z: set() for z in td.nodes[current_node]["bag"]}

            flag = True
            for w in key.keys():
                if w == v:
                    continue
                if v_label == F_sigma.SI:
                    if not (ord(w[0]) in td.original_graph.neighbors(ord(v[0]))):
                        options_for_new_labels[w].add(key[w])
                    elif key[w].in_rho:
                        for j in range(3):
                            if max(0, j - 1) == key[w] - F_rho.R0:
                                options_for_new_labels[w].add(trans_dict["R" + str(j)])
                                # TODO: can't there be break here?
                    else:
                        flag = False
                        break
                elif v_label == F_sigma.S0:
                    if not (ord(w[0]) in td.original_graph.neighbors(ord(v[0]))) or (F_sigma.S0 <= key[w] <= F_sigma.S1):
                        options_for_new_labels[w].add(key[w])
                    elif key[w].in_rho:
                        for j in range(3):
                            if max(0, j - 1) == key[w] - F_rho.R0:
                                options_for_new_labels[w].add(trans_dict["R" + str(j)])
                    elif key[w] == F_omega.W0:
                        options_for_new_labels[w].add(F_omega.W1)
                    else:
                        flag = False
                        break
                elif v_label == F_sigma.S1:
                    if not (ord(w[0]) in td.original_graph.neighbors(ord(v[0]))) or (F_sigma.S0 <= key[w] <= F_sigma.S1):
                        options_for_new_labels[w].add(key[w])
                    elif key[w].in_rho:
                        for j in range(3):
                            if max(0, j - 1) == key[w] - F_rho.R0:
                                options_for_new_labels[w].add(trans_dict["R" + str(j)])
                    elif key[w] == F_omega.W0:
                        options_for_new_labels[w].add(F_omega.W1)
                    else:
                        flag = False
                        break
                elif v_label == F_omega.W1:
                    if not (ord(w[0]) in td.original_graph.neighbors(ord(v[0]))) or not key[w].in_sigma:
                        options_for_new_labels[w].add(key[w])
                    else:
                        flag = False
                        break
                elif v_label == F_omega.W0:
                    if not (ord(w[0]) in td.original_graph.neighbors(ord(v[0]))) or not key[w].in_sigma:
                        options_for_new_labels[w].add(key[w])
                    elif F_sigma.S0 <= key[w] <= F_sigma.S1:
                        options_for_new_labels[w].add(F_sigma.S1)
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
                    td.nodes[current_node]["factor"].set_value(frozendict(opt), 1)


    else:
        child_node_1 = children[0]
        child_node_2 = children[1]

        # taking a simple conjunction of the two factors
        # Here is a kind of natural join algorithm on the two factors, where we check if eah node has the same label
        # in both factors, and if so on all the nodes in that key, we add the combine key it to the new factor.

        for key_1 in td.nodes[child_node_1]["factor"].get_all_keys():
            if td.nodes[child_node_1]["factor"].get_value(key_1) == 0:
                continue
            for key_2 in td.nodes[child_node_2]["factor"].get_all_keys():
                if td.nodes[child_node_2]["factor"].get_value(key_2) == 0:
                    continue
                flag = True
                new_key = {z: 0 for z in td.nodes[current_node]["bag"]}
                for w in key_1.keys():
                    label1 = key_1[w]
                    w2 = w[:-1] + str((int(w[-1]) + 1) % 2)
                    if w2 in key_2.keys():
                        label2 = key_2[w2]
                        if ((label1.same_class(label2)) or (label1.in_omega and label2.in_rho) or
                                (label1.in_rho and label2.in_omega)):
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
                    td.nodes[current_node]["factor"].set_value(frozendict(new_key), 1)

    if type_of_node == NodeType.ROOT:
        return


def calculate_factors_for_mds_enum_iterative(td: RootedDisjointBranchNiceTreeDecomposition, options_for_labels=False):
    """
    This is a for loop version of calculate_factors_for_mds_enum, using a stack.
    """

    def calculate_k_range(vertex: str, label_lower: int, label_upper: int, assignment: dict):

        # TODO: why not local neighbors?
        neighbors = {chr(n) for n in td.original_graph.neighbors(ord(vertex[0]))}
        v_label_set = {w1[0] for w1, l1 in assignment.items() if label_lower <= l1 <= label_upper}

        return len(neighbors.intersection(v_label_set))

    def calculate_k(vertex: str, label_a: str, assignment: dict):

        if label_a == "S":
            return calculate_k_range(vertex, F_sigma.SI, F_sigma.S1, assignment)
        elif label_a == "W":
            return calculate_k_range(vertex, F_omega.W0, F_omega.W1, assignment)
        elif label_a == "R":
            return calculate_k_range(vertex, F_rho.R0, F_rho.R2, assignment)
        else:
            return calculate_k_range(vertex, trans_dict[label_a], trans_dict[label_a], assignment)

    stack = [td.get_root()]
    while len(stack) != 0:

        current_node = stack[-1]

        type_of_node = td.nodes[current_node]["type"]
        children = list(td.successors(current_node))

        # Check if children have been processed
        children_processed = True
        if type_of_node != NodeType.LEAF:
            for child in children:
                if "processed" not in td.nodes[child] or not td.nodes[child]["processed"]:
                    children_processed = False
                    stack.append(child)
            child_node = children[0]

        if not children_processed:
            continue  # Continue processing children

        # All children processed (or leaf node), process current node
        stack.pop()  # Remove the current node from the stack

        if type_of_node == NodeType.LEAF:
            td.nodes[current_node]["factor"].set_value(frozendict(), 1)
        elif type_of_node == NodeType.BIG_JOIN_INTRODUCE:

            set_of_real_vertices = {v[0] + td.nodes[current_node]["br"] for v in td.nodes[current_node]["bag"]}
            dict_of_copies = {v: [] for v in set_of_real_vertices}

            for v1 in sorted(td.nodes[current_node]["bag"]):
                for v2 in set_of_real_vertices:
                    if v2 != v1 and v2.startswith(v1[0]):
                        if td.is_semi_nice and len(v1) < len(v2):
                            del dict_of_copies[v2]
                            dict_of_copies[v1] = []
                            dict_of_copies[v1].append(v1)
                        else:
                            if v2 in dict_of_copies.keys():
                                dict_of_copies[v2].append(v1)
                            else:
                                new_key = None
                                for key in dict_of_copies.keys():
                                    if key.startswith(v1[0]) and len(key) < len(v2):
                                        new_key = key
                                        break
                                dict_of_copies[new_key].append(v1)
            for key, list1 in dict_of_copies.items():
                if len(list1) == 3:
                    dict_of_copies[key].remove(key)

            if td.is_semi_nice:
                set_of_real_vertices = set(dict_of_copies.keys())

            # the new vertex that was introduced
            v = td.nodes[current_node]["bag"].difference(td.nodes[child_node]["bag"]).pop()

            for key in td.nodes[child_node]["factor"].get_all_keys():

                phi = dict(key)
                phi[v] = "N"
                flag = True

                for original_vertex in dict_of_copies.keys():
                    if flag is False:
                        break
                    if len(dict_of_copies[original_vertex]) == 0 or len(dict_of_copies[original_vertex]) == 1:
                        continue

                    first_copy = dict_of_copies[original_vertex][0]
                    second_copy = dict_of_copies[original_vertex][1]

                    label_0 = phi[first_copy]
                    label_1 = phi[second_copy]

                    flag, label_original_vertex = join_labels(label_0, label_1)
                    if not flag:
                        continue
                    phi[original_vertex] = label_original_vertex
                if options_for_labels and not (phi[v] in td.original_graph.nodes[ord(v[0])]["options"]):
                    flag = False
                if flag:
                    td.nodes[current_node]["factor"].set_value(frozendict(phi), 1)

        elif type_of_node == NodeType.JOIN_INTRODUCE:

            # the new vertex that was introduced
            v = td.nodes[current_node]["bag"].difference(td.nodes[child_node]["bag"]).pop()

            for key in td.nodes[child_node]["factor"].get_all_keys():
                old_value = td.nodes[child_node]["factor"].get_value(key)
                if old_value == 0:
                    continue
                for label in F:
                    if options_for_labels and not (label in td.original_graph.nodes[ord(v[0])]["options"]):
                        continue
                    phi = dict(key)
                    phi[v] = label
                    td.nodes[current_node]["factor"].set_value(frozendict(phi), 1)

        elif type_of_node == NodeType.INTRODUCE:

            # the new vertex that was introduced
            v = td.nodes[current_node]["bag"].difference(td.nodes[child_node]["bag"]).pop()

            for key in td.nodes[child_node]["factor"].get_all_keys():
                if td.nodes[child_node]["type"] != NodeType.LEAF and key == frozendict():  # TODO: why is it even possible?
                    print(f"Error: empty key at node {td.nodes[child_node]['bag']}")
                    exit(-1)
                old_value = td.nodes[child_node]["factor"].get_value(key)
                if old_value == 0:
                    continue
                phi = dict(key)
                for label in [F_sigma.SI, F_rho.R0, F_omega.W0, F_sigma.S0]:
                    if options_for_labels and not (label in td.original_graph.nodes[ord(v[0])]["options"]):
                        continue
                    phi[v] = label
                    if label == F_sigma.SI:
                        k_v_s = calculate_k(v, "S", key)
                        if k_v_s == 0:
                            td.nodes[current_node]["factor"].set_value(frozendict(phi), 1)
                        else:
                            td.nodes[current_node]["factor"].set_value(frozendict(phi), 0)
                    elif label == F_sigma.S0:
                        k_v_si = calculate_k(v, "SI", key)
                        if k_v_si == 0:
                            td.nodes[current_node]["factor"].set_value(frozendict(phi), 1)
                        else:
                            td.nodes[current_node]["factor"].set_value(frozendict(phi), 0)
                    else:
                        td.nodes[current_node]["factor"].set_value(frozendict(phi), 1)

        elif type_of_node == NodeType.JOIN_FORGET:

            v = td.nodes[child_node]["bag"].difference(td.nodes[current_node]["bag"]).pop()

            for key in td.nodes[child_node]["factor"].get_all_keys():
                # Just copy the factor
                old_value = td.nodes[child_node]["factor"].get_value(key)
                if old_value == 0:
                    continue
                new_key = dict(key)
                del new_key[v]
                new_key = frozendict(new_key)
                td.nodes[current_node]["factor"].set_value(new_key, 1)

        elif type_of_node == NodeType.FORGET or type_of_node == NodeType.ROOT:

            # print("current node bag: ", {ord(xx[0]) for xx in td.nodes[current_node]["bag"]})
            # print("child node bag: ", {ord(xx[0]) for xx in td.nodes[child_node]["bag"]})
            # print("Type of child node: ", td.nodes[child_node]["type"].name)

            v = td.nodes[child_node]["bag"].difference(td.nodes[current_node]["bag"]).pop()
            # print("v :", ord(v[0]))

            # All possible keys pass
            for key in td.nodes[child_node]["factor"].get_all_keys():
                old_value = td.nodes[child_node]["factor"].get_value(key)
                if old_value == 0:
                    continue

                v_label = key[v]
                k_v_si = calculate_k(v, "SI", key)
                k_v_s = calculate_k(v, "S", key)
                k_v_w = calculate_k(v, "W", key)
                k_v_w1 = calculate_k(v, "W1", key)
                k_v_w0 = calculate_k(v, "W0", key)

                if v_label == F_sigma.SI:
                    if not (k_v_s == 0 and k_v_w == 0):
                        td.nodes[child_node]["factor"].set_value(key, 0)
                        continue
                elif v_label.in_rho:
                    j = v_label - F_rho.R0
                    if not (j + k_v_s >= 2):
                        td.nodes[child_node]["factor"].set_value(key, 0)
                        continue
                elif v_label == F_omega.W1:
                    if not (k_v_s == 0):
                        td.nodes[child_node]["factor"].set_value(key, 0)
                        continue
                elif v_label == F_omega.W0:
                    if not (k_v_si == 0 and k_v_s == 1):
                        td.nodes[child_node]["factor"].set_value(key, 0)
                        continue
                elif v_label == F_sigma.S1:
                    if not (k_v_si == 0 and k_v_w1 == 0):
                        td.nodes[child_node]["factor"].set_value(key, 0)
                        continue
                elif v_label == F_sigma.S0:
                    if not (k_v_si == 0 and k_v_w1 == 0 and k_v_w0 >= 1):
                        td.nodes[child_node]["factor"].set_value(key, 0)
                        continue

                options_for_new_labels = {z: set() for z in td.nodes[current_node]["bag"]}

                flag = True
                for w in key.keys():
                    if w == v:
                        continue
                    if v_label == F_sigma.SI:
                        if not (ord(w[0]) in td.original_graph.neighbors(ord(v[0]))):
                            options_for_new_labels[w].add(key[w])
                        elif key[w].in_rho:
                            for j in range(3):
                                if max(0, j - 1) == key[w] - F_rho.R0:
                                    options_for_new_labels[w].add(trans_dict["R" + str(j)])
                                    # TODO: can't there be break here?
                        else:
                            flag = False
                            break
                    elif v_label == F_sigma.S0:
                        if not (ord(w[0]) in td.original_graph.neighbors(ord(v[0]))) or (F_sigma.S0 <= key[w] <= F_sigma.S1):
                            options_for_new_labels[w].add(key[w])
                        elif key[w].in_rho:
                            for j in range(3):
                                if max(0, j - 1) == key[w] - F_rho.R0:
                                    options_for_new_labels[w].add(trans_dict["R" + str(j)])
                        elif key[w] == F_omega.W0:
                            options_for_new_labels[w].add(F_omega.W1)
                        else:
                            flag = False
                            break
                    elif v_label == F_sigma.S1:
                        if not (ord(w[0]) in td.original_graph.neighbors(ord(v[0]))) or (F_sigma.S0 <= key[w] <= F_sigma.S1):
                            options_for_new_labels[w].add(key[w])
                        elif key[w].in_rho:
                            for j in range(3):
                                if max(0, j - 1) == key[w] - F_rho.R0:
                                    options_for_new_labels[w].add(trans_dict["R" + str(j)])
                        elif key[w] == F_omega.W0:
                            options_for_new_labels[w].add(F_omega.W1)
                        else:
                            flag = False
                            break
                    elif v_label == F_omega.W1:
                        if not (ord(w[0]) in td.original_graph.neighbors(ord(v[0]))) or not key[w].in_sigma:
                            options_for_new_labels[w].add(key[w])
                        else:
                            flag = False
                            break
                    elif v_label == F_omega.W0:
                        if not (ord(w[0]) in td.original_graph.neighbors(ord(v[0]))) or not key[w].in_sigma:
                            options_for_new_labels[w].add(key[w])
                        elif F_sigma.S0 <= key[w] <= F_sigma.S1:
                            options_for_new_labels[w].add(F_sigma.S1)
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
                        td.nodes[current_node]["factor"].set_value(frozendict(opt), 1)


        elif type_of_node == NodeType.JOIN:
            child_node_1 = children[0]
            child_node_2 = children[1]

            # taking a simple conjunction of the two factors
            # Here is a kind of natural join algorithm on the two factors, where we check if eah node has the same label
            # in both factors, and if so on all the nodes in that key, we add the combine key it to the new factor.

            for key_1 in td.nodes[child_node_1]["factor"].get_all_keys():
                if td.nodes[child_node_1]["factor"].get_value(key_1) == 0:
                    continue
                for key_2 in td.nodes[child_node_2]["factor"].get_all_keys():
                    if td.nodes[child_node_2]["factor"].get_value(key_2) == 0:
                        continue
                    flag = True
                    new_key = {z: 0 for z in td.nodes[current_node]["bag"]}
                    for w in key_1.keys():
                        label1 = key_1[w]
                        if not td.is_semi_nice and len(w[1:]) == len(td.nodes[current_node]["br"]) + 1:
                            w2 = w[:-1] + str((int(w[-1]) + 1) % 2)
                            if w2 in key_2.keys():
                                label2 = key_2[w2]
                                if ((label1.same_class(label2)) or (label1.in_omega and label2.in_rho) or
                                        (label1.in_rho and label2.in_omega)):
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
                        td.nodes[current_node]["factor"].set_value(frozendict(new_key), 1)

        # Mark the current node as processed
        td.nodes[current_node]["processed"] = True
