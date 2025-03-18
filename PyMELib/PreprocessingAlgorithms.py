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


def calculate_factors_for_mds_enum(self, current_node, with_options=False):
    """
    This is a dynamic programming algorithm that calculates the factors of the TD.
    In order to implement the algorithm for enumeration of dominating sets.
    This what we call in our paper: pre-processing phase.

    :param current_node: The current node that we are on TD.
    :return: None
    """

    def calculate_k_range(vertex: str, label_lower: int, label_upper: int, assignment: dict):

        # TODO: why not local neighbors?
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
            self.calculate_factors_for_mds_enum(child)
        child_node = children[0]

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


def calculate_factors_for_mds_enum_iterative(self, with_options=False):
    """
    This is a for loop version of calculate_factors_for_mds_enum, using a stack.
    """

    def calculate_k_range(vertex: str, label_lower: int, label_upper: int, assignment: dict):

        # TODO: why not local neighbors?
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



