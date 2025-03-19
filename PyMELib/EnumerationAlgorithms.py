from labels import *
from PyMELib.utils.labels_utils import *


def EnumMDS(self, theta: Dict[str, Label], i=0, debug_flag=False) -> None:
    """
    This algorithm means to enumerate all the minimal dominating sets of the graph.
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
            if original_c.in_rho:
                if c.in_rho and first_c.in_rho and original_c - Label.F_rho.R0 == c - Label.F_rho.R0 + first_c - Label.F_rho.R0:
                    if first_c == Label.F_rho.R1 and c == Label.F_rho.R1:
                        new_theta = self.IncrementLabeling2(theta, i, c)
                    elif first_c == Label.F_rho.R0 and c == Label.F_rho.R0:
                        new_theta = self.IncrementLabeling2(theta, i, c)
                    elif first_c == Label.F_rho.R0 and c == Label.F_rho.R1:
                        new_theta = self.IncrementLabeling2(theta, i, c)
                    else:
                        if debug_flag:
                            print("Not Valid Labeling")
                            print("-" * 20)
                        continue
                elif original_c == Label.F_rho.R1 and first_c == Label.F_rho.R1 and c == Label.F_omega.W0:
                    new_theta = self.IncrementLabeling2(theta, i, c)
                elif original_c == Label.F_rho.R2 and first_c == Label.F_omega.W0 and c == Label.F_rho.R2:
                    new_theta = self.IncrementLabeling2(theta, i, c)
                elif original_c == Label.F_rho.R2 and first_c == Label.F_rho.R2 and c == Label.F_omega.W0:
                    new_theta = self.IncrementLabeling2(theta, i, c)
                else:
                    if debug_flag:
                        print("Not Valid Labeling")
                        print("-" * 20)
                    continue
            elif original_c.in_sigma and first_c.in_sigma and c.in_sigma:
                if original_c == Label.F_sigma.SI and first_c == Label.F_sigma.SI and c == Label.F_sigma.SI:
                    new_theta = self.IncrementLabeling2(theta, i, c)
                elif original_c == Label.F_sigma.S0 and first_c == Label.F_sigma.S0 and c == Label.F_sigma.S0:
                    new_theta = self.IncrementLabeling2(theta, i, c)
                elif original_c == Label.F_sigma.S1 and first_c == Label.F_sigma.S1 and c == Label.F_sigma.S0:
                    new_theta = self.IncrementLabeling2(theta, i, c)
                elif original_c == Label.F_sigma.S1 and first_c ==Label.F_sigma.S0 and c == Label.F_sigma.S1:
                    new_theta = self.IncrementLabeling2(theta, i, c)
                elif original_c == Label.F_sigma.S1 and first_c == Label.F_sigma.S1 and c == Label.F_sigma.S1:
                    new_theta = self.IncrementLabeling2(theta, i, c)
                else:
                    if debug_flag:
                        print("Not Valid Labeling")
                        print("-" * 20)
                    continue
            elif original_c.in_omega and first_c.in_omega and c.in_omega:
                if original_c == Label.F_omega.W0 and first_c == Label.F_omega.W0 and c == Label.F_omega.W0:
                    new_theta = self.IncrementLabeling2(theta, i, c)
                elif original_c == Label.F_omega.W1 and first_c == Label.F_omega.W0 and c == Label.F_omega.W1:
                    new_theta = self.IncrementLabeling2(theta, i, c)
                elif original_c == Label.F_omega.W1 and first_c == Label.F_omega.W1 and c == Label.F_omega.W0:
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
                yield from self.EnumMDS(option, i + 1, debug_flag=debug_flag)


def EnumMHS(self, theta: Dict[str, Label], i=0, debug_flag=False) -> None:
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
            if original_c.in_rho:
                if c.in_rho and first_c.in_rho and original_c - Label.F_rho.R0 == c - Label.F_rho.R0 + first_c - Label.F_rho.R0:
                    if first_c == Label.F_rho.R1 and c == Label.F_rho.R1:
                        new_theta = self.IncrementLabeling2(theta, i, c)
                    elif first_c == Label.F_rho.R0 and c == Label.F_rho.R0:
                        new_theta = self.IncrementLabeling2(theta, i, c)
                    elif first_c == Label.F_rho.R0 and c == Label.F_rho.R1:
                        new_theta = self.IncrementLabeling2(theta, i, c)
                    else:
                        if debug_flag:
                            print("Not Valid Labeling")
                            print("-" * 20)
                        continue
                elif original_c == Label.F_rho.R1 and first_c == Label.F_rho.R1 and c == Label.F_omega.W0:
                    new_theta = self.IncrementLabeling2(theta, i, c)
                elif original_c == Label.F_rho.R2 and first_c == Label.F_omega.W0 and c == Label.F_rho.R2:
                    new_theta = self.IncrementLabeling2(theta, i, c)
                elif original_c == Label.F_rho.R2 and first_c == Label.F_rho.R2 and c == Label.F_omega.W0:
                    new_theta = self.IncrementLabeling2(theta, i, c)
                else:
                    if debug_flag:
                        print("Not Valid Labeling")
                        print("-" * 20)
                    continue
            elif original_c.in_sigma and first_c.in_sigma and c.in_sigma:
                if original_c == Label.F_sigma.SI and first_c == Label.F_sigma.SI and c == Label.F_sigma.SI:
                    new_theta = self.IncrementLabeling2(theta, i, c)
                elif original_c == Label.F_sigma.S0 and first_c == Label.F_sigma.S0 and c == Label.F_sigma.S0:
                    new_theta = self.IncrementLabeling2(theta, i, c)
                elif original_c == Label.F_sigma.S1 and first_c == Label.F_sigma.S1 and c == Label.F_sigma.S0:
                    new_theta = self.IncrementLabeling2(theta, i, c)
                elif original_c == Label.F_sigma.S1 and first_c == Label.F_sigma.S0 and c == Label.F_sigma.S1:
                    new_theta = self.IncrementLabeling2(theta, i, c)
                elif original_c == Label.F_sigma.S1 and first_c == Label.F_sigma.S1 and c == Label.F_sigma.S1:
                    new_theta = self.IncrementLabeling2(theta, i, c)
                else:
                    if debug_flag:
                        print("Not Valid Labeling")
                        print("-" * 20)
                    continue
            elif original_c.in_omega and first_c.in_omega and c.in_omega:
                if original_c == Label.F_omega.W0 and first_c == Label.F_omega.W0 and c == Label.F_omega.W0:
                    new_theta = self.IncrementLabeling2(theta, i, c)
                elif original_c == Label.F_omega.W1 and first_c == Label.F_omega.W0 and c == Label.F_omega.W1:
                    new_theta = self.IncrementLabeling2(theta, i, c)
                elif original_c == Label.F_omega.W1 and first_c == Label.F_omega.W1 and c == Label.F_omega.W0:
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
                yield from self.EnumMHS(option, i + 1, debug_flag=debug_flag)


def EnumMHS_iterative(self, debug_flag=False) -> None:
    """
    This is a for loop version of EnumMHS, using a stack.
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
                if original_c.in_rho:
                    if c.in_rho and first_c.in_rho and original_c - Label.F_rho.R0 == c - Label.F_rho.R0 + first_c - Label.F_rho.R0:
                        if first_c == Label.F_rho.R1 and c == Label.F_rho.R1:
                            new_theta = self.IncrementLabeling2(theta, i, c)
                        elif first_c == Label.F_rho.R0 and c == Label.F_rho.R0:
                            new_theta = self.IncrementLabeling2(theta, i, c)
                        elif first_c == Label.F_rho.R0 and c == Label.F_rho.R1:
                            new_theta = self.IncrementLabeling2(theta, i, c)
                        else:
                            if debug_flag:
                                print("Not Valid Labeling")
                                print("-" * 20)
                            continue
                    elif original_c == Label.F_rho.R1 and first_c == Label.F_rho.R1 and c == Label.F_omega.W0:
                        new_theta = self.IncrementLabeling2(theta, i, c)
                    elif original_c == Label.F_rho.R2 and first_c == Label.F_omega.W0 and c == Label.F_rho.R2:
                        new_theta = self.IncrementLabeling2(theta, i, c)
                    elif original_c == Label.F_rho.R2 and first_c == Label.F_rho.R2 and c == Label.F_omega.W0:
                        new_theta = self.IncrementLabeling2(theta, i, c)
                    else:
                        if debug_flag:
                            print("Not Valid Labeling")
                            print("-" * 20)
                        continue
                elif original_c.in_sigma and first_c.in_sigma and c.in_sigma:
                    if original_c == Label.F_sigma.SI and first_c == Label.F_sigma.SI and c == Label.F_sigma.SI:
                        new_theta = self.IncrementLabeling2(theta, i, c)
                    elif original_c == Label.F_sigma.S0 and first_c == Label.F_sigma.S0 and c == Label.F_sigma.S0:
                        new_theta = self.IncrementLabeling2(theta, i, c)
                    elif original_c == Label.F_sigma.S1 and first_c == Label.F_sigma.S1 and c == Label.F_sigma.S0:
                        new_theta = self.IncrementLabeling2(theta, i, c)
                    elif original_c == Label.F_sigma.S1 and first_c == Label.F_sigma.S0 and c == Label.F_sigma.S1:
                        new_theta = self.IncrementLabeling2(theta, i, c)
                    elif original_c == Label.F_sigma.S1 and first_c == Label.F_sigma.S1 and c == Label.F_sigma.S1:
                        new_theta = self.IncrementLabeling2(theta, i, c)
                    else:
                        if debug_flag:
                            print("Not Valid Labeling")
                            print("-" * 20)
                        continue
                elif original_c.in_omega and first_c.in_omega and c.in_omega:
                    if original_c == Label.F_omega.W0 and first_c == Label.F_omega.W0 and c == Label.F_omega.W0:
                        new_theta = self.IncrementLabeling2(theta, i, c)
                    elif original_c == Label.F_omega.W1 and first_c == Label.F_omega.W0 and c == Label.F_omega.W1:
                        new_theta = self.IncrementLabeling2(theta, i, c)
                    elif original_c == Label.F_omega.W1 and first_c == Label.F_omega.W1 and c == Label.F_omega.W0:
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
    Procedure IncrementLabeling receives as input a labeling which we assume to be extendable (see EnumMDS),
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


def IncrementLabeling(self, theta: Dict[str, Label], i, c: Label, V_label_S, V_label_W):
    """
    Procedure IncrementLabeling receives as input a labeling which we assume to be extendable (see EnumMDS),
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
    K_i = {x + self.nodes[self.first_appear[current_vertex]]["br"] for x in K_i}
    N_i = K_i.intersection(V_label_S)
    W_i = K_i.intersection(V_label_W)

    flag_of_two = False

    if c.in_sigma:
        for v in K_i:
            if theta[v].in_rho:
                new_theta[v] = max(Label.F_rho.R0, theta[v] - 1)

    if c == Label.F_sigma.SI and (len(N_i) != 0 or len(W_i) != 0):
        return False
    if Label.F_sigma.S0 <= c <= Label.F_sigma.S1:
        if len([w for w in K_i if theta[w] in {Label.F_sigma.SI,Label.F_omega.W0}]) != 0 or \
                (c == Label.F_sigma.S0 and len(W_i) == 0):
            return False
        else:
            for w in W_i:
                if theta[w] == Label.F_omega.W1:
                    new_theta[w] = Label.F_omega.W0
    if c.in_omega:
        if len([w for w in N_i if theta[w] == Label.F_sigma.SI]) != 0 or \
                len(N_i) >= 2 or \
                (len(N_i) == 0 and c == Label.F_omega.W0) or \
                (len(N_i) != 0 and c == Label.F_omega.W1):
            return False
        elif c == Label.F_omega.W0:
            v = N_i.pop()
            if theta[v] == Label.F_sigma.S0:
                return False
            flag_of_two = v
    if c and max(0, 2 - len(N_i)) != c - Label.F_rho.R0:
        return False
    if flag_of_two:
        new_theta[flag_of_two] = Label.F_sigma.S0
        new_theta2 = dict(new_theta)
        new_theta2[flag_of_two] = Label.F_sigma.S1
        return [new_theta, new_theta2]
    return [new_theta]