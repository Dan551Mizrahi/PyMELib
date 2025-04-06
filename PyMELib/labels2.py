from typing import Dict

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

def same_class(l0: int, l1: int):
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
reverse_trans_dict = {v:k for k,v in trans_dict.items()}

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
