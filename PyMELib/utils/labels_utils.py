from PyMELib.labels import *
from typing import Tuple, Dict
from PyMELib.utils.comb_utils import reduce_dict_by_function

def join_labels(l0: Label, l1: Label) -> Tuple[bool, Label]:
    """
    Join two labels of vertex's copies to get the original vertex label.
    :param l0: Label of the first vertex copy.
    :param l1: Label of the second vertex copy.
    :return: boolean indicating whether the labels can be joined and the label of the original vertex if so.
    """
    flag = True
    label_original_vertex = None

    sum_of_labels = l0 + l1

    if l0.same_class(l1):
        if l0.in_sigma:
            if sum_of_labels == 0:
                label_original_vertex = F_sigma.SI
            elif sum_of_labels == 1:
                flag = False
            elif sum_of_labels == 2:
                if (l0 == F_sigma.SI and l1 == F_sigma.S1) or (l0 == F_sigma.S1 and l1 == F_sigma.SI):
                    label_original_vertex = F_sigma.S1
                else:
                    label_original_vertex = F_sigma.S0
            else:
                label_original_vertex = F_sigma.S1
        elif l0.in_omega:
            if sum_of_labels > 7:
                flag = False
            elif sum_of_labels == 7:
                label_original_vertex = F_omega.W1
            else:
                label_original_vertex = F_omega.W0
        elif l0.in_rho:
            if sum_of_labels == 10:
                label_original_vertex = F_rho.R0
            elif sum_of_labels == 11:
                label_original_vertex = F_rho.R1
            else:
                label_original_vertex = F_rho.R2

    elif l0.in_rho and l1 == F_omega.W0:
        label_original_vertex = l0
    elif l0 == F_omega.W0  and l1.in_rho:
        label_original_vertex = l1
    else:
        flag = False

    return flag, label_original_vertex

def V_label(label_a: str, theta: Dict[str, Label]) -> set:
    """
    Calculate the set of vertices with a specific label.
    :param label_a: Name of a label to filter the vertices (could be a class of labels: sigma ("S"), omega ("W"), rho ("R")).
    :param theta: Dictionary of vertices and their labels.
    :return: Set of vertices with the specified label.
    """
    if label_a == "S":
        return reduce_dict_by_function(theta, Label.in_sigma)
    elif label_a == "W":
        return reduce_dict_by_function(theta, Label.in_omega)
    elif label_a == "R":
        return reduce_dict_by_function(theta, Label.in_rho)
    else:
        return reduce_dict_by_function(theta, lambda x: x.name == label_a)

def V_label_S_W(theta: Dict[str, Label]):
    """
    Optimization of V_label in case we want to execute it for "S" and "W" labels.
    :param theta: Dictionary of vertices and their labels.
    :return: Tuple of sets of vertices from sigma (first set) amd omega (second set).
    """
    return_set_S = set()
    return_set_W = set()
    for vertex, label in theta.items():
        if label.in_sigma:
            return_set_S.add(vertex)
        elif label.in_omega:
            return_set_W.add(vertex)
    return return_set_S, return_set_W
