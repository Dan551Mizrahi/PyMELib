from PyMELib.PreprocessingAlgorithms import *
from PyMELib.EnumerationAlgorithms import *
from PyMELib.utils.readHypergraphFromFile import read_hypergraph
from PyMELib.utils.addConstraints import add_constraints_on_graph
import os

script_dir = os.path.dirname(os.path.abspath(__file__))

HG1_path = os.path.join(script_dir, "exampleHG.dat")
HG2_path = os.path.join(script_dir, "exampleHG2.dat")

example_HG1 = read_hypergraph(HG1_path)
add_constraints_on_graph(example_HG1, include_in_ds=[5], exclude_from_ds=[6])
example_MHS1 = {
                frozenset({2, 5, 7}),
                frozenset({1, 5, 7}),
                }

example_HG2 = read_hypergraph(HG2_path)
add_constraints_on_graph(example_HG2, include_in_ds=[1], exclude_from_ds=[4])
example_MHS2 = {
        frozenset({1, 3}),
        frozenset({1, 2}),
    }

def test_EnumMDS_recursive():
    td = RootedDisjointBranchNiceTreeDecomposition(example_HG2, debug_flag=False, semi_dntd=False)
    create_factors(td)
    calculate_factors_for_mds_enum(td, td.get_root(), options_for_labels=True)
    assert set(EnumMHS(td, debug_flag=False)) == example_MHS2

def test_EnumMDS_iterative():
    td = RootedDisjointBranchNiceTreeDecomposition(example_HG1, debug_flag=False, semi_dntd=False)
    create_factors(td)
    calculate_factors_for_mds_enum_iterative(td, options_for_labels=True)
    assert set(EnumMHS_iterative(td, debug_flag=False)) == example_MHS1

    td = RootedDisjointBranchNiceTreeDecomposition(example_HG2, debug_flag=False, semi_dntd=False)
    create_factors(td)
    calculate_factors_for_mds_enum_iterative(td, options_for_labels=True)
    assert set(EnumMHS_iterative(td, debug_flag=False)) == example_MHS2


def test_EnumMDS_iterative_semi():
    td = RootedDisjointBranchNiceTreeDecomposition(example_HG1, debug_flag=False, semi_dntd=True)
    create_factors(td)
    calculate_factors_for_mds_enum_iterative(td, options_for_labels=True)
    assert set(EnumMHS_iterative(td, debug_flag=False)) == example_MHS1

    td = RootedDisjointBranchNiceTreeDecomposition(example_HG2, debug_flag=False, semi_dntd=True)
    create_factors(td)
    calculate_factors_for_mds_enum_iterative(td, options_for_labels=True)
    assert set(EnumMHS_iterative(td, debug_flag=False)) == example_MHS2

