import networkx as nx
from PyMELib.PreprocessingAlgorithms import *
from PyMELib.EnumerationAlgorithms import *

example_graph1 = nx.Graph()
example_graph1.add_edges_from([(1, 2),
                             (2, 3),])
example_mds1 = {frozenset({2}), frozenset({1, 3})}

example_graph2 = nx.Graph()
example_graph2.add_edges_from([(97, 98),
                             (97, 99),
                             (97, 102),
                             (97, 103),
                             (104, 102),
                             (102, 103),
                             (103, 105),
                             (98, 100),
                             (99, 100),
                             (100, 101),
                             ])
example_mds2 = {frozenset({100, 102, 103}),
                frozenset({98, 99, 101, 103, 104}),
                frozenset({98, 99, 101, 104, 105}),
                frozenset({98, 99, 101, 102, 105}),
                frozenset({98, 99, 101, 102, 103}),
                frozenset({104, 100, 103}),
                frozenset({105, 100, 102}),
                frozenset({104, 105, 99, 100}),
                frozenset({104, 105, 98, 100}),
                frozenset({104, 97, 101, 105}),
                frozenset({104, 97, 101, 103}),
                frozenset({97, 101, 102, 105}),
                frozenset({97, 101, 102, 103}),
                frozenset({104, 97, 100, 105}),
                }

example_graph3 = nx.Graph()
example_graph3.add_edges_from([(97, 98),
                             (97, 99),
                             (97, 102),
                             (97, 103),
                             (104, 102),
                             (102, 103),
                             (103, 105),
                             (98, 100),
                             (99, 100),
                             (100, 101),
                             (106, 97),
                             (106, 98),
                             (106, 99),
                             (106, 100),
                             (106, 101),
                             (106, 102),
                             (106, 103),
                             (106, 104),
                             (106, 105),
                             ])
example_mds3 = {frozenset({106}),
                frozenset({100, 102, 103}),
                frozenset({98, 99, 101, 103, 104}),
                frozenset({98, 99, 101, 104, 105}),
                frozenset({98, 99, 101, 102, 105}),
                frozenset({98, 99, 101, 102, 103}),
                frozenset({104, 100, 103}),
                frozenset({105, 100, 102}),
                frozenset({104, 105, 99, 100}),
                frozenset({104, 105, 98, 100}),
                frozenset({104, 97, 101, 105}),
                frozenset({104, 97, 101, 103}),
                frozenset({97, 101, 102, 105}),
                frozenset({97, 101, 102, 103}),
                frozenset({104, 97, 100, 105}),
                }

def test_EnumMDS_recursive():
    td = RootedDisjointBranchNiceTreeDecomposition(example_graph1, debug_flag=False, semi_dntd=False)
    create_factors(td)
    calculate_factors_for_mds_enum(td, td.get_root())
    assert set(EnumMDS(td, debug_flag=False)) == example_mds1

    td = RootedDisjointBranchNiceTreeDecomposition(example_graph2, debug_flag=False, semi_dntd=False)
    create_factors(td)
    calculate_factors_for_mds_enum(td, td.get_root())
    assert set(EnumMDS(td, debug_flag=False)) == example_mds2

    td = RootedDisjointBranchNiceTreeDecomposition(example_graph3, debug_flag=False, semi_dntd=False)
    create_factors(td)
    calculate_factors_for_mds_enum(td, td.get_root())
    assert set(EnumMDS(td, debug_flag=False)) == example_mds3

def test_EnumMDS_iterative():
    td = RootedDisjointBranchNiceTreeDecomposition(example_graph1, debug_flag=False, semi_dntd=False)
    create_factors(td)
    calculate_factors_for_mds_enum_iterative(td)
    assert set(EnumMDS_iterative(td, debug_flag=False)) == example_mds1

    td = RootedDisjointBranchNiceTreeDecomposition(example_graph2, debug_flag=False, semi_dntd=False)
    create_factors(td)
    calculate_factors_for_mds_enum_iterative(td)
    assert set(EnumMDS_iterative(td, debug_flag=False)) == example_mds2

    td = RootedDisjointBranchNiceTreeDecomposition(example_graph3, debug_flag=False, semi_dntd=False)
    create_factors(td)
    calculate_factors_for_mds_enum_iterative(td)
    assert set(EnumMDS_iterative(td, debug_flag=False)) == example_mds3

def test_EnumMDS_iterative_semi():
    td = RootedDisjointBranchNiceTreeDecomposition(example_graph1, debug_flag=False, semi_dntd=True)
    create_factors(td)
    calculate_factors_for_mds_enum_iterative(td)
    assert set(EnumMDS_iterative(td, debug_flag=False)) == example_mds1

    td = RootedDisjointBranchNiceTreeDecomposition(example_graph2, debug_flag=False, semi_dntd=True)
    create_factors(td)
    calculate_factors_for_mds_enum_iterative(td)
    assert set(EnumMDS_iterative(td, debug_flag=False)) == example_mds2

    td = RootedDisjointBranchNiceTreeDecomposition(example_graph3, debug_flag=False, semi_dntd=True)
    create_factors(td)
    calculate_factors_for_mds_enum_iterative(td)
    assert set(EnumMDS_iterative(td, debug_flag=False)) == example_mds3