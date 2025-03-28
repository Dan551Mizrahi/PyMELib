<h1>PyMELib (Python Minimal Enumeration Library) - A Python Library for Fixed-Parameter-Linear Delay Enumeration of Minimal Dominating Sets and Minimal Hitting Sets</h1>

<p>This repository houses a Python implementation of the algorithm presented in the research paper "Enumeration of Minimal Hitting Sets Parameterized by Treewidth" by Kenig and Mizrahi (2024) [1]. The library aims to provide an efficient way to enumerate minimal hitting sets (and minimal dominating sets) of hypergraphs (of graphs), leveraging the theoretical guarantees outlined in the paper, particularly focusing on fixed-parameter tractability with respect to treewidth. This library also includes an implementation of the novel disjoint branch nice tree decomposition data structure.</p>

<h2>Introduction</h2>

<p>Minimal Hitting Set (MHS) enumeration is a fundamental problem with applications in various domains like databases, AI, bioinformatics, and constraint satisfaction. Given a hypergraph (a collection of sets), the goal is to find all minimal subsets of the hypergraph's vertex set that have a non-empty intersection with every hyperedge (set in the collection). While generally it is not known if the problem has computationally tractable delay-time solution, the algorithm by Kenig and Mizrahi provides an efficient approach for hypergraphs of bounded treewidth, achieving fixed-parameter-linear delay after an FPT preprocessing phase. This library implements that approach.</p>

<h2>Files</h2>

<ul>
    <li>PyMELib/
    <ul>
    <li><code>utils/</code>: Directory containing helper functions for:
        <ul>
            <li><code>readHypergraphFromFile.py</code>: Reading input files and transforming them to the reduction graph.</li>
            <li><code>addConstraints.py</code>: Adding constraints on vertices inclusion/exclusion in the output.</li>
            <li><code>graph_utils.py</code>: General graph utilities.</li>
            <li><code>comb_utils.py</code>: Combinatorial utilities.</li>
            <li><code>labels_utils.py</code>: Utilities related to labels.</li>
        </ul>
    </li>
    <li><code>EnumerationAlgorithms.py</code>: Contains the enumeration algorithms.</li>
    <li><code>TreeDecompositions.py</code>: Handles tree decomposition related tasks.</li>
    <li><code>PreprocessingAlgorithms.py</code>: Contains preprocessing algorithms.</li>
    <li><code>Factors.py</code>: Contains classes used as factors.</li>
    <li><code>labels.py</code>: Handling Labels' definition and logic as defined in [1].</li>
    </ul>
    </li>
    <li><code>tests/</code>: Contains unit tests for various components (e.g., <code>test_basicEnumDSAlgo.py</code>, <code>test_basicEnumHSAlgo.py</code>).</li>
    <li><code>setup.py</code>: Build and installation script.</li>
</ul>

<h2>Functionality</h2>
The following list include the main functionalities of the library, more detailed descriptions can be found in the rest of this README and inside the docstrings of the code:
<ul>
    <li><b>Tree Decompositions:</b> Includes different classes for different types of tree decompositions and functions to work with them. You can look at <code>TreeDecompositions.py</code> which includes:
    <ol>
    <li><code>RootedTreeDecomposition</code>: Class that represents some rooted junction tree of a graph. 
    Its <code>__init__</code> gets a <code>NetworkX</code> graph and construct a rooted tree decomposition of it using <code>NetworkX.junction_tree</code>.
    This class also includes a basic visualisation option of the tree decomposition, <code>self.draw</code>.</li>
    <li><code>RootedNiceTreeDecomposition</code>: Class that represents a nice tree decomposition (it inherits from <code>RootedNiceTreeDecomposition</code>). It has two modes for "niceness": 1. regular "niceness" - when <code>semi_nice=False</code> 
    it transforms the regular rooted tree decomposition to a nice one using the function <code>self.transform_to_nice_rec</code>. 2. semi-nice "niceness" - when <code>semi_nice=True</code>  it transforms the regular rooted tree decomposition to another one using <code>self.transform_to_semi_nice_rec</code>, if the rooted tree decomposition can be
    transformed into a disjoint branch one, it will transform it, if not, the function will transform it to the closest tree to a disjoint one, i.e. the tree won;t have the same vertex in both branches if it doesn't have to be there. Additionally, this class also provide a complete order of the vertices <code>Q</code> as described in [1].</li>
    <li><code>NodeType</code>: This <code>IntEnum</code> represents the types of nodes (bags) in nice tree decomposition, i.e. Forget, Join, ...</li>
    <li><code>RootedDisjointBranchNiceTreeDecomposition</code>: This new data structure introduced in [1] is implemented using this class. Like before it includes a regular option using <code>semi_dntd = False</code> and a "semi" one using <code>semi_dntd = True</code>.</li>
    </ol></li>
    <li><b>Preprocessing:</b> Implements necessary preprocessing algorithm as described in [1], constructing the appropriate factors before enumeration for the RootedDisjointBranchNiceTreeDecomposition. You can look at <code>PreprocessingAlgorithms.py</code> which includes:
    <ol>
    <li><code>create_factors(td: RootedDisjointBranchNiceTreeDecomposition)</code>: Creates empty factors for the nodes of the TD before the preprocessing phase.</li>
    <li><code>calculate_factors_for_mds_enum(td: RootedDisjointBranchNiceTreeDecomposition, current_node: int, options_for_labels=False)</code>: This is a dynamic programming algorithm that calculates the factors of the TD as described in [1].</li>
    <li><code>calculate_factors_for_mds_enum_iterative(td: RootedDisjointBranchNiceTreeDecomposition, options_for_labels=False)</code>: This is a for loop version of <code>calculate_factors_for_mds_enum</code>, using a stack. Much better in terms of both memory and time consumption.</li>
    </ol>
    </li>
    <li><b>Core Enumeration Algorithms:</b> The main implementation of the minimal dominating set enumeration with fixed-parameter-linear delay. You can look at <code>EnumerationAlgorithms.py</code> which includes:
    <ol>
    <li><code>IsExtendable(td: RootedDisjointBranchNiceTreeDecomposition, theta, i)</code>: IsExtendable procedure described in [1].</li>
    <li><code>IncrementLabeling(td: RootedDisjointBranchNiceTreeDecomposition, theta: Dict[str, Label], i, c: Label, V_label_S, V_label_W)</code>: An efficient version of IncrementLabeling procedure described in [1].</li>
    <li><code>EnumMDS(td: RootedDisjointBranchNiceTreeDecomposition, theta: Dict[str, Label] = dict(), i=0, debug_flag=False, options_for_labels=False)</code>: The algorithm for minimal dominating sets enumeration described in [1]. Theoretically the time delay is bounded by $O(nw)$, $n$ being the number of nodes in the graph and $w$ is the graph's treewidth.</li>
    <li><code>EnumMDS_iterative(td: RootedDisjointBranchNiceTreeDecomposition, debug_flag=False, options_for_labels=False)</code>: This is a for loop version of <code>EnumMDS</code>, using a stack.</li>
    <li><code>EnumMHS(td: RootedDisjointBranchNiceTreeDecomposition, theta: Dict[str, Label] = dict(), i=0, debug_flag=False)</code>: This is a version of <code>EnumMDS</code> designated for minimal hitting set enumeration.</li>
    <li><code>EnumMHS_iterative(td: RootedDisjointBranchNiceTreeDecomposition, debug_flag=False)</code>: This is a for loop version of <code>EnumMHS</code>, using a stack.</li>
    </ol></li>
    <li><b>Input:</b> Utilities for reading hypergraph data in the input format of the dualization repository [2], constructing the proper reduction to enumeration of minimal dominating sets described in [1] and adding other inclusion/exclusion constraints on the problem. You can look at <code>utils/readHypergraphFromFile.py</code> and <code>utils/addConstraints.py</code> which include:
    <ol>
    <li><code>add_constraints_on_graph(G: nx.Graph, include_in_ds: Iterable = [], exclude_from_ds: Iterable = [])</code>: This function is a helper function to add constraints on a graph, before the running of the preprocessing phase or the enumeration phase. There are two kinds of constraints: include_in_ds and exclude_from_ds. After using this function, remember to call the preprocessing phase and the enumeration phase with options_for_labels=True.</li>
    <li><code>read_hypergraph(path_to_file: str) -> nx.Graph</code>: Reads a hypergraph from a file in the designated format (see Data - Input Format) and transform it to the reduction graph described in [1].</li>
    </ol></li>
</ul>

<h2>Getting Started</h2>

<h3>Prerequisites</h3>
<ul>
    <li>Python 3.9.6 or higher (as specified by you)</li>
    <li>Dependencies:
        <ul>
            <li>NetworkX</li>
            <li>Matplotlib</li>
            <li>EoN (Epidemics on Networks) [3]</li>
            <li>Plotly</li>
            <li>tqdm</li>
        </ul>
    </li>
</ul>

<h3>Installation</h3>
Using pip: <code>pip install PyMELib</code>
<p>or clone the repository and install it locally:</p>
<ol>
    <li>Clone the repository: <br> <code>git clone https://github.com/Dan551Mizrahi/PyMELib.git</code><br><code>cd PyMELib</code></li>
    <li>Install the package (preferably in a virtual environment):<br> <code>pip install .</code><br> or for development:<br> <code>pip install -e .</code></li>
</ol>

<h3>Running</h3>
<p>Import the library contents:</p>
<pre>
<code class="language-python">
from PyMELib.TreeDecompositions import RootedDisjointBranchNiceTreeDecomposition<br>
from PyMELib.PreprocessingAlgorithms import create_factors, calculate_factors_for_mds_enum_iterative<br>
from PyMELib.EnumerationAlgorithms import EnumMHS<br>
from PyMELib.utils.readHypergraphFromFile import read_hypergraph
from PyMELib.utils.addConstraints import add_constraints_on_graph
</code>
</pre>
<p>Read your hypergraph:</p>
<pre>
<code class="language-python">
H = read_hypergraph("path/to/your/hypergraph.hg") # Replace with the actual path to your hypergraph file
</code>
</pre>
<p>Add constraints:</p>
<pre>
<code class="language-python">
include_in_ds = [1, 2] # Example vertices to include
exclude_from_ds = [3, 4] # Example vertices to exclude
add_constraints_on_graph(H, include_in_ds, exclude_from_ds)
</code>
</pre>
<p>Build tree decomposition:</p>
<pre>
<code class="language-python">
td = RootedDisjointBranchNiceTreeDecomposition(H, semi_dntd=True)
</code>
</pre>
<p>Create factors and run the preprocessing phase:</p>
<pre>
<code class="language-python">
create_factors(td)
calculate_factors_for_mds_enum_iterative(td, options_for_labels=True)
</code>
</pre>
<p>Run the enumeration algorithm:</p>
<pre>
<code class="language-python">
results = EnumMHS_iterative(td)
</code>
</pre>
<p><code>results</code> is now a Python generator capable of generating all the minimal hitting sets of the hypergraph.</p>

<h2>Examples</h2>

<pre>
<code class="language-python">
from PyMELib.TreeDecompositions import RootedDisjointBranchNiceTreeDecomposition
from PyMELib.PreprocessingAlgorithms import create_factors, calculate_factors_for_mds_enum_iterative
from PyMELib.EnumerationAlgorithms import EnumMHS
from PyMELib.utils.readHypergraphFromFile import read_hypergraph
from PyMELib.utils.addConstraints import add_constraints_on_graph

# 1. Load hypergraph
# Replace with the actual path to your hypergraph file
H = read_hypergraph("path/to/your/hypergraph.hg") 

# 2. create tree decomposition
td = RootedDisjointBranchNiceTreeDecomposition(H, semi_dntd=True)

# 3. Preprocess hypergraph
create_factors(td)
calculate_factors_for_mds_enum_iterative(td, options_for_labels=True)

# 4. Run enumeration
results = EnumMHS_iterative(td)

# 5. Print results
for mhs in results:
    print(mhs)
</code>
</pre>

<h2>Data - Input Format</h2>
<p>If you want to run the code in order to enumerate minimal dominating sets of a regular graph you can simply pass to the constructor of the tree decomposition a <code>NetworkX</code> graph object. <br>
In the case of hypergraphs, we are consistent with the input format of the dualization repository [2]. The input file should be in the following format, each line representing a hyperedge (set) of the hypergraph, vertices are non-negative integers and are separated by spaces:</p>
<pre><code>vertex1 vertex2 vertex3
vertex2 vertex4
...
</code></pre>

<h2>Results - Output Format</h2>

<p>The minimal hitting sets/dominating sets return as Python frozensets:</p>
<pre>
<code class="language-python">
# Example output
frozenset({1, 2, 3})
frozenset({2, 4})
frozenset({1, 4})
...
</code>
</pre>

<h2>Limitations</h2>

<ul>
<li>The efficiency relies on the input hypergraph having a small treewidth. The preprocessing phase complexity depends exponentially on the treewidth.</li>
</ul>

<h2>Authors</h2>
<ul>
<li>Dan S. Mizrahi</li>
<li>Batya Kenig</li>
</ul>

<h2>AI usage</h2>
This code was created with the help of GitHub Copilot and Gemini.

<h2>Further Reading</h2>

[1] Kenig, Batya, and Dan Shlomo Mizrahi. "Enumeration of Minimal Hitting Sets Parameterized by Treewidth." arXiv preprint arXiv:2408.15776 (2024).

[2] Keisuke Murakami & Takeaki Uno (uno@nii.jp). Hypergraph Dualization Repository - Program Codes and Instances for Hypergraph Dualization (minimal hitting set enumeration). <a href="https://research.nii.ac.jp/~uno/dualization.html">https://research.nii.ac.jp/~uno/dualization.html</a>.

[3] Miller et al., (2019). EoN (Epidemics on Networks): a fast, flexible Python package for simulation, analytic approximation, and analysis of epidemics on networks. Journal of Open Source Software, 4(44), 1731, https://doi.org/10.21105/joss.01731