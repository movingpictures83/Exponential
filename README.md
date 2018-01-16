# Exponential
# Language: Python
# Input: CSV (network)
# Output: NOA (central nodes and centrality values)

PluMA plugin that computes Exponential Centrality (Benzi and Klymko, 2013),
a method for finding the important nodes in a network.

This plugin accepts input as a CSV file where both rows and columns
represent nodes in the network, and entry (i, j) represents the weight
of the edge from i to j.

The plugin then produces output as a NOde Attribute (NOA) file for
Cytoscape, containing a ranked list of nodes and centrality values.
This file can then be imported into Cytoscape, giving each node a new
property Centrality that can be used for further downstream analysis,
or visualization.
