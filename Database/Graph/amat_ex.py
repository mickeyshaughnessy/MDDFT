#!/usr/bin/env python
import pydot

mat = [[0,2,4,1],
       [1,0,1,1],
       [5,0,0,1],
       [1,0,1,0]]
 
#graph = pydot.Dot('cluster_graph', graph_type='graph') 
#subg = pydot.Subgraph.graph_from_adjaceny_matrix(amat,node_prefix="C")
#graph.add_subgraph(subg)
graph = pydot.graph_from_adjacency_matrix(mat)
graph.write("./graph.ps",prog="dot",format="ps")
