#!/usr/bin/env python
import pydot

n =10
graph = pydot.Dot('cluster_graph', graph_type='graph') 
subg = pydot.Subgraph('', rank='same') 
graph.add_subgraph(subg)
for i in range(n):
  subg.add_node(pydot.Node("C"+str(i))) 
for i in range(n-1):
  subg.add_edge(pydot.Edge("C"+str(0),"C"+str(i+1))) 
#print graph.create()
graph.write("./graph.png",prog="dot",format="png")
