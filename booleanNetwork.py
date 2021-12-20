import itertools
import numpy as np
import matplotlib
import networkx as nx  # because we want it for Graph generation
import matplotlib.pyplot as plt


# function to generate a list of interaction between the genes
def generate_edges(graph_edges):
    edge = []
    for node in graph_edges:
        for neighbor in graph_edges[node]:
            edge.append((node, neighbor))
    return edge


# function for state upgrading with boolean rules
def bool_rule(ini_state):
    PR1 = ini_state[0]
    PDF12 = ini_state[1]
    WRKY70 = ini_state[2]
    WRKY54 = ini_state[3]
    WRKY33 = ini_state[4]
    MYC2 = ini_state[5]
    ERF1 = ini_state[6]
    LOX2 = ini_state[7]

    PR1_u = int((WRKY70 or ERF1 or PR1) and (not WRKY54))
    PDF12_u = int((WRKY54 or WRKY33 or MYC2 or ERF1) and (not LOX2 or not WRKY70))
    WRKY70_u = int(ERF1 and (not MYC2 or not WRKY54))
    WRKY54_u = int((WRKY33 or ERF1) and (not MYC2 or not WRKY54))
    WRKY33_u = int(PR1 and ((not PDF12) or (not WRKY33) or (not LOX2)))
    MYC2_u = int(not ERF1)
    ERF1_u = int((not WRKY33) or (not WRKY54))
    LOX2_u = int((not ERF1) or (not PDF12) or (not LOX2))

    state_n = [PR1_u, PDF12_u, WRKY70_u, WRKY54_u, WRKY33_u, MYC2_u, ERF1_u, LOX2_u]
    return state_n




# create a dictionarie with different node and their connections
connection = {'PR1': {'WRKY70', 'ERF1', 'PR1', 'WRKY54'},
              'PDF12': {'WRKY54', 'WRKY33', 'MYC2', 'ERF1', 'LOX2', 'WRKY70'},
              'WRKY70': {'ERF1', 'MYC2', 'WRKY54'}, 'WRKY54': {'WRKY33', 'ERF1', 'MYC2', 'WRKY54'},
              'WRKY33': {'PR1', 'PDF12', 'WRKY33', 'LOX2'}, 'MYC2': {'ERF1'}, 'ERF1': {'WRKY33', 'WRKY54'},
              'LOX2': {'ERF1', 'PDF12', 'LOX2'}}

# create a dictionary that contain only negativ regulation connection
red_connection = {'PR1': {'WRKY54'}, 'PDF12': {'LOX2', 'WRKY70'}, 'WRKY70': {'MYC2', 'WRKY54'},
                  'WRKY54': {'MYC2', 'WRKY54'},
                  'WRKY33': {'PDF12', 'WRKY33', 'LOX2'}, 'MYC2': {'ERF1'}, 'ERF1': {'WRKY33', 'WRKY54'},
                  'LOX2': {'ERF1', 'PDF12', 'LOX2'}}
#  create Node list
listofnode = ['PR1', 'PDF12', 'WRKY70', 'WRKY54', 'WRKY33', 'MYC2', 'ERF1', 'LOX2']

# create a list of all the edges
edges = generate_edges(connection)
print(edges)

# recuperate only negativ regulation
red_edges = generate_edges(red_connection)
print("la liste des regulations negative est ", red_edges)

# recuperate only positive regulation interaction
pos_edge = []
for e in edges:
    if e not in red_edges:
        pos_edge.append(e)
print(pos_edge)

# Draw regulation graph
gene_interact = nx.MultiDiGraph()  
gene_interact.add_nodes_from(listofnode)  
gene_interact.add_edges_from(pos_edge) 
gene_interact.add_edges_from(red_edges)  

 # draw Nodes
nx.draw(gene_interact.reverse(), with_labels=True, node_size=1000, node_color="skyblue", nodelist=listofnode,
        pos=nx.shell_layout(gene_interact)) 

 # draw pos edge colored in black
nx.draw_networkx_edges(gene_interact.reverse(), pos=nx.shell_layout(gene_interact), edgelist=pos_edge, alpha=0.5,
                       width=2, edge_color="black") 

# draw negative edge cololored in black
nx.draw_networkx_edges(gene_interact.reverse(), pos=nx.shell_layout(gene_interact), edgelist=red_edges, alpha=0.5,
                       width=2, edge_color="red")  

plt.show()


states = list(map(list, itertools.product([0, 1], repeat=8)))
print(states)
states_n = list(map(bool_rule, states))
print(states_n)



