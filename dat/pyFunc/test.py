import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

# Fixing random state for reproducibility
# np.random.seed(19680801)

# # Compute areas and colors
# N = 150
# r = 2 * np.random.rand(N)

# theta = 2 * np.pi * np.random.rand(N)
# area = 200 * r**2
# colors = r
# area = np.array([np.pi*14.0]*N)

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='polar')
# c = ax.scatter(theta, r, c=colors, s=area, cmap='hsv', alpha=0.75)
# plt.show()

riskmatx = np.array([[0, 2, 3, 0, 0, 0, 0], 
 		             [2, 0, 9, 1, 2, 0, 0], 
 		             [3, 9, 0, 1, 0, 1, 0], 
 		             [0, 1, 1, 0, 1, 1, 0], 
 		             [0, 2, 0, 1, 0, 5, 2],
 		             [0, 0, 1, 1, 5, 0, 1],
 		             [0, 0, 0, 0, 2, 1, 0]])
   
G = nx.Graph(riskmatx)
#print(nx.single_source_shortest_path_length(G,6))  # return a dict
S1 = [1,2,3,4,5,6]
S2 = [2,4,5]
comp_set = list(set(S1)-set(S2))
print(comp_set)
# subG = G.subgraph(comp_set)
# print(subG.edges)
# for e in G.neighbors(1):
# 	print(e)
# print(list(G.neighbors(0)))
# print(G.nodes)
# print(G.edges[0,1]['weight'])


# path = [1, 6]
# E = list(G.edges(path))
# e = E[0]
# G[e[0]][e[1]]['weight'] =10
# print(G[e[0]][e[1]]['weight'])



====>>>> Algorithm iniialization: least risky path to sink 
... Propogate from source 3 to node 34
... Propogate from source 3 to node 64
... Propogate from source 3 to node 72
... Propogate from source 3 to node 74
 =====>>>> Time: 11.5711
 =====>>> total reward: 8.712, Demand: 3.4848, Reward on min-risk path: 0.837
 =====>>> optimal obj value: 3.12174, collected reward: 3.52
          optimal path [3 74 34 47 24 72 21 52 62 41 45 49 67 16 28 43 54 31 61 37 17 68 65 63 32 15 ]
