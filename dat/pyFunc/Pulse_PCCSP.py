from sys import argv
import re
import math
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import time
from find_admispaths import solve_CSP

class Pulse():
	def __init__(self, graphx, budget_pct, minrisk_to_sink, OPrewards):
		self._graph = graphx
		self._nb_nodes = len(self._graph)
		self._source = 0
		self._sink = len(self._graph.nodes) -1
		self._total_reward = sum(OPrewards)
		self._budget = budget_pct * sum(OPrewards)
		# best known objective value
		self._curbest_objval = float("inf") 
		self._curbest_path = []
		self._leastrisk_to_sink = minrisk_to_sink
		self._OPrewards = OPrewards
		self._buckets = []
		for i in range(len(graphx.nodes)):
			self._buckets.append([])
		self._iterations = 0
		self._time1 = 0
		self._time2 = 0
		# self.print_inst_info()


	def print_inst_info(self):
		print("source node: ", self._source)
		print("sink node: ", self._sink)
		print("total reward: ", self._total_reward)
		print("required budget: ", self._budget)
	def print_optimal_sol(self):
		self.print_inst_info()
		print("optimal legnth: ", self._curbest_objval, '\t optimal path:', self._curbest_path)
		print("time: ", self._time1, self._time2)

	def recursive_search(self, curnode, path, pathreward, pathrisk, log_on):
		if log_on:
			print("\n")
			print(self._iterations, ". consider add node", curnode, "to current path:", path, " [cur path reward:", pathreward, "cur path risk:", pathrisk,"]")

		self._iterations += 1
		# if self._iterations > 20:
		# 	return None

		if curnode == self._sink:
			if log_on:
				print("ARRIVE AT SINK")		
			if pathreward >= self._budget and pathrisk + self._graph.edges[path[-1],self._sink]['weight'] < self._curbest_objval:
				self._curbest_path = path.append(self._sink)
				self._curbest_objval = pathrisk + self._graph.edges[path[-1],self._sink]['weight'] 
				if log_on:
					print("============>>>>!!!!find a better path at SINK !!! total risk:", self._curbest_objval)
					print("============>>>>!!!!find a better path at SINK !!! new opt path:", self._curbest_path)
			else:
				if log_on:
					print("  arrive at the sink but path is no better ")
			return None

		if curnode != self._source:
			lastnode = path[-1]
			if curnode not in path: # no cycle 
				if log_on:
					print("No CYCLE")
					print("plus the risk of going to curnode:", curnode, "new cumulative risk: ", pathrisk + self._graph.edges[lastnode,curnode]['weight'])					
				if pathrisk + self._graph.edges[lastnode,curnode]['weight'] < self._curbest_objval:
					if log_on:
						print("PARTIAL PATH: CHECK ITS BUDGET")
					if pathreward + self._OPrewards[curnode] >= self._budget:
						if log_on:
							print("SATISFY BUDGET")
						comp_set = [i for i in range(self._nb_nodes) if i not in path] 
						subgraph = self._graph.subgraph(comp_set)
						t0 = time.clock()
						if nx.is_connected(subgraph):
							self._time1 += time.clock() - t0;
							t0 = time.clock()
							rVal_sink, path_sink = nx.single_source_dijkstra(subgraph,curnode,self._sink)
							self._time2 += time.clock() - t0;

							# rVal_sink = nx.dijkstra_path_length(subgraph,curnode,self._sink)
							''' update the curbest obj value'''
							if pathrisk + self._graph.edges[lastnode,curnode]['weight'] + rVal_sink < self._curbest_objval:
								self._curbest_objval =pathrisk + self._graph.edges[lastnode,curnode]['weight'] + rVal_sink
								self._curbest_path = path + path_sink
								if log_on:
									print("============>>>>!!!!find a better path !!! total risk:", self._curbest_objval)
									print("============>>>>!!!!find a better path !!!"  "new opt path:", self._curbest_path)
							else:
								if log_on:
									print(" new LONGER feasible path is found. ")
						else:
							'''subgraph is disconnected '''
							if log_on:
								print("subgraph is disconnected")

					else:
						if log_on:
							print("NOT SATISFY BUDGET")

							print("ARRIVE AT DOMINANCE CHECK: all its neighbors", list(self._graph.neighbors(curnode)))
						for ngb in self._graph.neighbors(curnode):
							if log_on:							
								print('---->', end=' ')
							if ngb not in path:
								if log_on:
									print(ngb,'is not visited', end='\t')
								newpath = path+[curnode]
								if log_on:
									print("~~~=== ", path, newpath)
								newpathreward = pathreward + self._OPrewards[curnode]
								newpathrisk = pathrisk + self._graph.edges[lastnode,curnode]['weight']
								if log_on:
									print("***new path: ", newpath, 'newpathreward:', newpathreward, 'newpathrisk:', newpathrisk)
								self.recursive_search(ngb, newpath, newpathreward, newpathrisk, log_on)
							else:
								if log_on:
									print(ngb,'is already visited')
				else:
					if log_on:
						print("PARTIAL PATH IS PRUNED BY PRIMAL BOUND")
		else:
			for ngb in self._graph.neighbors(curnode):
				if log_on:
					print("*** Propogate from source to node ", ngb)
				self.recursive_search(ngb, [self._source], 0, 0, log_on)
			


if __name__ == "__main__":
	riskmatx = np.array([[0, 2, 3, 0, 0, 0, 0], 
						 [2, 0, 9, 1, 2, 0, 0], 
						 [3, 9, 0, 1, 0, 1, 0], 
						 [0, 1, 1, 0, 1, 1, 0], 
						 [0, 2, 0, 1, 0, 5, 2],
						 [0, 0, 1, 1, 5, 0, 1],
						 [0, 0, 0, 0, 2, 1, 0]])
	OBPREWtest = [0, 2, 3, 3, 4, 2, 0]

 
	t0 = time.clock()
	G = nx.Graph(riskmatx)
	# spaths = nx.shortest_path(G, source=len(G.nodes)-1)
	# print(spaths)
	# _minrisk = [0]*len(G.nodes)
	# for i in range(len(G.nodes)):
	# 	print(spaths[i])
	# 	if len(spaths[i]) >= 2:
	   #  	for j in range(len(spaths[i])-1):
	   #  		_minrisk[i] += riskmatx[spaths[i][j]][spaths[i][j+1]]
	minrisk_sink = nx.single_source_dijkstra_path_length(G,source=len(G.nodes)-1)
	budget_pct = 1.0
	print('--------------------------------------------------------')
	path_model = solve_CSP(riskmatx, OBPREWtest, budget_pct)
	print('--------------------------------------------------------')
	test = Pulse(G, budget_pct, minrisk_sink, OBPREWtest)
	# print(test._leastrisk_to_sink )
	test.recursive_search(0, [], 0.0, 0.0, False)
	test.print_optimal_sol()
	print("Running Time: ", time.clock() - t0, "seconds ")


	