#!/usr/bin/env python3.7
#import sys
#sys.path.insert(1,'/projects/academic/josewalt/caigao/gurobi-env/lib/python3.7/site-packages/')
import gurobipy as gp
from gurobipy import GRB
from itertools import combinations

from sys import argv
import re
import math
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import time
import cProfile
import matplotlib.tri as mtri








def find_all_circles(selected, n, entry, exit):
	visited = [False]*n
	# print(" matrix = ", selected)
	visited[entry] = True
	cur_node = entry
	# path = [cur_node]
	Cirset = []
	circle = [cur_node]
	while True:
		zeronext = True
		i = 0
		while i < n:
			if selected[cur_node, i] > 0.5:
				zeronext = False
			if selected[cur_node, i] > 0.5 and visited[i] == False:
				circle.append(i)
				visited[i] = True
				cur_node = i
				break
			i += 1			
		if i == n: # could be a cycle or a path. 
			# print(circle)
			if zeronext != True and len(circle) > 1:
				Cirset.append(circle)
			j = 0
			while j < n: # start from another unvisited node
				if visited[j] == False:
					cur_node = j
					visited[j] = True
					circle = [cur_node]
					break
				j += 1
			if j == n: # no unvisited node left; terminate
				break 
	# print(Cirset)
	return Cirset     

def cb_subtourelim(model, where):
	if where == GRB.Callback.MIPSOL:
		vals = model.cbGetSolution(model._vars)
		n = 0
		for i, j in model._vars.keys():
			n = max(i,j)
		n += 1
		# print(" n = ", n)
		if False:
			for i in range(0, n):
				for j in range(0, n):
					print(int(vals[i,j]), end= ' ')
				print('')    
		
		# Cirset = find_circle(vals, n, model._entry, model._exit)
		Cirset = find_all_circles(vals, n, model._entry, model._exit)

		if (len(Cirset) != 0):
			for circle in Cirset:
				constr = 0
				for i in range(0,len(circle)-1):
					constr += model._vars[circle[i],circle[i+1]]
				constr += model._vars[circle[-1],circle[0]]
				model.cbLazy(constr <= len(circle)-1)

''' ----------------------------------------- '''
''' Risk & Reward Asset Routing Problem Model '''
''' ----------------------------------------- '''
def solve_CSP(rGgraph, graphsize, nb_tars,):
	n = nb_tars + 2
	try:
		# Create a new model
		model = gp.Model("RRARP")
		model.setParam(GRB.Param.OutputFlag, 1)
		model.setParam(GRB.Param.TimeLimit, 600.0)

		x = model.addVars(n, n, vtype=GRB.BINARY,name="x")
		y = model.addVars(graphsize, graphsize, vtype=GRB.BINARY, name="y")

		constr = 0
		for i in range(0,graphsize):
			for j in range(0,graphsize):
				if rG[i][j] > 0:
					constr += G[i][j] * y[i,j]
		model.setObjective(constr, GRB.MINIMIZE)
		# model.modelSense = GRB.MINIMIZE

		for i in range(0, n):		
			constr1 = 0
			constr2 = 0 
			model.addConstr(x[i,i] == 0, 'dig_'+str(i))
			for j in range(0, n):
				constr1 += x[i,j]
				constr2 += x[j,i]         
			model.addConstr(constr1 == 1, 'row_'+str(i))
			model.addConstr(constr2 == 1, 'col_'+str(i))
			model.update()

		model.addConstr(x[0, nb_tars+1] == 0, 's_t')


		constr = 0 
		for i in range(0, n):
			if G[i][exit] > 0:
				constr += x[i,exit]         
		model.addConstr(constr == 1, 't1') 
		# constr = 0 
		# for i in range(0, n):
		# 	if G[exit][i] > 0 and i != exit:
		# 		constr += x[exit,i]       
		# model.addConstr(constr == 0, 't2') 

		for i in range(0,n):
			if i == entry or i == exit:
				continue
			constr = 0 
			for j in range(0,n):
				if G[i][j] > 0:
					constr += x[i,j] 
			for j in range(0,n):
				if G[j][i] > 0:
					constr -= x[j,i] 
			model.addConstr(constr == 0, "flow_"+str(i))
			model.update()

		for i in range(0,n):
			# if i == entry or i == exit:
			# 	continue
			constr = 0 
			for j in range(0,n):
				if G[i][j] > 0:
					constr += x[i,j]
				if G[j][i] > 0:
					constr += x[j,i] 
			model.addConstr(constr <= 2, "deg_"+str(i))
			model.update()
	 
		for i in range(0, n):
			for j in range(0, n):
				if abs(G[i][j]) < 0.000000001:
					model.addConstr(x[i,j] == 0, "force"+str(i)+str(j))
		model.update()

		if G[entry][exit] > 0:
			model.addConstr(x[entry,exit] == 0, "force"+str(entry)+str(exit))
			model.update()

		constr1 = 0
		for i in range(0,n):
			for j in range(0,n):
				if G[i][j] > 0:
					constr1 += x[i,j] * (obp_rews[i] + obp_rews[j])
		constr1 += obp_rews[entry] + obp_rews[exit]
		model.addConstr(constr1/2.0 >= demand, "rc")

		model.update()

		# model.write('model.lp')

		model._vars = x
		model.Params.lazyConstraints = 1
		model.optimize(cb_subtourelim)
		# model.optimize()

		''' ------------- model output  ----------------------'''
		if model.status == GRB.OPTIMAL or model.status == GRB.TIME_LIMIT:
			# for i in range(n):
			# 	for j in range(n):
			# 		print(int(x[i,j].x), end=' ')
			# 	print('\n', end='')
			visited = [False]*n
			visited[entry] = True
			cur_node = entry
			path = [cur_node]
			while True:
				for i in range(0, n):
					if x[cur_node,i].x > 0.5 and visited[i] == False:
						path.append(i)
						visited[i] = True
						cur_node = i
				if cur_node == exit:
					break
			print("the optimal path: ", path)
			collec_r = 0
			for e in path:
				collec_r += obp_rews[e]
			print('total: ',total_rew, 'collected rewards: ', collec_r, ' requirment:', demand)
			total_risk = 0
			for i in range(len(path)-1):
				total_risk += G[path[i]][path[i+1]]
			print('total risk over optimal path: %g' % total_risk)
			print('Gurobi Time: %g' % model.Runtime)
			print('BestBound: %g' % model.ObjBound)
			print('MIP Gap: %g' % model.MIPGap)


	except gp.GurobiError as e:
		print('Error code ' + str(e.errno) + ': ' + str(e))

	except AttributeError:
		print('Encountered an attribute error')

	# return path


def read_RRARP_instance(instancefile):
	file_ = open(instancefile, 'r')
	line_ = file_.readline()
	list_ = re.split("\t|\n", line_)
	nb_tars_ = int(list_[0])
	nb_edges_ = int(list_[1])
	#nodes
	line_ = file_.readline()
	list_ = re.split(" |\t|\n", line_)    
	x = []
	y = []
	w = []
	G = []
	for e in list_:
		if e != '':
			node_ = re.split(":", e)
			x.append(float(node_[1]))
			y.append(float(node_[2]))
			w.append(float(node_[3]))
	x = np.array(x)
	y = np.array(y)
	w = np.array(w)
	size = len(x)
	G = np.zeros((size, size))

	line_ = file_.readline()
	list_ = re.split(" |\t|\n", line_)    
	for e in list_:
		if e != '':
			edge_ = re.split(":", e)
			s = int(edge_[0])
			t = int(edge_[1])
			d = float(edge_[2])
			G[s,t] = d
			G[t,s] = d
	file_.close()
	return G, x, y, w




if __name__ == "__main__":
	nb_obps = int(argv[1])
	idx = int(argv[2])
	gamma = float(argv[3])
	dir = '/home/cai/Dropbox/Box_Research/Github/RRARP_LinKern/dat/Inst_Pulse/'
	filename = 'tG_'+ str(nb_obps) +'_2_' + str(idx) +'.dat'
	G, x, y, w = read_pccsp_instance(dir+filename)
	_graph = nx.Graph(G)
	rval_lrpath, path_sink = nx.single_source_dijkstra(_graph,0,1)
	print(path_sink)
	LBw = 0
	for e in path_sink:
		LBw += w[e]
	print(LBw)
	total_rew = np.sum(w)
	t0 = time.time()
	path = solve_CSP(G, w, 0, 1, gamma * (total_rew-LBw) + LBw)
	print("Python Time:",time.time() - t0)



