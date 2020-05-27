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

def find_circle(selected, n, entry, exit):
	visited = [False]*n
	# print(" matrix = ", selected)
	visited[entry] = True
	cur_node = entry
	''' find the feasible path from entry to exit '''
	path = [cur_node]
	while True:
		for i in range(0, n):
			if selected[cur_node, i] > 0.5 and visited[i] == False:
				path.append(i)
				visited[i] = True
				cur_node = i
		if cur_node == exit:
			break
	# print("path = ", path)
	''' check existence of some circle(s) '''
	flag_circle = False
	start_node = -1
	for i in range(1, n-1):
		if visited[i] == False:            
			for j in range(1, n-1):
				if selected[i,j] > 0.5 and visited[i] == False:
					start_node = i
					flag_circle = True
					break
			if flag_circle:
				break
	circle = []
	# print('the starting node of a circle is ', start_node)
	if flag_circle == True:
		cur_node = start_node
		circle.append(cur_node)
		while True:
			# print('current node ', cur_node)
			for i in range(1, n-1):
				if selected[cur_node, i] > 0.5 and visited[i] == False:
					# print('next visit node ', i)
					visited[i] = True
					cur_node = i
					break
			# print('cur_node = ', cur_node, 'circle = ', circle, ' || ', cur_node in circle)
			if cur_node in circle:
				break
			else:
				circle.append(cur_node)
	# print("circle = ", circle)
	return circle       


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
''' solve Budgeted Prize-Collecting Constrained Shortest Path (PCCSP) model '''
''' ----------------------------------------- '''
def solve_CSP(G, obp_rews, entry, exit, demand, rscG, budget):
	n = G.shape[0]
	total_rew = np.sum(obp_rews)
	if obp_rews[entry] + obp_rews[exit] >= demand:
		return [entry, exit]
	try:
		# Create a new model
		model = gp.Model("PCCSP")
		model.setParam(GRB.Param.OutputFlag, 1)
		model.setParam(GRB.Param.TimeLimit, 100.0)
		model._entry = entry
		model._exit = exit
		x = model.addVars(n, n, vtype=GRB.BINARY,name="x")
		constr = 0
		for i in range(0,n):
			for j in range(0,n):
				if G[i][j] > 0:
					constr += G[i][j] * x[i,j]
		model.setObjective(constr, GRB.MINIMIZE)
		# model.modelSense = GRB.MINIMIZE

		# source node
		constr = 0 
		for i in range(0, n):		
			if G[entry][i] > 0:
				constr += x[entry,i]         
		model.addConstr(constr == 1, 's1') 
		# constr = 0 
		# for i in range(0, n):
		# 	if G[i][entry] > 0 and i != entry:
		# 		constr += x[i,entry]         
		# model.addConstr(constr == 0, 's2') 
		# sink node
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


		constr1 = 0
		for i in range(0,n):
			for j in range(0,n):
				if G[i][j] > 0:
					constr1 += x[i,j] * rscG[i][j]
		model.addConstr(constr1 <= budget, "rsc")
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
			total_used_rsc = 0
			total_risk = 0
			for i in range(len(path)-1):
				total_risk += G[path[i]][path[i+1]]
				total_used_rsc += rscG[path[i]][path[i+1]]
			print('Budget: ', budget, 'Used resource: ', total_used_rsc)

			print('total risk over optimal path: %g' % total_risk)
			print('Gurobi Time: %g' % model.Runtime)
			print('BestBound: %g' % model.ObjBound)
			print('MIP Gap: %g' % model.MIPGap)


	except gp.GurobiError as e:
		print('Error code ' + str(e.errno) + ': ' + str(e))

	except AttributeError:
		print('Encountered an attribute error')

	return path


def read_budgeted_pccsp_instance(instancefile):
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
	rscG = np.zeros((size, size))


	line_ = file_.readline()
	list_ = re.split(" |\t|\n", line_)    
	for e in list_:
		if e != '':
			edge_ = re.split(":", e)
			s = int(edge_[0])
			t = int(edge_[1])
			d = float(edge_[2])
			rsc = float(edge_[3])
			G[s,t] = d
			G[t,s] = d
			rscG[s,t] = rsc
			rscG[t,s] = rsc
	file_.close()
	# plot_ObP_graph([1,1,1], x, y, w)
	return G, rscG, x, y, w

def plot_target_area(center, obp_x, obp_y, obp_reward, path=None):
	fig = plt.figure()
	ax = fig.add_subplot(111)
	plt.axis('off')
	ax.set_xlim([center[0]-center[2]-0.1, center[0]+center[2]+0.1])
	ax.set_ylim([center[1]-center[2]-0.1, center[1]+center[2]+0.1])
	ax.set_aspect('equal', adjustable='box')
	# big circle 
	circle = plt.Circle((center[0], center[1]), radius=center[2], edgecolor='k',fill=False, alpha=0.5)
	ax.add_artist(circle)

	#observation points
	colors = obp_reward
	area = np.array([np.pi*20.0]*len(obp_x))
	c = ax.scatter(obp_x, obp_y, c=colors, s=area, cmap='hsv', alpha=0.85)
	for i, txt in enumerate(range(1,len(obp_x)+1)):
		# ax.annotate(txt, (obp_x[i]+0.02, obp_y[i]+0.02),size=10)
		ax.annotate(txt,  (obp_x[i]+0.02, obp_y[i]-0.02), horizontalalignment='right', verticalalignment='top',size=10)
	if path != None:
		# minimum risk path with profit
		pathX = [entryP[0]]
		# pathX = []
		for e in path[1:len(path)-1]:
			pathX.append(obp_x[e-1])
		pathX.append(exitP[0])

		pathY = [entryP[1]]
		# pathY = []
		for e in path[1:len(path)-1]:
			pathY.append(obp_y[e-1])
		pathY.append(exitP[1])
		plt.plot(pathX, pathY,'k')
	# entry and exit point
	# circle = plt.Circle((entryP[0], entryP[1]), radius=0.02, edgecolor='k',facecolor='k')
	# ax.add_artist(circle)
	# circle = plt.Circle((exitP[0], exitP[1]), radius=0.02, edgecolor='k',facecolor='k')
	# ax.add_artist(circle)

	plt.show()

def plot_ObP_graph(center, obpx, obpy, obprew):
	x = obpx
	y = obpy
	fig = plt.figure()
	ax = fig.add_subplot(111)
	# plt.axis('off')
	ax.set_xlim([(center[0]-center[2]-0.1), (center[0]+center[2]+0.1)])
	ax.set_ylim([(center[1]-center[2]-0.1), (center[1]+center[2]+0.1)])
	ax.set_aspect('equal', adjustable='box')
	''' big circle ''' 
	circle = plt.Circle((center[0], center[1]), radius=center[2], edgecolor='k',fill=False)
	ax.add_artist(circle)
	''' plog points, including boundary points and obser. points, with triangulation '''
	triang = mtri.Triangulation(x, y)
	ax.triplot(triang, c="#D3D3D3", marker='.', markerfacecolor="black", markeredgecolor="black", markersize=5)
	''' plot 16 boundary points '''
	if True:
		for i in range(len(obpx)):
			circle = plt.Circle((x[i], y[i]), radius=0.01, edgecolor='k',facecolor='k')
			ax.add_artist(circle)

	for i in range(0,len(x)):
		ax.annotate(str(i), (x[i]+0.02, y[i]-0.02), horizontalalignment='right', verticalalignment='top',size=10)
	plt.show()



if __name__ == "__main__":
	nb_obps = int(argv[1])
	idx = int(argv[2])
	gamma = float(argv[3])
	dir = '/home/cai/Dropbox/Box_Research/Github/RRARP_LinKern/dat/Inst_Pulse/'
	filename = 'tG2_'+ str(nb_obps) +'_16_' + str(idx) +'.dat'
	G, rscG, x, y, w = read_budgeted_pccsp_instance(dir+filename)
	_graph = nx.Graph(G)
	rval_lrpath, path_sink = nx.single_source_dijkstra(_graph,0,8)
	print(path_sink)
	_rscgraph = nx.Graph(rscG)
	rval_lrscpath, rscpath_sink = nx.single_source_dijkstra(_rscgraph,0,8)


	LBw = 0
	for e in path_sink:
		LBw += w[e]
	print('reward on the min-risk path: ', LBw)

	budget = 0
	for i in range(0,len(path_sink)-1):
		budget += rscG[path_sink[i]][path_sink[i+1]]
	budget *= 2.0

	total_rew = np.sum(w)
	t0 = time.time()
	path = solve_CSP(G,  w, 0, 8, (1.0+gamma) * LBw, rscG, budget)
	print("Python Time:",time.time() - t0)


	''' preprocess the inner path for Linkernighan Heuristic '''
	# nb_obps = 70
	# for idx in range(1,2):
	# 	gamma = 0.3
	# 	dir = '/home/cai/Dropbox/Box_Research/Github/RRARP_LinKern/dat/InnerGraphs/'
	# 	filename = 'tG_'+ str(nb_obps) +'_16_' + str(idx) +'.dat'
	# 	G, x, y, w = read_pccsp_instance(dir+filename)
	# 	plot_ObP_graph([1,1,1], x, y, w)
	# 	_graph = nx.Graph(G)
	# 	wilename = '/home/cai/Dropbox/Box_Research/Github/RRARP_LinKern/dat/InnerPaths/1kCSP-' + str(idx)+'.dat'
	# 	wile = open(wilename, "w")
	# 	for i in range(0,16):
	# 		for j in range(i+1,16):
	# 			rval_lrpath, path_sink = nx.single_source_dijkstra(_graph, i, j)
	# # print(path_sink)
	# 			LBw = 0
	# 			for e in path_sink:
	# 				LBw += w[e]
	# # print(LBw)
	# 			total_rew = np.sum(w)
	# 			# t0 = time.time()
	# 			path = solve_CSP(G, w, i, j, gamma * (total_rew-LBw) + LBw)
	# 			pathstring = ''
	# 			for e in path:
	# 				pathstring += str(e) + ','
	# 			wile.write(str(idx) + '\t' + str(i) + '\t' + str(j) + '\t' + pathstring + '\n')
	# 	wile.close()





# ... Propogate from source 0 to node 1
# ... Propogate from source 0 to node 15
# ... Propogate from source 0 to node 17
# ... Propogate from source 0 to node 34
# ... Propogate from source 0 to node 40
# ... Propogate from source 0 to node 70
#  =====>>>> Time: 3.54356
#  =====>>> total reward: 503, Demand: 114, Reward on min-risk path: 76
#  =====>>> Budget: 104
#  =====>>> Optimal path [0 17 73 58 98 45 82 90 99 59 64 33 48 42 16 66 36 62 18 8 ]
#  =====>>> optimal obj value: 2.502, collected reward: 114 used resource: 78
# 0-8 1?2.502 = 0 17 73 58 98 45 82 90 99 59 64 33 48 42 16 66 36 62 18 8 











	# G = np.array([[0, 2, 3, 0, 0, 0, 0], 
	#               [2, 0, 9, 1, 2, 0, 0], 
	#               [3, 9, 0, 1, 0, 1, 0], 
	#               [0, 1, 1, 0, 1, 1, 0], 
	#               [0, 2, 0, 1, 0, 5, 2],
	#               [0, 0, 1, 1, 5, 0, 1],
	#               [0, 0, 0, 0, 2, 1, 0]])
	# w = [1, 2, 3, 1, 4, 2, 1]

	# path = solve_CSP(G, w, 0, 1, 1.0)

