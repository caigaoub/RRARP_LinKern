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
		if False:
			for i in range(0, model._x_size):
				for j in range(0, model._x_size):
					if int(vals[i,j]) > 0.5:
						print(i,' ', j, end= ', ')
				# print('')    
			# print('\n')
		# Cirset = find_circle(vals, n, model._entry, model._exit)
		Cirset = find_all_circles(vals, model._x_size, model._entry, model._exit)
		if (len(Cirset) != 0):
			for circle in Cirset:
				constr = 0
				for i in range(0,len(circle)-1):
					constr += model._vars[circle[i],circle[i+1]]
				constr += model._vars[circle[-1],circle[0]]
				model.cbLazy(constr <= len(circle)-1)
				model._nb_SECs  += 1

''' ----------------------------------------- '''
''' solve Prize-Collection Constrained Shortest Path (PCCSP) model '''
''' ----------------------------------------- '''
def solve_CSP(G, obp_rews, entry, exit, demand):
	n = G.shape[0]
	total_rew = np.sum(obp_rews)
	if obp_rews[entry] + obp_rews[exit] >= demand:
		return [entry, exit]
	try:
		# Create a new model
		model = gp.Model("PCCSP")
		model.setParam(GRB.Param.OutputFlag, 1)
		model.setParam(GRB.Param.TimeLimit, 200.0)
		model._entry = entry
		model._exit = exit
		model._x_size = n
		model._nb_SECs = 0
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

		# model.write('model.lp')

		model._vars = x
		model.Params.lazyConstraints = 1
		model.optimize(cb_subtourelim)
		# model.optimize()

		''' ------------- model output  ----------------------'''
		if model.status == GRB.OPTIMAL or model.status == GRB.TIME_LIMIT:
			# vals = np.zeros((n,n))
			# for i in range(n):
			# 	# print(i, '= ', end=' ')
			# 	for j in range(n):
			# 		# print(int(x[i,j].x), end=' ')
			# 		vals[i,j] = int(x[i,j].x)
			# 		if int(x[i,j].x) > 0.5:
			# 			print(i, ' ', j, end=', ')
			# print('dddddddd')
			# vals = model.cbGetSolution(model._vars)
			# find_all_circles(vals, n, entry, exit)
			# print('\n')
			visited = [False]*n
			visited[entry] = True
			cur_node = entry
			path = [cur_node]
			
			while True:
				for i in range(0, n):
					# print(cur_node, x[cur_node,i].x, visited[i], end=', ' )
					if x[cur_node,i].x > 0.5 and visited[i] == False:
						# print(path, cur_node)
						path.append(i)
						visited[i] = True
						cur_node = i
				# print('\n')
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
			# print(' & ' + str(int(model.Runtime *100)/100) + ' & ' + str(int(total_risk*100)/100) + ' & ' + str(int(model.ObjBound*100)/100) + ' & ' + str(int(model.MIPGap*100)/100))
			# print(' & ' + str("{:.2f}".format(model.Runtime)) + ' & ' + str(model._nb_SECs) + ' & ' + str("{:.2f}".format(model.MIPGap*100)))

	except gp.GurobiError as e:
		print('Error code ' + str(e.errno) + ': ' + str(e))

	except AttributeError:
		print('Encountered an attribute error')

	return total_risk, path, "{:.2f}".format(model.MIPGap*100)

def read_pccsp_instance(instancefile):
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
	# plot_ObP_graph([1,1,1], x, y, w)
	return G, x, y, w

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


def preprocess():
	nb_instances = 10
	dir = '/home/cai/Dropbox/Box_Research/Github/RRARP_LinKern/dat/InfoRegions/'
	EEpairs = []
	for i in range(0,8):
		for j in range(i+1,8):
			EEpairs.append([i,j])
	for v in range(1,501):
		filename = 'tG_70_8_' + str(v) +'.dat'
		G, x, y, w = read_pccsp_instance(dir+filename)
		_graph = nx.Graph(G)
		# maxw_spath = 0
		# for ee in EEpairs:
		# 	rval_lrpath, path_sink = nx.single_source_dijkstra(_graph,ee[0],ee[1])
		# 	LBw = 0
		# 	for e in path_sink:
		# 		LBw += w[e]
		# 	if LBw > maxw_spath:
		# 		maxw_spath = LBw
		# demand = 1 * maxw_spath
		demand = 0.2 * np.sum(w)
		file = open("/home/cai/Dropbox/Box_Research/Github/RRARP_LinKern/dat/InfoRegions/preprocess_results/prepro_"+str(v)+'.out', "w")
		# k = 0
		# print(EEpairs)
		for ee in EEpairs:
			# rval_lrpath, path_sink = nx.single_source_dijkstra(_graph,ee[0],ee[1])
			print(ee)
			UB, path, mipgap = solve_CSP(G, w, ee[0], ee[1], demand)
			file.write(str(v) + ' ' + str(ee[0]) + ' ' + str(ee[1]) + ' ' +  str(mipgap) + '?' +  str(UB) + " = ")
			for e in path:
				file.write(str(e) + ',')
			file.write('\n')
			# if k > 2:
			# 	break
			# k += 1

		file.close()


def run_one_instance(nb_obps, idx, gamma):
	dir = '/home/cai/Dropbox/Box_Research/Github/RRARP_LinKern/dat/Inst_Pulse/'
	filename = 'tG_'+ str(nb_obps) +'_2_' + str(idx) +'.dat'
	t0 = time.time()

	G, x, y, w = read_pccsp_instance(dir+filename)
	_graph = nx.Graph(G)
	rval_lrpath, path_sink = nx.single_source_dijkstra(_graph,0,1)
	# print(path_sink)
	LBw = 0
	for e in path_sink:
		LBw += w[e]
	
	total_rew = np.sum(w)
	print("max gamma: ", total_rew/LBw)

	path = solve_CSP(G, w, 0, 1, gamma * LBw)
	path = solve_CSP(G, w, 0, 1, gamma * LBw)


if __name__ == "__main__":
	# nb_obps = int(argv[1])
	# idx = int(argv[2])
	# gamma = float(argv[3])
	# run_one_instance(nb_obps, idx, gamma)
	preprocess()


	# print("Python Time:",time.time() - t0)


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

















	# G = np.array([[0, 2, 3, 0, 0, 0, 0], 
	#               [2, 0, 9, 1, 2, 0, 0], 
	#               [3, 9, 0, 1, 0, 1, 0], 
	#               [0, 1, 1, 0, 1, 1, 0], 
	#               [0, 2, 0, 1, 0, 5, 2],
	#               [0, 0, 1, 1, 5, 0, 1],
	#               [0, 0, 0, 0, 2, 1, 0]])
	# w = [1, 2, 3, 1, 4, 2, 1]

	# path = solve_CSP(G, w, 0, 1, 1.0)

