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
import recursive_search_CSP as rs
import time
import cProfile

def find_circle(selected, n):
	visited = [False]*n
	# print(" matrix = ", selected)
	visited[0] = True
	cur_node = 0
	''' find the feasible path from 0 to n-1 '''
	path = [cur_node]
	while True:
		for i in range(0, n):
			if selected[cur_node, i] > 0.5 and visited[i] == False:
				path.append(i)
				visited[i] = True
				cur_node = i
		if cur_node == n-1:
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
		
		circle = find_circle(vals, n)
		if (len(circle) != 0):
			constr = 0
			for i in range(0,len(circle)-1):
				constr += model._vars[circle[i],circle[i+1]]
			constr += model._vars[circle[-1],circle[0]]
			model.cbLazy(constr <= len(circle)-1)


''' ----------------------------------------- '''
''' solve the constrained shortest path model '''
''' ----------------------------------------- '''

def solve_CSP(riskmatx, obp_rews, rew_lb):
	# print(riskmatx)
	n = riskmatx.shape[0]
	# print(n)
	total_rew = sum(obp_rews)
	# make rewards on source and sink zero
	# obp_rews.append(0)
	# obp_rews.insert(0, 0)    
	try:
		# Create a new model
		model = gp.Model("CSP")
		model.setParam(GRB.Param.OutputFlag, 1)
		model.setParam(GRB.Param.TimeLimit, 3600.0)

		x = model.addVars(n, n, vtype=GRB.BINARY,name="x")
		constr = 0
		for i in range(0,n):
			for j in range(0,n):
				if riskmatx[i][j] > 0:
					# print(riskmatx[i][j], end=' ')
					constr += riskmatx[i][j] * x[i,j]
			# print('\n')
		model.setObjective(constr, GRB.MINIMIZE)
		# model.modelSense = GRB.MINIMIZE

		# source node
		constr = 0 
		for i in range(1, n-1):
			if riskmatx[0][i] > 0:
				constr += x[0,i]         
		model.addConstr(constr == 1, 's1') 
		constr = 0 
		for i in range(1, n-1):
			if riskmatx[i][0] > 0:
				constr += x[i,0]         
		model.addConstr(constr == 0, 's2') 
		# sink node
		constr = 0 
		for i in range(1, n-1):
			if riskmatx[i][n-1] > 0:
				constr += x[i,n-1]         
		model.addConstr(constr == 1, 't1') 
		constr = 0 
		for i in range(1, n-1):
			if riskmatx[n-1][i] > 0:
				constr += x[n-1,i]       
		model.addConstr(constr == 0, 't2') 

		for i in range(1,n-1):
			constr = 0 
			for j in range(0,n):
				if riskmatx[i][j] > 0:
					constr += x[i,j] 
			for j in range(0,n):
				if riskmatx[j][i] > 0:
					constr -= x[j,i] 
			model.addConstr(constr == 0, "flow_"+str(i))
			model.update()

		for i in range(1,n-1):
			constr = 0 
			for j in range(0,n):
				if riskmatx[i][j] > 0:
					constr += x[i,j]
				if riskmatx[j][i] > 0:
					constr += x[j,i] 
			model.addConstr(constr <= 2, "deg_"+str(i))
			model.update()

		# for i in range(0,n):
		#     model.addConstr(x[i,i] == 0, "_"+str(i))
		# model.update()
	 
		for i in range(0, n):
			for j in range(0, n):
				if abs(riskmatx[i][j]) < 0.000000001:
					model.addConstr(x[i,j] == 0, "force"+str(i))
		model.update()


		constr = 0 
		for i in range(0,n):
			for j in range(0,n):
					if riskmatx[i][j] > 0:
						constr += x[i,j] * (obp_rews[i] + obp_rews[j])
		model.addConstr(constr/2.0 >= total_rew * rew_lb , "rc")
		# model.addConstr(constr >= 14.1, "rc")

		model.update()

		# model.write('model.lp')

		model._vars = x
		model.Params.lazyConstraints = 1
		model.optimize(cb_subtourelim)

		''' ------------- model output  ----------------------'''
		if model.status == GRB.OPTIMAL or model.status == GRB.TIME_LIMIT:
			# for i in range(0, n):
			#     for j in range(0, n):
			#         print(int(x[i,j].x), end = ' ')
			#     print('')
			visited = [False]*n
			visited[0] = True
			cur_node = 0
			path = [cur_node]
			while True:
				for i in range(0, n):
					if x[cur_node,i].x > 0.5 and visited[i] == False:
						path.append(i)
						visited[i] = True
						cur_node = i
				if cur_node == n-1:
					break
			print("the optimal path: ", path)
			collec_r = 0
			for e in path:
				collec_r += obp_rews[e]
			print('collected rewards: ', collec_r, ' requirment:', total_rew * rew_lb)
			total_risk = 0
			for i in range(len(path)-1):
				total_risk += riskmatx[path[i]][path[i+1]]
			print('total risk over optimal path: %g' % total_risk)

	except gp.GurobiError as e:
		print('Error code ' + str(e.errno) + ': ' + str(e))

	except AttributeError:
		print('Encountered an attribute error')

	return path


def read_instance(instancefile):
	file_ = open(instancefile, 'r')
	line_ = file_.readline()
	list_ = re.split("\t|\n", line_)
	nb_tars_ = int(list_[0])
	#departure
	line_ = file_.readline()
	list_ = re.split("\t|\n", line_)    
	departure = [float(list_[0])]
	departure.append(float(list_[1]))
	#arrival
	line_ = file_.readline()
	list_ = re.split("\t|\n", line_)    
	arrival = [float(list_[0])]
	arrival.append(float(list_[1]))
	#target circles
	itr = 1
	line_ = file_.readline()
	tar_locs = []
	while itr <= nb_tars_:
		list_ = re.split("\t|\n", line_)
		tar_loc_x = float(list_[0])
		tar_loc_y = float(list_[1])
		tar_rad = float(list_[2])
		tar_locs.append([tar_loc_x, tar_loc_y, tar_rad])

		itr += 1   
		line_ = file_.readline()
	''' rewards shreshold '''
	list_ = re.split("\t|\n", line_)
	tar_rews = [0]
	for e in list_:
		if e != '':
			tar_rews.append(float(e))
	tar_rews.append(0)
	''' print observation point '''
	line_ = file_.readline()
	if True:
		tar_idx = 1
		OBPX = []
		OBPY = []
		OBPREW = []
		while tar_idx <= nb_tars_:
			list_ = re.split("\t|\n", line_)
			obpX = []
			obpY = []
			obpRew = []
			# print(list_)
			for e in list_:
				if e != '':
					list2_ = re.split(":", e)
					r = float(list2_[0])
					angle = float(list2_[1])
					obp_x = r * math.cos(angle) + tar_locs[tar_idx-1][0]
					obp_y = r * math.sin(angle) + tar_locs[tar_idx-1][1]
					obp_rew = float(list2_[2])
					obpX.append(obp_x)
					obpY.append(obp_y)
					obpRew.append(obp_rew)
			tar_idx += 1
			OBPX.append(obpX)
			OBPY.append(obpY)
			OBPREW.append(obpRew)
			line_ = file_.readline()

	file_.close()
	return departure, arrival, tar_locs, tar_rews, OBPX, OBPY,OBPREW

def plot_target_area(center, obp_x, obp_y, obp_reward, entryP, exitP, path=None):
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
	circle = plt.Circle((entryP[0], entryP[1]), radius=0.02, edgecolor='k',facecolor='k')
	ax.add_artist(circle)
	circle = plt.Circle((exitP[0], exitP[1]), radius=0.02, edgecolor='k',facecolor='k')
	ax.add_artist(circle)

	plt.show()

def read_config(config):
	file_ = open(config, 'r')
	line_ = file_.readline()
	return line_

if __name__ == "__main__":
	# ccr 
	configfile = argv[1]
	position = argv[2]
	path = ''
	if int(position) == 1:
		path = '/home/cai/Dropbox/Box_Research/Github/RRARP_LinKern/'
	else:
		path = '/projects/academic/josewalt/caigao/RRARP_LinKern/'
	instancefile = path + 'dat/' + read_config(configfile)

	# instancefile = argv[1]
	departure, arrival, tar_locs, tar_rews, OBPX, OBPY, OBPREW = read_instance(instancefile)
	# print(departure, arrival, tar_locs, tar_rews, OBPX, OBPY)

	nb_obserPts =30
	size_K = 16
	subangle = 2.0*np.pi/size_K;
	''' build the risk matrix '''
	riskmatx = np.zeros(shape=(len(tar_locs), nb_obserPts+2, nb_obserPts+2))
	for t in range(0, len(tar_locs)):
		for i in range(0, nb_obserPts):
			for j in range(i, nb_obserPts):
				if i != j:
					if math.sqrt((OBPX[t][i] - OBPX[t][j])**2 + (OBPY[t][i] - OBPY[t][j])**2) < 0.5:
						riskmatx[t][i+1][j+1] = math.sqrt((OBPX[t][i] - OBPX[t][j])**2 + (OBPY[t][i] - OBPY[t][j])**2)
						riskmatx[t][j+1][i+1] = riskmatx[t][i+1][j+1]
					# riskmatx[t][i][j] = (OBPX[t][i-1] - OBPX[t][j-1])**2 + (OBPY[t][i-1] - OBPY[t][j-1])**2
				# else:
				#     riskmatx[t][i][j] = 0.0


	# print(riskmatx[0])

	# print(type(instancefile))
	instancename = re.split(".dat", re.split("/", instancefile)[-1])[0]
	file = open(path + 'dat/InnerPaths/' + instancename+'.path', "w")
	file.write('target' + '\t' + 'entry' + '\t' + 'exit' + '\t' + 'riskval' + '\t' + 'path' +'\n')
	for t in range(4, len(tar_locs)): 
		print(" working on target ", t)
		for i in range(0, size_K):
			for j in range(2, size_K):
				print("entry: ", i, "exit: ", j)
				''' update risk matrix '''
				entryP = [tar_locs[t][0] + tar_locs[t][2]*math.cos(subangle*i), tar_locs[t][1] + tar_locs[t][2]*math.sin(subangle*i)];
				exitP = [tar_locs[t][0] + tar_locs[t][2]*math.cos(subangle*j), tar_locs[t][1] + tar_locs[t][2]*math.sin(subangle*j)];
				for ii in range(1, nb_obserPts+1):
					if math.sqrt((entryP[0] - OBPX[t][ii-1])**2 + (entryP[1] - OBPY[t][ii-1])**2) < 1.0:
						riskmatx[t][0][ii] = math.sqrt((entryP[0] - OBPX[t][ii-1])**2 + (entryP[1] - OBPY[t][ii-1])**2)
						riskmatx[t][ii][0] = riskmatx[t][0][ii]
					if math.sqrt((exitP[0] - OBPX[t][ii-1])**2 + (exitP[1] - OBPY[t][ii-1])**2) < 1.0:
						riskmatx[t][ii][nb_obserPts+1] = math.sqrt((exitP[0] - OBPX[t][ii-1])**2 + (exitP[1] - OBPY[t][ii-1])**2)
						riskmatx[t][nb_obserPts+1][ii] = riskmatx[t][ii][nb_obserPts+1]
				# path = solve_CSP(riskmatx[t], OBPREW[t][0:nb_obserPts], tar_rews[t])

				budget_pct = 0.5
				OP_rews = OBPREW[t][0:nb_obserPts]
				OP_rews.append(0)
				OP_rews.insert(0, 0)   
				t0 = time.time()
				# print(riskmatx[t])
				path = solve_CSP(riskmatx[t], OP_rews, budget_pct)
				print("Model Time:",time.time() - t0)

				# G = nx.Graph(riskmatx[t])
				# minrisk_sink = nx.single_source_dijkstra_path_length(G,source=len(G.nodes)-1)
				
				# test = rs.Pulse(G, budget_pct, minrisk_sink, OP_rews)
				# t0 = time.time()
				# test.recursive_search(0, [], 0.0, 0.0,False)
				# test.print_optimal_sol()
				# print("Pulse Algorithm Time:",time.time() - t0)
				total_risk = 0
				for idx in range(0, len(path)-1):
					total_risk += riskmatx[t][path[idx]][path[idx+1]]
				total_risk = round(total_risk*1000.0)/1000.0
				path_str = ''
				for e in path[1:-1]:
					path_str += str(e)+','    
				tar_idx = t	
				file.write(str(t) + '\t' + str(i) + '\t' + str(j) + '\t' + str(total_risk) + '\t' + path_str + '\n')
				plot_target_area(tar_locs[tar_idx], OBPX[tar_idx][0:nb_obserPts], OBPY[tar_idx][0:nb_obserPts],OBPREW[tar_idx][0:nb_obserPts], entryP, exitP, path)
				# plot_target_area(tar_locs[tar_idx], OBPX[tar_idx][0:nb_obserPts], OBPY[tar_idx][0:nb_obserPts],OBPREW[tar_idx][0:nb_obserPts], entryP, exitP)
				
				exit()

	file.close()






	# risktest = np.array([[0, 2, 3, 0, 0, 0, 0], 
	#                      [2, 0, 9, 1, 2, 0, 0], 
	#                      [3, 9, 0, 1, 0, 1, 0], 
	#                      [0, 1, 1, 0, 1, 1, 0], 
	#                      [0, 2, 0, 1, 0, 5, 2],
	#                      [0, 0, 1, 1, 5, 0, 1],
	#                      [0, 0, 0, 0, 2, 1, 0]])
	# OBPREWtest = [2, 3, 1, 4, 2]
	# solve_CSP(risktest, OBPREWtest, 0.6)



