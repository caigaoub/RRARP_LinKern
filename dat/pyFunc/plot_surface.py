import gurobipy as gp
from gurobipy import GRB
from itertools import combinations

from sys import argv
import re
import math
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import matplotlib.tri as mtri
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Circle, PathPatch
from matplotlib.text import TextPath
from matplotlib.transforms import Affine2D
import mpl_toolkits.mplot3d.art3d as art3d


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
	# for i in range(len(obp_x)):
	# 	for j in range(len(obp_y)):
	# 		if (obp_x[i]-obp_x[j])**2+(obp_y[i]-obp_y[j])**2 > 0:
	# 			dist=math.sqrt((obp_x[i]-obp_x[j])**2+(obp_y[i]-obp_y[j])**2)
	# 			if dist < 0.5:
	# 				plt.plot([obp_x[i],obp_x[j]], [obp_y[i],obp_y[j]],'b', alpha=0.2)


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


def plot_target_area3D_general():
	precision = 1000.0
	center = [2,2,1.0]
	r = []
	angle = []
	rewp = []
	nbp = 500
	x = []
	y = []
	for i in range(0, nbp):
		r.append(round(math.sqrt(np.random.uniform(0.0,0.99))*precision)/precision)
		angle.append(round(np.random.uniform(0,2.0*np.pi)*precision)/precision)
		rewp.append(round(np.random.uniform(1,10)*precision)/precision)
		x.append(r[i] * math.cos(angle[i]) + center[0])
		y.append(r[i] * math.sin(angle[i]) + center[1])

	x = np.array(x) * 50
	y = np.array(y) * 50
	z = np.sinc((1-x)/100*3.14) + np.sinc((1-y)/100*3.14)+0.3


	fig = plt.figure()
	ax = fig.add_subplot(111,projection='3d')
	# ax = fig.add_subplot(111)
	plt.axis('off')
	ax.set_xlim([50*(center[0]-center[2]-0.1), 50*(center[0]+center[2]+0.1)])
	ax.set_ylim([50*(center[1]-center[2]-0.1), 50*(center[1]+center[2]+0.1)])
	ax.set_zlim(min(z), max(z))

	# ax.set_aspect('equal', adjustable='box')
	# big circle 
	p = Circle((100, 100), 50,fill=None, alpha = 0.5)
	ax.add_patch(p)
	art3d.pathpatch_2d_to_3d(p, z=0)

	# circle = plt.Circle((50*center[0], 50*center[1]), radius=center[2], edgecolor='k',fill=False, alpha=0.5)
	# ax.add_artist(circle)
	#observation points
	colors = rewp
	# area = np.array([np.pi*20.0]*len(x))
	# c = ax.scatter(x, y, c=colors, s=area, cmap='hsv', alpha=0.85)
	# for i, txt in enumerate(range(1,len(x)+1)):
	# 	# ax.annotate(txt, (obp_x[i]+0.02, obp_y[i]+0.02),size=10)
	# 	ax.annotate(txt,  (x[i]+0.02, y[i]-0.02), horizontalalignment='right', verticalalignment='top',size=10)
	# for i in range(len(obp_x)):
	# 	for j in range(len(obp_y)):
	# 		if (obp_x[i]-obp_x[j])**2+(obp_y[i]-obp_y[j])**2 > 0:
	# 			dist=math.sqrt((obp_x[i]-obp_x[j])**2+(obp_y[i]-obp_y[j])**2)
	# 			if dist < 0.5:
	# 				plt.plot([obp_x[i],obp_x[j]], [obp_y[i],obp_y[j]],'b', alpha=0.2)

	triang = mtri.Triangulation(x, y)
	ax.triplot(triang, c="#D3D3D3", marker='.', markerfacecolor="black", markeredgecolor="black", markersize=5)
	# ax.triplot(triang, c="#D3D3D3", marker='.', markerfacecolor="black", markeredgecolor="black", markersize=5)

	ax.plot_trisurf(triang, z, cmap='jet')
	ax.scatter(x,y,z, marker='.', s=20, c="black", alpha=0.5)
	ax.view_init(elev=60, azim=-45)

	plt.show()


def plot_area2D_randObPs():
	precision = 1000.0
	center = [2,2,1.0]
	r = []
	angle = []
	rewp = []
	nbp = 500
	x = []
	y = []
	for i in range(0, nbp):
		r.append(round(math.sqrt(np.random.uniform(0.0,0.99))*precision)/precision)
		angle.append(round(np.random.uniform(0,2.0*np.pi)*precision)/precision)
		rewp.append(round(np.random.uniform(1,10)*precision)/precision)
		x.append(r[i] * math.cos(angle[i]) + center[0])
		y.append(r[i] * math.sin(angle[i]) + center[1])
	x.append(math.cos(0) + center[0])
	y.append(math.sin(0) + center[1])

	x.append(math.cos(1.25*np.pi) + center[0])
	y.append(math.sin(1.25*np.pi) + center[1])

	x = np.array(x) * 50
	y = np.array(y) * 50

	fig = plt.figure()
	ax = fig.add_subplot(111)
	plt.axis('off')
	ax.set_xlim([50*(center[0]-center[2]-0.1), 50*(center[0]+center[2]+0.1)])
	ax.set_ylim([50*(center[1]-center[2]-0.1), 50*(center[1]+center[2]+0.1)])
	ax.set_aspect('equal', adjustable='box')
	# big circle 
	circle = plt.Circle((100, 100), radius=50, edgecolor='k',fill=False, alpha=0.5)
	ax.add_artist(circle)
	triang = mtri.Triangulation(x, y)
	ax.triplot(triang, c="#D3D3D3", marker='.', markerfacecolor="black", markeredgecolor="black", markersize=5)
	circle = plt.Circle((x[-1], y[-1]), radius=1, edgecolor='b',facecolor='b')
	ax.add_artist(circle)
	circle = plt.Circle((x[-2], y[-2]), radius=1, edgecolor='b',facecolor='b')
	ax.add_artist(circle)
	print(x)
	print(triang.edges)
	plt.show()


def plot_ObP_graph(center, obpx, obpy, obprew):
	x = np.array(obpx)
	y = np.array(obpy)
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
		for i in range(16):
			circle = plt.Circle((x[i], y[i]), radius=0.02, edgecolor='k',facecolor='k')
			ax.add_artist(circle)

	plt.show()


if __name__ == "__main__":
	instance = argv[1]

	tar_idx = 1
	'''two depots, n target areas, reward pct, all observation points, and reward values on Obser. point '''
	departure, arrival, tar_locs, tar_rews, OBPX, OBPY, OBPREW = read_instance(instance)

	plot_ObP_graph(tar_locs[tar_idx], OBPX[tar_idx], OBPY[tar_idx], OBPREW[tar_idx])



	''' test '''
	# plot_area2D_randObPs()







	#https://fabrizioguerrieri.com/blog/2017/9/7/surface-graphs-with-irregular-dataset
