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

	for i, txt in enumerate(range(1,len(x)+1)):
		ax.annotate(txt-1,  (x[i]+0.02, y[i]-0.02), horizontalalignment='right', verticalalignment='top',size=10)
	plt.show()


def write_triGraph(instance, filename, nb_ObPs):
	departure, arrival, tar_locs, tar_rews, OBPX, OBPY, OBPREW = read_instance(instance)
	file = open(filename, "w")
	''' write observation points '''
	for i in range(len(tar_locs)):
		triang = mtri.Triangulation(OBPX[i][0:nb_ObPs], OBPY[i][0:nb_ObPs])
		file.write(str(len(triang.edges)) + '\t')
		for e in triang.edges:
			file.write(str(e[0]) + ':' + str(e[1]) + '\t')
		file.write('\n')
	
	file.close()

if __name__ == "__main__":
	instance = argv[1]
	nb_ObPs = int(argv[2])
	nb_ObPs += 16
	'''two depots, n target areas, reward pct, all observation points, and reward values on Obser. point '''
	departure, arrival, tar_locs, tar_rews, OBPX, OBPY, OBPREW = read_instance(instance)
	tar_idx = 4
	plot_ObP_graph(tar_locs[tar_idx], OBPX[tar_idx][0:nb_ObPs], OBPY[tar_idx][0:nb_ObPs], OBPREW[tar_idx][0:nb_ObPs])
	# triang = mtri.Triangulation(OBPX[tar_idx], OBPY[tar_idx])
	# print(triang.edges)




	# filename = '/home/cai/Dropbox/Box_Research/Github/RRARP_LinKern/dat/InnerGraphs/n_6_h_1.graph_200'
	# write_triGraph(instance, filename, 16+200)
	# filename = '/home/cai/Dropbox/Box_Research/Github/RRARP_LinKern/dat/InnerGraphs/n_6_h_1.graph_150'
	# write_triGraph(instance, filename, 16+150)
	# filename = '/home/cai/Dropbox/Box_Research/Github/RRARP_LinKern/dat/InnerGraphs/n_6_h_1.graph_120'
	# write_triGraph(instance, filename, 16+120)
	# filename = '/home/cai/Dropbox/Box_Research/Github/RRARP_LinKern/dat/InnerGraphs/n_6_h_1.graph_100'
	# write_triGraph(instance, filename, 16+100)
	# filename = '/home/cai/Dropbox/Box_Research/Github/RRARP_LinKern/dat/InnerGraphs/n_6_h_1.graph_80'
	# write_triGraph(instance, filename, 16+80)
	# filename = '/home/cai/Dropbox/Box_Research/Github/RRARP_LinKern/dat/InnerGraphs/n_6_h_1.graph_60'
	# write_triGraph(instance, filename, 16+60)
	# filename = '/home/cai/Dropbox/Box_Research/Github/RRARP_LinKern/dat/InnerGraphs/n_6_h_1.graph_50'
	# write_triGraph(instance, filename, 16+50)
	# filename = '/home/cai/Dropbox/Box_Research/Github/RRARP_LinKern/dat/InnerGraphs/n_6_h_1.graph_40'
	# write_triGraph(instance, filename, 16+40)
	# filename = '/home/cai/Dropbox/Box_Research/Github/RRARP_LinKern/dat/InnerGraphs/n_6_h_1.graph_30'
	# write_triGraph(instance, filename, 16+30)
	# filename = '/home/cai/Dropbox/Box_Research/Github/RRARP_LinKern/dat/InnerGraphs/n_6_h_1.graph_10'
	# write_triGraph(instance, filename, 16+10)



	# 0.158 0.105 0.172 0.156 0.157 0.115 0.16 0.12 0.108 0.157
