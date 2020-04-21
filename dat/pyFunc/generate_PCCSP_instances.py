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
from operator import itemgetter 
from scipy.ndimage.filters import gaussian_filter1d


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

def generate_trigraph_random(dir):
	
	total_insts = 1000
	nb_obps = 100
	filename = 'trigraph_k16_'
	for i in range(1,total_insts+1):
		print("trigraph ", i)
		file = open(dir + str(i) + '.trigraph_k16', "w")

		''' write observation points '''
		K = 16
		precision = 1000.0
		ObPx = []
		ObPy = []
		ObPw = [] # reward
		R = []
		A = []
		center =  [1.0, 1.0, 1.0]
		'''write K boundary points '''
		for j in range(0, K):
			r = 1.0
			angle = round(2.0*np.pi/K * j *precision)/precision
			rewp = 0
			R.append(r)
			A.append(angle)
			x = r * math.cos(angle) + center[0]
			y = r *math.sin(angle) + center[1]
			ObPx.append(x)
			ObPy.append(y)
			ObPw.append(rewp)

		for j in range(0, nb_obps-K - 1):
			r = round(math.sqrt(np.random.uniform(0.,0.9))*precision)/precision
			angle = round(np.random.uniform(0,2.0*np.pi)*precision)/precision
			rewp = round(np.random.uniform(0.1,0.2)*precision)/precision
			R.append(r)
			A.append(angle)
			x = r * math.cos(angle) + center[0]
			y = r * math.sin(angle) + center[1]
			ObPx.append(x)
			ObPy.append(y)
			ObPw.append(rewp)

		r = round(math.sqrt(np.random.uniform(0.0,0.9))*precision)/precision
		angle = round(np.random.uniform(0,2.0*np.pi)*precision)/precision
		rewp = round(np.random.uniform(0.1,0.2)*precision)/precision
		R.append(r)
		A.append(angle)
		x = r * math.cos(angle) + center[0]
		y = r * math.sin(angle) + center[1]
		ObPx.append(x)
		ObPy.append(y)
		ObPw.append(rewp)		
		
		x = np.array(ObPx)
		y = np.array(ObPy)
		triang = mtri.Triangulation(x, y)
		file.write(str(nb_obps)+'\t' + str(len(triang.edges)) + '\n')

		''' write reward on nodes '''
		for i in range(len(ObPw)-1):
			file.write(str(i)+':'+str(R[i])+':' +str(A[i])+':'+str(ObPw[i])+'\t')
		file.write(str(len(ObPw)-1)+':'+str(R[-1])+':' +str(A[-1])+':'+str(ObPw[-1])+'\n')

		''' write weight info on edges '''
		for e in triang.edges:
			# print(e[0], e[1])
			dx = ObPx[e[0]] - ObPx[e[1]]
			dy = ObPy[e[0]] - ObPy[e[1]]
			dist = math.sqrt(dx**2 + dy**2)
			dist = round(dist*precision)/precision
			file.write(str(e[0]) + ':' + str(e[1]) + ':' + str(dist) + '\t')
		file.close()
		# plot_ObP_graph(center, x, y, ObPw)

if __name__ == "__main__":
	dir = '/home/cai/Dropbox/Box_Research/Github/RRARP_LinKern/dat/InnerGraphs/'
	generate_trigraph_random(dir)

