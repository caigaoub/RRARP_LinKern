
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
import sys


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def read_one_console(Rilename, Wilename):
	file_ = open(Rilename, 'r')
	line_ = file_.readline()
	nb_lines = 120

	wile = open(Wilename, 'w')
	for i in range(nb_lines):
		line_ = file_.readline()
		list_ = re.split(" ", line_)
		nodetuple = re.split("-", list_[0])
		# print(nodetuple)
		wile.write(nodetuple[0] + '\t' + nodetuple[1] + '\t' + nodetuple[2] + '\t')
		optimalstatus =list_[1].split('?')[0]
		# print(optimalstatus)
		if optimalstatus == '0':
			eprint("not optimal path!!")
			return

		objval = list_[1].split('?')[1]
		wile.write(objval + '\t')
		optpath = []
		for j in range(3, len(list_)-1):
			optpath.append(list_[j])
			wile.write(list_[j]+',')
		# print(list_)
		wile.write('\n')	
	file_.close()


if __name__ == "__main__":
	Rilename = '/home/cai/Dropbox/Box_Research/Github/RRARP_LinKern/dat/console/1kCSP-'
	Wilename = '/home/cai/Dropbox/Box_Research/Github/RRARP_LinKern/dat/InnerPaths/1kCSP-'
	for i in range(1, 1001):
		print(i)
		if i not in [25, 106, 151, 184,198,202,211,254]:
			Rname = '/home/cai/Dropbox/Box_Research/Github/RRARP_LinKern/dat/console/1kCSP-'+str(i)+'.out'
			Wname = '/home/cai/Dropbox/Box_Research/Github/RRARP_LinKern/dat/InnerPaths/1kCSP-'+str(i)+'.dat'
			read_one_console(Rname, Wname)