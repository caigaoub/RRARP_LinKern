import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as pch
import re
import math
import time
from sys import argv

def set_plot_attributes(plt, ax):
    plt.grid(alpha=.5)
    # plt.xticks(np.arange(-10,1000,1))
    # plt.yticks(np.arange(-10,1000,1))
    ax.set(xlim=(0,1000), ylim=(0,1000))
    # ax.set_aspect('equal', adjustable='box')


def plot_instance(instancefile):
    file_ = open(instancefile, 'r')
    line_ = file_.readline()
    list_ = re.split("\t|\n", line_)
    nb_tars_ = int(list_[0])
    fig, ax  = plt.subplots()
    # set_plot_attributes(plt, ax)
    ax.set_title(instancefile)
    minx = 100000
    maxx = -100
    miny = 100000
    maxy = -100
    #departure
    line_ = file_.readline()
    list_ = re.split("\t|\n", line_)    
    departure = [float(list_[0])]
    departure.append(float(list_[1]))
    minx = min(minx, departure[0]-3)
    maxx = max(maxx, departure[0]+3)
    miny = min(miny, departure[1]-3)
    maxy = max(maxy, departure[1]+3)
    circle = pch.Circle((departure[0], departure[1]), radius=0.1, edgecolor='r', facecolor='r', alpha=1)
    ax.add_artist(circle)
    ax.annotate("("+str(0)+": " + str(departure[0])+","+str(departure[1])+")", xy=(departure[0], departure[1]), xytext=(departure[0]+0.3, departure[1]+0.3))

    #arrival
    line_ = file_.readline()
    list_ = re.split("\t|\n", line_)    
    arrival = [float(list_[0])]
    arrival.append(float(list_[1]))
    minx = min(minx, arrival[0]-3)
    maxx = max(maxx, arrival[0]+3)
    miny = min(miny, arrival[1]-3)
    maxy = max(maxy, arrival[1]+3)
    circle = pch.Circle((arrival[0], arrival[1]), radius=0.1, edgecolor='r', facecolor='r', alpha=1)
    ax.add_artist(circle)
    ax.annotate("("+str(nb_tars_+1)+": "+str(arrival[0])+","+str(arrival[1])+")", xy=(arrival[0], arrival[1]), xytext=(arrival[0]+0.3, arrival[1]+0.3))

    #target circles
    itr = 1
    line_ = file_.readline()
    while itr <= nb_tars_:
        list_ = re.split("\t|\n", line_)
        tar_loc_x = float(list_[0])
        tar_loc_y = float(list_[1])
        tar_rad = float(list_[2])
        # print(tar_loc_x,tar_loc_y,tar_rad)
        minx = min(minx, tar_loc_x-tar_rad)
        maxx = max(maxx, tar_loc_x+tar_rad)
        miny = min(miny, tar_loc_y-tar_rad)
        maxy = max(maxy, tar_loc_y+tar_rad)
        circle = plt.Circle((tar_loc_x, tar_loc_y), radius=tar_rad, edgecolor='k',facecolor='k',alpha=0.2)
        ax.add_artist(circle)
        circle = plt.Circle((tar_loc_x, tar_loc_y), radius=0.01, edgecolor='r',facecolor='r',alpha=0.8)
        ax.add_artist(circle)
        ax.annotate("("+str(itr)+": "+str(tar_loc_x)+","+str(tar_loc_y)+")", xy=(tar_loc_x, tar_loc_y), xytext=(tar_loc_x+0.3, tar_loc_y+0.3))
        itr += 1   
        line_ = file_.readline()
    file_.close()
    ax.set(xlim=(minx,maxx), ylim=(miny,maxy))
    ax.set_aspect('equal', adjustable='box')
    # plt.xticks(np.arange(minx,maxx,1))
    # plt.yticks(np.arange(miny,maxy,1))
    plt.grid(alpha=.5)
    plt.show()

################################################################################################

#  transformation files
# instancefile = '../ret/inst_n_10/n_10_h_6.txt'
# instancefile = '../ret/cluster_n_30/n_30_c_5_1.txt'
# instancefile = '/home/cai/Dropbox/Box_Research/Github/RRARP_BD/BendersDecomp/RRARP-BD/dat/test_n_6.txt'
# instancefile = '/home/caigao/Dropbox/Box_Research/Github/RRARP_BD/BendersDecomp/RRARP-BD/dat/test_n_6.txt'

instancefile = argv[1]
plot_instance(instancefile)
