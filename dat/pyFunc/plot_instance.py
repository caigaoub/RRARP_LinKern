import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as pch
import re
import math
import time
from sys import argv
from operator import itemgetter 

def set_plot_attributes(plt, ax):
    plt.grid(alpha=.5)
    # plt.xticks(np.arange(-10,1000,1))
    # plt.yticks(np.arange(-10,1000,1))
    ax.set(xlim=(0,1000), ylim=(0,1000))
    # ax.set_aspect('equal', adjustable='box')


def plot_instance(instancefile, instsize, count):
    file_ = open(instancefile, 'r')
    line_ = file_.readline()
    list_ = re.split("\t|\n", line_)
    nb_tars_ = int(list_[0])
    fig, ax  = plt.subplots()
    # set_plot_attributes(plt, ax)
    # ax.set_title(instancefile)
    minx = 100000
    maxx = -100
    miny = 100000
    maxy = -100
    #departure depot
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
    # ax.annotate("("+str(0)+": " + str(departure[0])+","+str(departure[1])+")", xy=(departure[0], departure[1]), xytext=(departure[0]+0.3, departure[1]+0.3))

    #arrival depot
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
    # ax.annotate("("+str(nb_tars_+1)+": "+str(arrival[0])+","+str(arrival[1])+")", xy=(arrival[0], arrival[1]), xytext=(arrival[0]+0.3, arrival[1]+0.3))

    #target circles
    itr = 1
    line_ = file_.readline()
    tar_locs = []
    while itr <= nb_tars_:
        list_ = re.split("\t|\n", line_)
        tar_loc_x = float(list_[0])
        tar_loc_y = float(list_[1])
        tar_locs.append([tar_loc_x, tar_loc_y])
        tar_rad = float(list_[2])
        # print(tar_loc_x,tar_loc_y,tar_rad)
        minx = min(minx, tar_loc_x-tar_rad)
        maxx = max(maxx, tar_loc_x+tar_rad)
        miny = min(miny, tar_loc_y-tar_rad)
        maxy = max(maxy, tar_loc_y+tar_rad)
        circle = plt.Circle((tar_loc_x, tar_loc_y), radius=tar_rad, edgecolor='k',fill=False,alpha=0.3,lw=1)
        ax.add_artist(circle)
        # circle = plt.Circle((tar_loc_x, tar_loc_y), radius=0.01, edgecolor='r',facecolor='r',alpha=1)
        ax.add_artist(circle)
        ax.annotate(str(itr), xy=(tar_loc_x, tar_loc_y), xytext=(tar_loc_x-1.2, tar_loc_y+1),size=6)

        # ax.annotate("("+str(itr)+": "+str(tar_loc_x)+","+str(tar_loc_y)+")", xy=(tar_loc_x, tar_loc_y), xytext=(tar_loc_x+0.3, tar_loc_y+0.3))
        itr += 1   
        line_ = file_.readline()

    if True:
        dir = '/home/cai/Dropbox/Box_Research/Github/RRARP_LinKern/dat/InfoRegions/'
        # list_ = re.split(',',pickedtrigraphs)
        list_ = range(1, instsize+1);
        allX = []
        allY = []
        for tar in list_:
            if tar != '':
                filename = dir + 'tG_70_8_' + str(tar) + '.dat'
                file_ = open(filename, 'r')
                line_ = file_.readline()
                list_ = re.split(' |\t', line_)
                nb_obps = int(list_[0])
                line_ = file_.readline()
                list_ = re.split(' |\t|\n', line_)
                R = []
                A = []
                W = []
                X = []
                Y = []
                # print(list_)
                for el in list_:
                    if el != '':
                        list2_ = re.split(':', el)
                        x = float(list2_[1])
                        y = float(list2_[2])
                        # print(r, angle, end=', ')
                        # R.append(r)
                        # A.append(angle)
                        W.append(float(list2_[3]))
                        # x = r * math.cos(angle) + tar_locs[int(tar)-1][0]
                        # y = r * math.sin(angle) + tar_locs[int(tar)-1][1]
                        x = x + tar_locs[int(tar)-1][0] - 1
                        y = y + tar_locs[int(tar)-1][1] - 1

                        circle = plt.Circle((x, y), radius=0.06, lw= 0.01,edgecolor='w',facecolor='b',alpha=1)
                        ax.add_artist(circle)
                        X.append(x)
                        Y.append(y)
                allX.append(X)
                allY.append(Y)
                # colors = np.random.rand(len(allX))
                # plt.scatter(allX, allY,s=1, c=colors, cmap='jet')
        '''plot inner path'''
        filename = '/home/cai/Dropbox/Box_Research/Github/RRARP_LinKern/dat/pyFunc/optimalpath_' + str(count) + '.txt'
        fileopath_ = open(filename, 'r')
        pathx = []
        pathy = []
        seq = [0]
        for i in range(2*nb_tars_+2):
            line_ = fileopath_.readline()
            if i % 2 == 1:
                seq.append(int(line_.split('\t')[0]))
            pathx.append(float(line_.split('\t')[1]))
            pathy.append(float(line_.split('\t')[2]))
        for i in range(len(pathx)):
            if i % 2 == 0:
                plt.plot(pathx[i:i+2], pathy[i:i+2],color='k',alpha=1,linewidth=1)
        # print(pathx)
        # print(pathy)
        # print([pathx[-2],pathx[-1]])
        # print([pathy[-2],pathy[-1]])
        plt.plot([pathx[-2],pathx[-1]], [pathy[-2],pathy[-1]],'k',alpha=1,  linewidth=1)

        # print(seq)
        seq.pop(0)
        seq.pop()
        for v in seq:
            line_ = fileopath_.readline()
            # print(re.split(',|\n', line_))
            innerpath = [int(e) for e in re.split(',|\n', line_) if e != '']
            # print()
            # print(list(itemgetter(*innerpath)(allY[i])))
            plt.plot(list(itemgetter(*innerpath)(allX[v-1])), list(itemgetter(*innerpath)(allY[v-1])), 'k', alpha=1, linewidth=1)

        fileopath_.close()

    file_.close()
    ax.set(xlim=(minx,maxx), ylim=(miny,maxy))
    ax.set_aspect('equal', adjustable='box')
    # plt.xticks(np.arange(minx,maxx,1))
    # plt.yticks(np.arange(miny,maxy,1))
    ax.axis('off')
    # plt.grid(alpha=.5)
    plt.savefig('/home/cai/Dropbox/Box_Research/Github/RRARP_LinKern/dat/pyFunc/temp/t_' + str(count) + '.jpg', format='jpg',dpi=500)
    # plt.show()

################################################################################################

#  transformation files
# instancefile = '../ret/inst_n_10/n_10_h_6.txt'
# instancefile = '../ret/cluster_n_30/n_30_c_5_1.txt'
# instancefile = '/home/cai/Dropbox/Box_Research/Github/RRARP_BD/BendersDecomp/RRARP-BD/dat/test_n_6.txt'
# instancefile = '/home/caigao/Dropbox/Box_Research/Github/RRARP_BD/BendersDecomp/RRARP-BD/dat/test_n_6.txt'

if __name__ == "__main__":
    instancefile = argv[1]
    instsize = int(argv[2])
    minidx = int(argv[3])
    maxidx = int(argv[4])
    for i in range(minidx, maxidx+1):
        plot_instance(instancefile, instsize, i)
