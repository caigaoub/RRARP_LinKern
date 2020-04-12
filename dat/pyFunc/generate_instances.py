import numpy as np
import math
import random
import matplotlib.pyplot as plt
import matplotlib.patches as pch
from decimal import *

def generate_instance(nb_targets, radius, densitycoef):
    precision = 1000.0
    # densitycoef = 0.5
    area = (densitycoef**2) * nb_targets * radius**2 * math.pi
    width = math.sqrt(area)
    height = width
    tar_locs = []
    rewards_ratio = []
    #depot 1 and 2
    while True:
        x1 =  round(random.uniform(0, width)*precision)/precision
        y1 =  round(random.uniform(0, height)*precision)/precision
        x2 =  round(random.uniform(0, width)*precision)/precision
        y2 =  round(random.uniform(0, height)*precision)/precision
        if x1 != x2 or y1 != y2:
            tar_locs.append((x1,y1,0))
            tar_locs.append((x2,y2,0))
            break
    # target location 
    tidx = 1
    itr = 0
    while True:
        x =  round(random.uniform(0, width)*precision)/precision
        y =  round(random.uniform(0, height)*precision)/precision
        # print(x, y)
        itr += 1
        # if itr > 10000:
        #     print('endless loop')
        #     exit()
        covered = False
        for e in tar_locs:
            dist = math.sqrt((e[0]-x)**2 + (e[1]-y)**2) - e[2] - radius - 0.5
            if dist < 0:
                covered = True
                break
        if covered == False:
            tar_locs.append((x,y,radius))
            tidx += 1
        if tidx > nb_targets:
            break
    # print(tar_locs)
    rewards_ratio = []
    for i in range(0, nb_targets):
        rewards_ratio.append(round(random.uniform(0.2, 0.3)*precision)/precision)
    return tar_locs, rewards_ratio, width, height


def write_instance(tar_locs, rewards_ratio, filename):
    nb_obps = 200

    file = open(filename, "w")
    file.write(str(len(tar_locs)-2) + '\t' + str(nb_obps) + '\n')
    file.write(str(tar_locs[0][0]) + '\t' + str(tar_locs[0][1]) + '\n')
    file.write(str(tar_locs[1][0]) + '\t' + str(tar_locs[1][1]) + '\n')
    for i in range(2,len(tar_locs)):
        file.write(str(tar_locs[i][0]) + '\t' + str(tar_locs[i][1]) + '\t' + str(tar_locs[i][2]) + '\n')
    # print(rewards_ratio)
    for i in range(0,len(rewards_ratio)-1):
        file.write(str(rewards_ratio[i]) + '\t')
    file.write(str(rewards_ratio[-1]) + '\n')

    
    ''' write observation points '''
    K = 16
    precision = 1000.0
    for i in range(0, len(tar_locs)-2):
        '''write K boundary points '''
        for j in range(0, K):
            r = 1.0
            angle = 2.0*np.pi/K * j;
            rewp = 0
            file.write(str(r) + ':' + str(angle) + ':' + str(rewp) + '\t')
        # file.write('T' + str(i+1) + '\t')
        for j in range(0, nb_obps-1):
            r = round(math.sqrt(np.random.uniform(0.1,0.9))*precision)/precision
            angle = round(np.random.uniform(0,2.0*np.pi)*precision)/precision
            rewp = round(np.random.uniform(0.1,0.2)*precision)/precision
            file.write(str(r) + ':' + str(angle) + ':' + str(rewp) + '\t')
        r = round(math.sqrt(np.random.uniform(0,1))*precision)/precision
        angle = round(np.random.uniform(0,2.0*np.pi)*precision)/precision
        rewp = round(np.random.uniform(0.1,0.2)*precision)/precision
        file.write(str(r) + ':' + str(angle) + ':' + str(rewp) + '\n')
    
    file.close()

if __name__ == "__main__":
    ''' basic setting '''
    if False:
        radius = 1.0
        size_range = [6, 200]
        path = '/home/cai/Dropbox/Box_Research/Github/RRARP_LinKern/dat/'
        for nb_targets in range(size_range[0], size_range[1]+1):
            # densitycoef = 5
            print("currently generate instance with size ", nb_targets)
            # for i in range(1,11):
            #     print('e',i)
            #     tar_locs, rewards_ratio, width, height = generate_instance(nb_targets, radius, densitycoef)
            #     write_instance(tar_locs, rewards_ratio, path +'n_'+ str(nb_targets) +'_' + 'e_' + str(i) + '.dat')

            # densitycoef = 4
            # for i in range(1,11):
            #     print('m',i)

            #     tar_locs, rewards_ratio, width, height = generate_instance(nb_targets, radius, densitycoef)
            #     write_instance(tar_locs, rewards_ratio,  path +'n_' + str(nb_targets) +'_' + 'm_' + str(i) + '.dat')

            densitycoef = 4
            for i in range(1,11):
                print('h',i)
                tar_locs, rewards_ratio, width, height = generate_instance(nb_targets, radius, densitycoef)
                write_instance(tar_locs, rewards_ratio,  path +'n_'+ str(nb_targets) +'_' + 'h_' + str(i) + '.dat')





node 0 to sink 8: 0 28 16 45 41 18 8 rewards on it: 0.755
node 1 to sink 8: 1 16 45 41 18 8 rewards on it: 0.565
node 2 to sink 8: 2 24 21 30 18 8 rewards on it: 0.525
node 3 to sink 8: 3 30 18 8 rewards on it: 0.302
node 4 to sink 8: 4 30 18 8 rewards on it: 0.302
node 5 to sink 8: 5 42 18 8 rewards on it: 0.295
node 6 to sink 8: 6 22 40 18 8 rewards on it: 0.5
node 7 to sink 8: 7 18 8 rewards on it: 0.172
node 8 to sink 8: 8 rewards on it: 0
node 9 to sink 8: 9 27 8 rewards on it: 0.165
node 10 to sink 8: 10 33 27 8 rewards on it: 0.344
node 11 to sink 8: 11 33 27 8 rewards on it: 0.344
node 12 to sink 8: 12 44 38 33 27 8 rewards on it: 0.618
node 13 to sink 8: 13 19 38 33 27 8 rewards on it: 0.62
node 14 to sink 8: 14 23 25 35 20 27 8 rewards on it: 0.712
node 15 to sink 8: 15 32 17 25 35 20 27 8 rewards on it: 0.812
node 16 to sink 8: 16 45 41 18 8 rewards on it: 0.565
node 17 to sink 8: 17 25 35 20 27 8 rewards on it: 0.697
node 18 to sink 8: 18 8 rewards on it: 0.172
node 19 to sink 8: 19 38 33 27 8 rewards on it: 0.62
node 20 to sink 8: 20 27 8 rewards on it: 0.322
node 21 to sink 8: 21 30 18 8 rewards on it: 0.417
node 22 to sink 8: 22 40 18 8 rewards on it: 0.5
node 23 to sink 8: 23 25 35 20 27 8 rewards on it: 0.712
node 24 to sink 8: 24 21 30 18 8 rewards on it: 0.525
node 25 to sink 8: 25 35 20 27 8 rewards on it: 0.592
node 26 to sink 8: 26 42 18 8 rewards on it: 0.402
node 27 to sink 8: 27 8 rewards on it: 0.165
node 28 to sink 8: 28 16 45 41 18 8 rewards on it: 0.755
node 29 to sink 8: 29 20 27 8 rewards on it: 0.492
node 30 to sink 8: 30 18 8 rewards on it: 0.302
node 31 to sink 8: 31 36 35 20 27 8 rewards on it: 0.68
node 32 to sink 8: 32 17 25 35 20 27 8 rewards on it: 0.812
node 33 to sink 8: 33 27 8 rewards on it: 0.344
node 34 to sink 8: 34 21 30 18 8 rewards on it: 0.589
node 35 to sink 8: 35 20 27 8 rewards on it: 0.435
node 36 to sink 8: 36 35 20 27 8 rewards on it: 0.555
node 37 to sink 8: 37 25 35 20 27 8 rewards on it: 0.701
node 38 to sink 8: 38 33 27 8 rewards on it: 0.464
node 39 to sink 8: 39 20 27 8 rewards on it: 0.444
node 40 to sink 8: 40 18 8 rewards on it: 0.34
node 41 to sink 8: 41 18 8 rewards on it: 0.306
node 42 to sink 8: 42 18 8 rewards on it: 0.295
node 43 to sink 8: 43 16 45 41 18 8 rewards on it: 0.723
node 44 to sink 8: 44 38 33 27 8 rewards on it: 0.618
node 45 to sink 8: 45 41 18 8 rewards on it: 0.407
