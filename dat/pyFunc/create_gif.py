import imageio
import os
import fnmatch
dirpath = "/home/cai/Dropbox/Box_Research/Github/RRARP_LinKern/dat/pyFunc/temp/"
nb_jpg = len(fnmatch.filter(os.listdir(dirpath), '*.jpg'))
with imageio.get_writer(dirpath+ 'lk.gif', mode='I', duration = 0.5) as writer:
    for i in range(0, nb_jpg):
    # for i in range(143, 172):
        image = imageio.imread(dirpath + "t_"+str(i)+".jpg")
        writer.append_data(image)
        # os.remove("temp_"+str(i)+".jpg") # delete image