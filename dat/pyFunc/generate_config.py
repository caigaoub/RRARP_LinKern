import numpy as np

"""  generate configuration files """

# instsize_min = 11
# instsize_max = 210
# diff_levels = ['e','m','h']
# subidx_min = 1 
# subidx_max = 10


# path = '/home/cai/Dropbox/Box_Research/Github/RRARP_LinKern/ret/'
# idx_config = 1
# for instsize in range(instsize_min, instsize_max+1):
# 	for sidx in range(subidx_min, subidx_max+1):
# 		print(idx_config)
# 		file = open(path + "/configs/config_" + str(idx_config), "w")
# 		file.write("n_"+str(instsize) + '_' + 'e' + '_' + str(sidx)+".dat")
# 		file.close()
# 		idx_config += 1
# 	print("done with EASY instance of", instsize, "targets")

# for instsize in range(instsize_min, instsize_max+1):
# 	for sidx in range(subidx_min, subidx_max+1):
# 		print(idx_config)
# 		file = open(path + "/configs/config_" + str(idx_config), "w")
# 		file.write("n_"+str(instsize) + '_' + 'm' + '_' + str(sidx)+".dat")
# 		file.close()
# 		idx_config += 1
# 	print("done with MEDIUM instance of", instsize, "targets")


# for instsize in range(instsize_min, instsize_max+1):
# 	for sidx in range(subidx_min, subidx_max+1):
# 		print(idx_config)
# 		file = open(path + "/configs/config_" + str(idx_config), "w")
# 		file.write("n_"+str(instsize) + '_' + 'h' + '_' + str(sidx)+".dat")
# 		file.close()
# 		idx_config += 1
# 	print("done with HARD instance of", instsize, "targets")


instsize_min = 6
instsize_max = 10
diff_levels = ['e','m','h']
subidx_min = 1 
subidx_max = 10


path = '/home/cai/Dropbox/Box_Research/Github/RRARP_LinKern/dat'
idx_config = 1
# for instsize in range(instsize_min, instsize_max+1):
# 	for sidx in range(subidx_min, subidx_max+1):
# 		print(idx_config)
# 		file = open(path + "/configs/config2_" + str(idx_config), "w")
# 		file.write("n_"+str(instsize) + '_' + 'e' + '_' + str(sidx)+".dat")
# 		file.close()
# 		idx_config += 1
# 	print("done with EASY instance of", instsize, "targets")

# for instsize in range(instsize_min, instsize_max+1):
# 	for sidx in range(subidx_min, subidx_max+1):
# 		print(idx_config)
# 		file = open(path + "/configs/config2_" + str(idx_config), "w")
# 		file.write("n_"+str(instsize) + '_' + 'm' + '_' + str(sidx)+".dat")
# 		file.close()
# 		idx_config += 1
# 	print("done with MEDIUM instance of", instsize, "targets")


for instsize in range(instsize_min, instsize_max+1):
	for sidx in range(subidx_min, subidx_max+1):
		print(idx_config)
		file = open(path + "/configs/config2_" + str(idx_config), "w")
		file.write("n_"+str(instsize) + '_' + 'h' + '_' + str(sidx)+".dat")
		file.close()
		idx_config += 1
	print("done with HARD instance of", instsize, "targets")