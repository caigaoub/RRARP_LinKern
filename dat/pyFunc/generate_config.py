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


# instsize_min = 6
# instsize_max = 10
# diff_levels = ['e','m','h']
# subidx_min = 1 
# subidx_max = 10


# path = '/home/cai/Dropbox/Box_Research/Github/RRARP_LinKern/dat'
# idx_config = 1
# # for instsize in range(instsize_min, instsize_max+1):
# # 	for sidx in range(subidx_min, subidx_max+1):
# # 		print(idx_config)
# # 		file = open(path + "/configs/config2_" + str(idx_config), "w")
# # 		file.write("n_"+str(instsize) + '_' + 'e' + '_' + str(sidx)+".dat")
# # 		file.close()
# # 		idx_config += 1
# # 	print("done with EASY instance of", instsize, "targets")

# # for instsize in range(instsize_min, instsize_max+1):
# # 	for sidx in range(subidx_min, subidx_max+1):
# # 		print(idx_config)
# # 		file = open(path + "/configs/config2_" + str(idx_config), "w")
# # 		file.write("n_"+str(instsize) + '_' + 'm' + '_' + str(sidx)+".dat")
# # 		file.close()
# # 		idx_config += 1
# # 	print("done with MEDIUM instance of", instsize, "targets")


# for instsize in range(instsize_min, instsize_max+1):
# 	for sidx in range(subidx_min, subidx_max+1):
# 		print(idx_config)
# 		file = open(path + "/configs/to_instances/cg_" + str(idx_config), "w")
# 		file.write("n_"+str(instsize) + '_' + 'h' + '_' + str(sidx)+".dat")
# 		file.close()
# 		idx_config += 1
# 	print("done with HARD instance of", instsize, "targets")










'''-------------------------------------------------------------------------------------------'''
instsize_min = 6
instsize_max = 10
diff_levels = ['e','m','h']
subidx_min = 1 
subidx_max = 10


path = '/home/cai/Dropbox/Box_Research/Github/RRARP_LinKern/dat'
idx_config = 1
graphsize = [30, 40, 50, 60, 80, 100, 120, 150, 180, 200]
alpha = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
for size in graphsize:
	for a in alpha:
		print(idx_config)
		file = open(path + "/configs/pulse/cg_" + str(idx_config), "w")
		# file.write("n_6_h_1.graph_"+str(size) + '\t' + str(a))
		file.write(str(size) + '\t' + str(a))

		file.close()
		idx_config += 1
