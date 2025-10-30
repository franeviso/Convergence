# mutual/config.ini
# simulation configuration file

N = 10000								# number of individuals per guild
E = 10000								# number of attempted mating events per generation 
G = 10000								# generations

update = "synchronous"					# "synchronous" or "asynchronous"

# define the spatial distance range
# range_type == "topological": range_value = number of 'neighbors' per guild
#							   distance range <- sqrt(range_value * N / pi)
# range_type == "euclidean":   distance range <- range_value
# range_type == "homogeneous": distance range N/A
range_type = "topological"				# distance metric "topological", "euclidean" or "homogeneous"
range_value = 283						# range value
spacing = "random"						# spacing, "random" or "grid"

animal.phenotype_dist = c(-1, 0)		# tolerable mutual phenotype distance
animal.threshold = 0.75					# threshold genetic compatibility
animal.mutation_rate = 0.001			# mutation rate per bit

animal.init.type = "random"				# "random", " strict" or "clone"
animal.init.p_bits = 50					# average initial 'phenotype' bits set
animal.init.c_bits = 90					# average initial 'compat' bits set
animal.init.n_bits = 128				# average initial 'neutral' bits set

plant.phenotype_dist = c(0, 1)			# tolerable mutual phenotype distance
plant.threshold = 0.75					# threshold genetic compatibility
plant.mutation_rate = 0.001				# mutation rate per bit

plant.init.type = "random"				# "random", "strict" or "clone"
plant.init.p_bits = 20					# average initial 'phenotype' bits set
plant.init.c_bits = 0					# average initial 'compat' bits set
plant.init.n_bits = 0					# average initial 'neutral' bits set

# data output
data.hist_interval = 10					# [generations]; 0: off
data.image_interval = 0				# [generation] deprecated, 0: off
