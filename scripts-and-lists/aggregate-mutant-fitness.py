###############################################################################
# Use data in all summary.txt files (from mutant fitness analysis) to find
# aggregate statistics. Calculate effects of changing environments and
# phenotypic plasticity on proportions of deleterious, neutral, and beneficial
# one-step mutations.
###############################################################################

import numpy

# Each tuple of ints is an environment characterized by the value of NAND and NOT, respectively.
static = ((1, 1),)
two_envs = ((1, -1), (-1, 1))

treatment_dict = {"Static": static, "Two-envs": two_envs, "Two-envs-sense": two_envs}
n_runs = 10

#with open("data/analysis/aggregate_summary.txt", "w") as out_file:
	#for treatment in treatment_dict:
		#
		## Header
		#out_file.write("---------- {} ----------\n\n".format(treatment))
		#
		#for env in treatment_dict[treatment]:
			## proportions for each run
			#p_lists = {name: [] for name in ["del", "neu", "ben", "pnand", "pnot"]}
			#
			#for run in range(1, n_runs + 1):
				#summary_filename = "data/analysis/{}_{}/env_nand_{}_not_{}/summary.txt".format(treatment, run, env[0], env[1])
				#with open(summary_filename, "r") as summary_file:
					## Skip first 3 lines
					#for i in range(3):
						#summary_file.readline()
					#for stat in p_lists:
						#line = summary_file.readline()
						#start, stop = line.index("("), line.index("%")
						#p = float(line[start + 1 : stop]) / 100 # Convert percentage as string to decimal
						#p_lists[stat].append(p)
				#
			## Find average and standard deviation for each p
			#mean_p = {stat: numpy.mean(p_lists[stat]) for stat in p_lists}
			#stdev_p = {stat: numpy.std(p_lists[stat]) for stat in p_lists}
			#
			## Write summary stats for treatment to output file
			#out_file.write("environment NAND {} NOT {}\n".format(env[0], env[1]))
			#out_file.write("stat  :  mean  ,  sdev\n")
			#for stat in sorted(mean_p):
				#out_file.write("{:5s} : {:6.2%} , {:6.2%}\n".format(stat, mean_p[stat], stdev_p[stat]))
			#out_file.write("\n")

def getPhenotypeFromDat(treatment, run, env, org):
	try:
		with open("data/analysis/{}_{}/env_nand_{}_not_{}/{}.dat".format(treatment, run, env[0], env[1], org), "r") as org_file:
			for i in range(12):
				org_file.readline()
			line_list = org_file.readline().split()
			fitness = float(line_list[0])
			pnand = bool(int(line_list[-2]))
			pnot = bool(int(line_list[-1]))
	except FileNotFoundError:
		return -1, -1, -1
	return fitness, pnand, pnot

n_tiers = 2
with open("data/analysis/aggregate_summary2.txt", "w") as out_file:
	for treatment in treatment_dict:
		
		# Header
		out_file.write("---------- {} ----------\n\n".format(treatment))
		
		for env in treatment_dict[treatment]:
			out_file.write("environment NAND {} NOT {}\n\n".format(env[0], env[1]))
			
			# Lists for storing phenotype info. 0th tier (0th element) is base organisms; nth tier is n-step mutants
			org_count = [0 for i in range(n_tiers)]
			nand_count = [0 for i in range(n_tiers)]
			not_count = [0 for i in range(n_tiers)]
			fitness_list = [[] for i in range(n_tiers)]

			# Aggregate phenotypes for all 1-step mutants
			for run in range(1, n_runs + 1):
				org = -1 # Name of base organism
				tier = 0 # 0th tier mutant
				while True:
					fitness, pnand, pnot = getPhenotypeFromDat(treatment, run, env, org)
					if fitness == -1: # FileNotFoundError
						print("{}:{}:{}:{}".format(treatment, run, env, org))
						break

					org_count[tier] += 1
					nand_count[tier] += pnand
					not_count[tier] += pnot
					fitness_list[tier].append(fitness)

					org += 1
					tier = 1
			mean_fitness = [numpy.mean(i) for i in fitness_list]
			std_fitness = [numpy.std(i) for i in fitness_list]
			p_nand = [count / org_count[i] for i, count in enumerate(nand_count)]
			p_not = [count / org_count[i] for i, count in enumerate(not_count)]
			
			# Write stats to file
			out_file.write("Organisms     :  Base   tier 1\n")
			out_file.write("mean fitness  : {:7.2f} {:7.2f}\n".format(mean_fitness[0], mean_fitness[1]))
			out_file.write("std fitness   : {:7.2f} {:7.2f}\n".format(std_fitness[0], std_fitness[1]))
			out_file.write("performed NAND: {:7.2%} {:7.2%}\n".format(p_nand[0], p_nand[1]))
			out_file.write("performed NOT : {:7.2%} {:7.2%}\n".format(p_not[0], p_not[1]))

			out_file.write("\n")
