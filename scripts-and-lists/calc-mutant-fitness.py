###############################################################################
# Generate counts of beneficial, deleterious, and neutral one-step mutations
# from final dominant genomes of each test run.
###############################################################################

import glob
import os
import subprocess

# instruction_set_basic = ["if-label", "mov-head", "jmp-head", "get-head", "set-flow", "shift-r", "shift-l", "inc", "dec", "push", "pop", "swap-stk", "swap", "add", "sub", "nand", "h-copy", "h-alloc", "h-divide", "IO", "h-search"]
# instruction_set_sense = instruction_set_basic + ["sense-react-NAND", "sense-react-NOT"]

instruction_set_basic = "abcdefghijklmnopqrstu"
instruction_set_sense = "abcdefghijklmnopqrstuvw"

def count_mutants(path, treatment, run, instruction_set):
	# Get genome
	filename = glob.glob(path + "env_nand_1_*/final_dominant.dat")[0] # dat file with string genome of final dominant
	with open(filename, "r") as genome_file:
		for i in range(16):
			genome_file.readline()
		genome_str = genome_file.readline().split()[4] # Genome as string of chars is 4th item on line.
	mutant_dict = generate_mutants(genome_str, instruction_set)
	
	# Generate .org files for each mutant
	# print_orgs(mutant_set, path)
	
	evaluate_genomes(treatment, run, genome_str, mutant_dict, instruction_set)

def evaluate_genomes(treatment, run, final_dom_gen, mutant_dict, instruction_set):
	# Name of this run
	run_name = treatment + "_" + str(run)
	
	# Find value of nand and not for each treatment
	if treatment == "Static":
		environments = ((1, 1),)
	else:
		environments = ((1, -1), (-1, 1))
	
	for env in environments:
		# Summary statistics for mutants of one final dominant in one environment
		fitness_list = []	# List of all fitness values
		nand_count = 0		# Number of orgs to perform nand
		not_count = 0		# Number of orgs to perform not
		
		# Get phenotype for base organism using -1 as index.
		analyze(instruction_set, -1, final_dom_gen, run_name, env[0], env[1])
		base_fitness, base_pnand, base_pnot = get_phenotype(run_name, env, -1)
		
		# Get phenotype for each mutant
		for i, org in enumerate(mutant_dict):
			analyze(instruction_set, i, org, run_name, env[0], env[1])
			fitness, pnand, pnot = get_phenotype(run_name, env, i)
			fitness_list.append(fitness)
			nand_count += bool(int(pnand))
			not_count += bool(int(pnot))
		
		# Count number of beneficial, deleterious mutations
		total_mutations = len(fitness_list)
		del_mutations, ben_mutations, neu_mutations = 0, 0, 0
		for fitness in fitness_list:
			if fitness > base_fitness:
				ben_mutations += 1
			elif fitness < base_fitness:
				del_mutations += 1
			else:
				neu_mutations += 1
		
		# Output to file
		out_filename = "data/analysis/{}/env_nand_{}_not_{}/summary.txt".format(run_name, env[0], env[1])
		with open(out_filename, "w") as summary_file:
			# Write base organism phenotype
			summary_file.write("Base organism performed NAND: {}; NOT: {}\n\n".format(base_pnand, base_pnot))
			
			# Write mutant summary stats
			summary_file.write("Out of {} 1-step mutations:\n".format(total_mutations))

			stats_tuple = (
				(del_mutations, "deleterious"),
				(neu_mutations, "neutral"),
				(ben_mutations, "beneficial"),
				(nand_count, "performed NAND"),
				(not_count, "performed NOT")
			)
			for stat in stats_tuple:
				summary_file.write("{:5d} ({:7.2%}) {:s}\n".format(stat[0], round(stat[0] / total_mutations, 4), stat[1]))

def analyze(instruction_set, *argv): #i, org, run_name, nand_val, not_val):
	# Set name of instruction set file
	if instruction_set == instruction_set_basic:
		inst_set = "instset-heads.cfg"
	else:
		inst_set = "instset-heads-sense.cfg"
	
	# Generate properly configured analyze.cfg file
	with open("analyze-mutant-temp.cfg", "r") as sample_file, open("analyze-mutant-current.cfg", "w") as analyze_file:
		i = 0
		for line in sample_file:
			if "%" in line:
				analyze_file.write(line.replace("%", str(argv[i])))
				i += 1
			else:
				analyze_file.write(line)
	
	# Run Avida in analyze mode
	subprocess.call("./avida -a -set ANALYZE_FILE analyze-mutant-current.cfg -def INST_SET " + inst_set + " -set EVENT_FILE events-static.cfg -set VERBOSITY 0", shell = True)

def get_phenotype(run_name, environment, n):
	with open("data/analysis/{:s}/env_nand_{:d}_not_{:d}/{:d}.dat".format(run_name, environment[0], environment[1], n), "r") as dat_file:
		for i in range(12):
			dat_file.readline()
		dat_line = dat_file.readline().split()
		fitness, length, seq, gestation, efficiency, pnand, pnot = dat_line
		return fitness, pnand, pnot

#def print_orgs(orgs_set, path):
	#if not os.path.exists(path + "mutants"):
		#os.makedirs(path + "mutants")
	#for i, genome in enumerate(mutant_set):
		#with open(path + "mutants/" + str(i) + ".gen", "w") as genome_file:
			#for inst in genome:
				#genome_file.write(inst)
	## Count total mutants
	#with open(path + "mutants-count.txt", "x") as mutants_file:
		#mutants_file.write(str(len(mutant_set)))

def generate_mutants(genome_str, instruction_set):
	# Generate list of one-step mutant genomes
	mutant_dict = {}

	for i in range(len(genome_str)):
		for inst in instruction_set:
			if inst != genome_str[i]:
				dict_add(mutant_dict, genome_str[:i] + inst + genome_str[i + 1:])	# Point
			dict_add(mutant_dict, genome_str[:i] + genome_str[i + 1:])			# Deletion
			dict_add(mutant_dict, genome_str[:i] + inst + genome_str[i:])			# Insertion

	# Handle insertions at end of genome
	for inst in instruction_set:
		dict_add(mutant_dict, genome_str + inst)

	return mutant_dict
	
def dict_add(d, x):
	if x in d:
		d[x] += 1
	else:
		d[x] = 1

treatments_list = ["Static", "Two-envs", "Two-envs-sense"]
n_runs = 10

# For each test run
for treatment in treatments_list:
	# Choose instruction set
	if treatment == "Two-envs-sense":
		instruction_set = instruction_set_sense
	else:
		instruction_set = instruction_set_basic
	
	for run in range(1, n_runs + 1):
		path = "../analysis/" + treatment + "_" + str(run) + "/final_dom/" # Path to directory for this treatment
		count_mutants(path, treatment, run, instruction_set)

