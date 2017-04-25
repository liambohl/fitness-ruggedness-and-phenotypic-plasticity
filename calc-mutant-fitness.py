###############################################################################
# Generate counts of beneficial, deleterious, and neutral one-step mutations
# from final dominant genomes of each test run.
###############################################################################

import glob
import os
import subprocess

instruction_set_basic = "abcdefghijklmnopqrstu"
instruction_set_sense = "abcdefghijklmnopqrstuvw"

def dict_add(d, x):
	'''
	Add x to dictionary d.
	'''
	if x in d:
		d[x] += 1
	else:
		d[x] = 1

def evaluate_genomes(treatment, run, final_dom_gen, mutant_dict, instruction_set):
	'''
	Evaluate the phenotype of base organism and mutant dict. Output summary to summary.txt.
	treatment: name of this treatment (string)
	run: index of this replicate (int)
	final_dom_gen: genome of base organism (string)
	mutant_dict: all 1-step mutants of base organism (dict of strings)
	instruction_set: name of instruction set. Used to determine which instruction set file to specify in Avida options. (string)
	'''
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
		analyze(instruction_set, -1, final_dom_gen, run_name, env)
		base_fitness, base_pnand, base_pnot = get_phenotype(run_name, env, -1)
		
		# Get phenotype for each mutant
		for i, org in enumerate(mutant_dict):
			analyze(instruction_set, i, org, run_name, env[0], env[1])
			fitness, pnand, pnot = get_phenotype(run_name, env, i)
			fitness_list.append(fitness)
			nand_count += pnand
			not_count + pnot
		
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
		out_filename = "../mutant-fitness/{}/env_nand_{}_not_{}/summary.txt".format(run_name, env[0], env[1])
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

def analyze(instruction_set, i, org, run_name, env):
	'''
	Run the given organism in the given environment in Avida analyze mode
	instruction_set: name of instruction set. Used to determine which instruction set file to specify in Avida options. (string)
	i: index of this organism (int)
	org: genome of this organism (string)
	run_name: treatment and replicate number (string)
	env: tuple of reaction values for each task (tuple of ints)
	'''
	# Set name of instruction set file
	if instruction_set == instruction_set_basic:
		inst_set = "instset-heads.cfg"
	else:
		inst_set = "instset-heads-sense.cfg"
	
	# Generate properly configured analyze.cfg file by replacing "%" with arguments
	arg_tuple = (org, env[0], env[1])
	with open("analyze-mutant-temp.cfg", "r") as sample_file, open("analyze-mutant-current.cfg", "w") as analyze_file:
		i = 0
		for line in sample_file:
			if "%" in line:
				analyze_file.write(line.replace("%", str(arg_tuple[i])))
				i += 1
			else:
				analyze_file.write(line)
	
	# Run Avida in analyze mode
	subprocess.call("./avida -a -set ANALYZE_FILE analyze-mutant-current.cfg -def INST_SET " + inst_set + " -set EVENT_FILE events-static.cfg -set VERBOSITY 0", shell = True)
	subprocess.call("mv data/dat ../mutant-fitness/{}/env_nand_{}_not_{}/{}.dat".format(run_name, env[0], env[1], i), shell = True)

def get_phenotype(run_name, environment, n):
	'''
	Get important aspects of the phenotype of an organism from its dat file.
	run_name: treatment and replicate number (string)
	environment: tuple of reaction values for each task (tuple of ints)
	returns:
		fitness (float)
		pnand (bool)
		pnot (bool)
	'''
	with open("../mutant-fitness/{}/env_nand_{}_not_{}/{}.dat".format(run_name, environment[0], environment[1], n), "r") as dat_file:
		for i in range(12):
			dat_file.readline()
		dat_line = dat_file.readline().split()
		fitness, length, seq, gestation, efficiency, pnand, pnot = dat_line
		return float(fitness), bool(pnand), bool(pnot)

def generate_mutants(genome_str, instruction_set):
	'''
	Generate all 1-step mutants of given genome, including point substitutions, deletions, and insertions.
	genome_str: base organism for generating mutants (string with each instruction as one character)
	instruction_set: set of characters used to represent instructions (string)
	return: dict of mutants and counts (dict)
	'''
	# Generate list of one-step mutant genomes
	mutant_dict = {}

	for i in range(len(genome_str)):
		for inst in instruction_set:
			if inst != genome_str[i]:
				dict_add(mutant_dict, genome_str[:i] + inst + genome_str[i + 1:])	# Point
			dict_add(mutant_dict, genome_str[:i] + genome_str[i + 1:])				# Deletion
			dict_add(mutant_dict, genome_str[:i] + inst + genome_str[i:])			# Insertion

	# Handle insertions at end of genome
	for inst in instruction_set:
		dict_add(mutant_dict, genome_str + inst)

	return mutant_dict
	
def main(treatments_list, n_runs, tasks_list):
	'''
	Cycle through all replicates of all treatments to find final dominant org, generate all 1-step mutants from these orgs,
	and calculate phenotypic match score for base org and mutants
	treatments_list: names of experimental and control treatments (list of strings)
	n_runs: number of replicates for each treatment (int)
	tasks_list: names of tasks rewarded and punished (list of strings)
	'''
	for treatment in treatments_list:
		# Choose instruction set
		if treatment == "Plastic":
			instruction_set = instruction_set_sense
		else:
			instruction_set = instruction_set_basic
		
		for run in range(1, n_runs + 1):
			# Get genome
			filename = glob.glob("../analysis/{}_{}/env_nand_1_*/final_dominant.dat".format(treatment, run))[0] # dat file with string genome of final dominant
			with open(filename, "r") as genome_file:
				for i in range(16):
					genome_file.readline()
				genome_str = genome_file.readline().split()[4] # Genome as string of chars is 4th item on line.
			
			# Generate mutants and evaluate their phenotypes for each environment
			mutant_dict = generate_mutants(genome_str, instruction_set)
			evaluate_genomes(treatment, run, genome_str, mutant_dict, instruction_set)

# All 1- or 2-input binary logic operations
all_tasks = ["NOT", "NAND", "AND", "ORN", "OR", "ANDN", "NOR", "XOR", "EQU"]

# Run with three treatments and 10 runs of each
treatments_list = ["Static", "Changing", "Plastic"]
tasks_list = all_tasks[:4]
main(treatments_list, 10, tasks_list)
