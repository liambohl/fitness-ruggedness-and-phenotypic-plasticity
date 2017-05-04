###############################################################################
# Iterate over all treatments and all replicates to find the 1-step mutants
# of the final dominant genotype. Evaluate the phenotypic match score of the
# final dominant genotypes and mutants across all environments they would
# encounter.
# Summarize mean and standard deviation of these scores for each replicate.
###############################################################################

import os
import subprocess
import numpy
import multiprocessing

instruction_set_basic = "abcdefghijklmnopqrstuvwxyz"
instruction_set_sense = "abcdefghijklmnopqrstuvwxyzAB"

def evaluate_genomes(treatment, run, final_dom_gen, mutant_list, instruction_set):
	'''
	Evaluate the phenotype of base organism and mutant list. Output summary to summary.txt.
	treatment: name of this treatment (string)
	run: index of this replicate (int)
	final_dom_gen: genome of base organism (string)
	mutant_list: all 1-step mutants of base organism (list of strings)
	instruction_set: name of instruction set. Used to determine which instruction set file to specify in Avida options. (string)
	'''
	# Name of this run
	run_name = treatment + "_" + str(run)
	
	# Find value of nand and not for each treatment
	if treatment == "Static":
		environments = ((1, 1),)
	else:
		environments = ((1, 1) ,(1, -1), (-1, 1), (-1, -1))
	
	# Get phenotypic match scores for base organism and all mutants.
	max_score = len(environments * len(environments[0]))
	base_score = score_org(instruction_set, -1, final_dom_gen, run_name, environments)
	score_list = []	# List of phenotypic match scores for mutants
	for i, org in enumerate(mutant_list):
		score = score_org(instruction_set, i, org, run_name, environments)
		score_list.append(score)
	
	# Summary Stats
	total_mutations = len(score_list)
	mean_score = numpy.mean(score_list)
	stdev_score = numpy.std(score_list)
	
	# Count number of beneficial, deleterious mutations
	del_mutations, ben_mutations, neu_mutations = 0, 0, 0
	for score in score_list:
		if score > base_score:
			ben_mutations += 1
		elif score < base_score:
			del_mutations += 1
		else:
			neu_mutations += 1
		
	# Output to file
	out_filename = "../mutant-fitness/{}/summary.txt".format(run_name)
	with open(out_filename, "w") as summary_file:
		# Write base organism phenotype
		summary_file.write("Max phenotypic match score: {:d}\n".format(max_score))
		summary_file.write("Base organism phenotypic match score: {:d}\n\n".format(base_score))
		
		# Write mutant summary stats
		summary_file.write("Out of {} 1-step mutations:\n".format(total_mutations))

		stats_tuple = (
			(del_mutations, "deleterious"),
			(neu_mutations, "neutral"),
			(ben_mutations, "beneficial"),
		)
		for stat in stats_tuple:
			summary_file.write("{:5d} ({:7.2%}) {:s}\n".format(stat[0], round(stat[0] / total_mutations, 4), stat[1]))
		summary_file.write("Mean score of 1-step mutants: {:.3f}\n".format(mean_score))
		summary_file.write("Standard dev score of 1-step mutants: {:.3f}\n\n".format(stdev_score))
		
		# Write all match scores
		score_lines = '\n'.join([str(score) for score in score_list])
		summary_file.write(score_lines)

def score_org(instruction_set, org, gen, run_name, environments):
	'''
	Run the given organism in each environment in Avida analyze mode. Score the match between this organism's behavior and the set of environments.
	instruction_set: name of instruction set. Used to determine which instruction set file to specify in Avida options. (string)
	i: index of this organism (int)
	org: genome of this organism (string)
	run_name: treatment and replicate number (string)
	environments: tuple of environments; each environment is a tuple of reaction values for each task (tuple of tuples of ints)
	returns: phenotypic match score - sum over all environments:
		+1 for tasks performed and rewarded
		+1 for tasks not performed and punished
		-1 for tasks performed and punished
		-1 for tasks not performed and rewarded
	'''
	# Set name of instruction set file
	if instruction_set == instruction_set_basic:
		inst_set = "instset-heads.cfg"
	else:
		inst_set = "instset-heads-sense.cfg"
	
	p_match_score = 0
	for env in environments:
		# Generate properly configured analyze.cfg file by replacing "%" with arguments
		arg_tuple = (gen, env[0], env[1])
		with open("analyze-mutant-temp.cfg", "r") as sample_file, open("analyze-mutant-current.cfg", "w") as analyze_file:
			i = 0
			for line in sample_file:
				if "%" in line:
					analyze_file.write(line.replace("%", str(arg_tuple[i])))
					i += 1
				else:
					analyze_file.write(line)
		
		# Run Avida in analyze mode
		print("Return code for {}/{}_{}/{}: ".format(run_name, env[0], env[1], org))
		print(subprocess.call("./avida -a -set ANALYZE_FILE analyze-mutant-current.cfg -def INST_SET " + inst_set + " -set EVENT_FILE events-static.cfg -set VERBOSITY 0", shell = True))
		
		# Save analyze.cfg and dat files
		# subprocess.call("cp analyze-mutant-current.cfg ../scoring/{}/{}_{}/{}.cfg".format(run_name, env[0], env[1], org), shell = True)
		# subprocess.call("cp data/dat ../scoring/{}/{}_{}/{}.dat".format(run_name, env[0], env[1], org), shell = True)
		
		# Score phenotypic match in this environment
		with open("data/dat", "r") as dat_file:
			for i in range(8):
				dat_file.readline()
			dat_line = dat_file.readline().split()[-2:]
			pnand, pnot = dat_line
			phenotype = tuple([-1 if i == "0" else 1 for i in (pnand, pnot)]) # 1 for tasks performed; -1 for tasks not performed
			for task, performed in enumerate(phenotype):
				if performed == env[task]: # "task x was performed" equals "task x was rewarded"
					p_match_score += 1
				else:
					p_match_score -= 1 # This is redundant, but keeps the range of possible scores centered at 0.
	return p_match_score

def generate_mutants(genome_str, instruction_set):
	'''
	Generate all 1-step mutants of given genome, including point substitutions, deletions, and insertions.
	genome_str: base organism for generating mutants (string with each instruction as one character)
	instruction_set: set of characters used to represent instructions (string)
	return: list of mutants (list)
	'''
	# Generate list of one-step mutant genomes
	mutant_list = []

	for i in range(len(genome_str)):
		for inst in instruction_set:
			if inst != genome_str[i]:
				mutant_list.append(genome_str[:i] + inst + genome_str[i + 1:])	# Point
			mutant_list.append(genome_str[:i] + inst + genome_str[i:])			# Insertion
		mutant_list.append(genome_str[:i] + genome_str[i + 1:])					# Deletion

	# Handle insertions at end of genome
	for inst in instruction_set:
		mutant_list.append(genome_str + inst)

	return mutant_list
	
def main(treatments_list, n_runs, tasks_list):
	'''
	Cycle through all replicates of all treatments to find final dominant org, generate all 1-step mutants from these orgs,
	and calculate phenotypic match score for base org and mutants
	treatments_list: names of experimental and control treatments (list of strings)
	tasks_list: names of tasks rewarded and punished (list of strings)
	'''
	
	# One process to evaluate base + mutant genotypes of each run
	process_list = []
		
	for start, treatment in enumerate(treatments_list):
		# Choose instruction set
		if treatment == "Plastic":
			instruction_set = instruction_set_sense
		else:
			instruction_set = instruction_set_basic
		
		for run in range(100 * start, 100 * start + n_runs):
			# Get genome
			filename = "../analysis/{}_{}/final_dom.dat".format(treatment, run)
			with open(filename, "r") as genome_file:
				for i in range(6):
					genome_file.readline()
				genome_str = genome_file.readline().strip()
			
			# Save base genome
			base_filename = "../mutant-fitness/{}_{}/base-genome".format(treatment, run)
			with open(base_filename, "w") as genome_file:
				genome_file.write(genome_str)
			
			# Generate mutants and evaluate their phenotypes for each environment
			mutant_list = generate_mutants(genome_str, instruction_set)
			
			# Output all mutant genomes to file
			mutant_filename = "../mutant-fitness/{}_{}/mutant-genomes".format(treatment, run)
			with open(mutant_filename, "w") as mutant_file:
				for i, gen in enumerate(mutant_list):
					mutant_file.write("{}\n".format(gen))
			
			# Begin processing
			process_list.append(subprocess.Process(target = evaluate_genomes, args = (treatment, run, genome_str, mutant_list, instruction_set)))
			process_list[-1].start()
	
	# Finish processing
	for process in process_list:
		process.join()

# All 1- or 2-input binary logic operations
all_tasks = ["NOT", "NAND", "AND", "ORN", "OR", "ANDN", "NOR", "XOR", "EQU"]

# Run with three treatments and 10 runs of each
treatments_list = ["Static", "Changing", "Plastic"]
n_runs = 10
tasks_list = all_tasks[:2]
main(treatments_list, n_runs, tasks_list)
