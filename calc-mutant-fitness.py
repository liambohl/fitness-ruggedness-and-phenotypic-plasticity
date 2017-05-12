###############################################################################
# Iterate over all treatments and all replicates to find the 1-step mutants
# of the final dominant genotype. Evaluate the phenotypic match score of the
# final dominant genotypes and mutants across all environments they would
# encounter.
# Summarize mean and standard deviation of these scores for each replicate.
###############################################################################

import subprocess
import numpy
import multiprocessing
import random

instruction_alphabet = "abcdefghijklmnopqrstuvwxyzABCDEFGHI"

def skip_legend(f):
	'''
	move to first line of real data in file f.
	'''
	legend = False
	while True:
		line = f.readline()
		if line[:2] == "# ":
			legend = True
		if line == "\n" and legend == True:
			return

def permute_environments(tasks_remaining, environments = None):
	'''
	Expand the given list of environments to include all possibilities of next task.
	environments: all permutations of values for previous tasks (list of tuples)
	tasks_remaining: which tasks are measured in this experiment (list of ints)
	returns: expanded list of environments
	'''
	if environments == None:
		environments = [[]]
	# Expand list by one task
	if tasks_remaining[0] == 0: # Next task is not measured
		new_environments = [env + [0] for env in environments]
	else:
		new_environments = []
		for env in environments:
			new_environments.append(env + [1])
			new_environments.append(env + [-1])

	# Call next step and return
	if len(tasks_remaining) > 1:
		return permute_environments(tasks_remaining[1:], new_environments)
	else:
		return new_environments

def evaluate_genomes(treatment, run, sensing, tasks_list, step_list):
	'''
	Evaluate the phenotype of base organism and mutant list. Output summary to summary.txt.
	treatment: name of this treatment (string)
	run: index of this replicate (int)
	sensing: True iff organisms evolved with sensing. Used to determine which instruction set file to specify in Avida options. (bool)
	tasks_list: which tasks are rewarded, punished, and measured (list of ints)
	step_list: Hamming distances from base organism at which to evaluate mutants (list of ints)
	'''
	# Name of this run
	run_name = treatment + "_" + str(run)
	
	# Build list of environments, find max score
	n_tasks = tasks_list.count(1)
	if treatment == "Static":
		environments = [tasks_list]
		max_score = n_tasks
	else:
		environments = permute_environments(tasks_list)
		max_score = n_tasks * 2 ** n_tasks # Number of measured tasks x number of environments
	
	# Get phenotypic match scores for base organism and all mutants.
	base_score = score_orgs(sensing, run_name, environments, 0)[0]
	mutant_scores = [score_orgs(sensing, run_name, environments, step) for step in step_list]
	
	# Output to file
	out_filename = "../mutant-fitness/{}/summary.txt".format(run_name)
	with open(out_filename, "w") as summary_file:
		# Max score and base score
		summary_file.write("{:12s}: {:d}\n".format("Max score", max_score))
		summary_file.write("{:12s}: {:d}\n".format("Base org", base_score))
		
		for i, step in enumerate(mutant_scores[1:]):
			# Header
			summary_file.write("\n---- Step {} -----\n".format(step_list[i + 1]))
			
			# Summary Stats
			total_mutations = len(step)
			mean_score = numpy.mean(step)
			stdev_score = numpy.std(step)
			summary_file.write("{:12s}: {}\n".format("Count", total_mutations))
			summary_file.write("{:12s}: {:.3f}\n".format("Mean score", mean_score))
			summary_file.write("{:12s}: {:.3f}\n\n".format("Standard dev", stdev_score))

			# Count number of beneficial, deleterious mutations
			del_mutations, ben_mutations, neu_mutations = 0, 0, 0
			for score in step:
				if score > base_score:
					ben_mutations += 1
				elif score < base_score:
					del_mutations += 1
				else:
					neu_mutations += 1
			summary_file.write("{:12}: {:5d} ({:7.2%})\n".format("Deleterious", del_mutations, round(del_mutations / total_mutations, 4)))
			summary_file.write("{:12}: {:5d} ({:7.2%})\n".format("Neutral", neu_mutations, round(neu_mutations / total_mutations, 4)))
			summary_file.write("{:12}: {:5d} ({:7.2%})\n".format("Beneficial", ben_mutations, round(ben_mutations / total_mutations, 4)))
			
		# Write all scores with header for each step
		summary_file.write("\n")
		score_lines = '\n'.join(["Level {}:\n{}".format(i, step) for i, step in enumerate(mutant_scores)])
		summary_file.write(score_lines)

def score_orgs(sensing, run_name, environments, step):
	'''
	Run the given organism in each environment in Avida analyze mode. Score the match between this organism's behavior and the set of environments.
	sensing: True iff organisms evolved with sensing. Used to determine which instruction set file to specify in Avida options. (bool)
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
	if sensing:
		inst_set = "instset-heads-sense.cfg"
	else:
		inst_set = "instset-heads.cfg"
	
	# Initialize list of phenotypic match scores
	mutant_genomes_filename = "../mutant-fitness/{}/gen-step-{}.spop".format(run_name, step)
	with open(mutant_genomes_filename, "r") as genomes:
		skip_legend(genomes)
		n_genomes = 0
		for line in genomes:
			n_genomes += 1
	p_match_list = [0 for i in range(n_genomes)]

	for env in environments:
		# Generate string to describe this environment
		env_str = ""
		for task in env:
			if task < 0:
				env_str += "-"
			elif task > 0:
				env_str += "+"
			else:
				env_str += "0"
		
		# Generate properly configured analyze.cfg file by replacing "%" with arguments
		output_filename = "{}/data/step-{}{}.dat".format(run_name, step, env_str)
		arg_tuple = (mutant_genomes_filename, output_filename) + tuple(env)

		sample_filename = "analyze-mutant-temp.cfg"
		analyze_filename = "data/{}/analyze/step-{}{}.cfg".format(run_name, step, env_str)
		with open(sample_filename, "r") as sample_file, open(analyze_filename, "w") as analyze_file:
			i = 0
			for line in sample_file:
				if "%" in line:
					analyze_file.write(line.replace("%", str(arg_tuple[i])))
					i += 1
				else:
					analyze_file.write(line)
		
		# Run Avida in analyze mode
		subprocess.call("./avida -a -set ANALYZE_FILE {} -def INST_SET {} -set EVENT_FILE events-static.cfg -set VERBOSITY 0".format(analyze_filename, inst_set), shell = True)
		
		# Score phenotypic match in this environment
		with open("data/{}".format(output_filename), "r") as dat_file:
			skip_legend(dat_file)
			for i, line in enumerate(dat_file):
				tasks = [bool(int(t)) for t in line.split()[-1]] # List of bools for whether the organism performed each respective task
				for task, performed in enumerate(tasks):
					if env[task] == 0: # Task is not measured in this environment.
						continue
					if (env[task] > 0 and performed) or (env[task] < 0 and not performed): # Task was correctly regulated for this environment
						p_match_list[i] += 1
					else:
						p_match_list[i] -= 1 # Task was misregulated. This is redundant, but keeps the range of possible scores centered at 0.

	return p_match_list

def recursive_exhaustive_mutants(base_list, instruction_set, steps_remaining):
	'''
	Recursively generate all 1-step mutants of all given genomes
	Include substitution, insertion, and deletion mutations.
	base_list: genomes to mutate from (list of strings)
	instruction_set: set of characters used to represent instructions (string)
	steps_remaining: number of times to mutate given genomes, including this time (int)
	return: list of mutants (list of strings)
	'''
	mutant_list = []
	for genome_str in base_list:
		for i in range(len(genome_str)):
			for inst in instruction_set:
				mutant_list.append(genome_str[:i] + inst + genome_str[i + 1:])	# Substitution
				mutant_list.append(genome_str[:i] + inst + genome_str[i:])		# Insertion
			mutant_list.append(genome_str[:i] + genome_str[i + 1:])				# Deletion

		# Handle insertions at end of genome
		for inst in instruction_set:
			mutant_list.append(genome_str + inst)
	if steps_remaining == 1: # This was the last step
		return mutant_list
	return recursive_exhaustive_mutants(mutant_list, instruction_set, steps_remaining - 1)

def generate_mutants(genome_str, instruction_set, step, n):
	'''
	Generate all or sample of mutants at given Hamming distance from given genome.
	Include substitution, insertion, and deletion mutations.
	Do not exclude redundant and degenerate mutations.
	genome_str: base organism for generating mutants (string with each instruction as one character)
	instruction_set: set of characters used to represent instructions (string)
	step: Hamming distance at which to generate mutants (int)
	n: Number of mutants to explore at this step, 0 for exhaustive (int)
	return: list of mutants (list of strings)
	'''
	if n == 0: # Exhaustive list - will use exorbitant amount of memory if used beyond step 1
		mutant_list = recursive_exhaustive_mutants([genome_str], instruction_set, step)
	else: # Random sample of size n
		mutant_list = []
		for i in range(n):
			new_genome = genome_str
			for s in range(step):
				l = len(new_genome)
				a = len(instruction_set)
				total_n = 2 * l * a + l + a		# Number of possible mutations
				p_sub = l * a / total_n			# Probability of substitution mutation
				p_ins = (l + 1) * a	/ total_n	# Probability of insertion mutation

				# Apply one random mutation
				mut_type = random.random()
				if mut_type < p_sub:
					# Substitution
					new_inst = instruction_set[int(random.random() * a)]
					point = int(random.random() * l)
					new_genome = new_genome[:point] + new_inst + new_genome[point + 1:]
				elif mut_type < p_sub + p_ins:
					# Insertion
					new_inst = instruction_set[int(random.random() * a)]
					point = int(random.random() * (l + 1))
					new_genome = new_genome[:point] + new_inst + new_genome[point:]
				else:
					# Deletion
					point = int(random.random() * l)
					new_genome = new_genome[:point] + new_genome[point + 1:]
			mutant_list.append(new_genome)
	return mutant_list
	
def main(treatments_list, n_runs, tasks_list, exploration):
	'''
	Cycle through all replicates of all treatments to find final dominant org, generate near mutants from these orgs,
	and calculate phenotypic match score for base org and mutants
	treatments_list: names of experimental and control treatments (list of strings)
	tasks_list: tasks rewarded or punished (list of ints)
	exploration: max number of mutants to explore at each Hamming distance from base organism, 0 for exhaustive (dict - ints: ints)
	'''
	# One process to evaluate landscape of each run
	process_list = []
		
	for start, treatment in enumerate(treatments_list):
		# Decide whether to use sensing instruction set
		sensing = (treatment == "Plastic")
		
		for run in range(100 * start, 100 * start + n_runs):
			# Get genome
			filename = "../analysis/{}_{}/final_dom.dat".format(treatment, run)
			with open(filename, "r") as genome_file:
				skip_legend(genome_file)
				genome_str = genome_file.readline().strip()
			
			# Save base genome
			base_filename = "../mutant-fitness/{}_{}/gen-step-0.spop".format(treatment, run)
			with open(base_filename, "w") as genome_file:
				genome_file.write("#filetype genotype_data\n#format sequence\n# 1: Genome Sequence\n\n") # Formatting info for Avida parser
				genome_file.write(genome_str)
			
			# Generate and save mutants at each step (Hamming distance) to explore
			for step in exploration:
				# Generate mutants at this step
				mutant_list = generate_mutants(genome_str, instruction_alphabet, step, exploration[step])
				
				# Output all mutant genomes to file
				mutant_filename = "../mutant-fitness/{}_{}/gen-step-{}.spop".format(treatment, run, step)
				with open(mutant_filename, "w") as mutant_file:
					mutant_file.write("#filetype genotype_data\n#format sequence\n# 1: Genome Sequence\n\n") # Formatting info for Avida parser
					for i, gen in enumerate(mutant_list):
						mutant_file.write("{}\n".format(gen))
			
			step_list = [0] + sorted(list(exploration)) # All steps, including zero, which are explored at all
			# Begin processing
			process_list.append(multiprocessing.Process(target = evaluate_genomes, args = (treatment, run, sensing, tasks_list, step_list)))
			process_list[-1].start()
	
	# Finish processing
	for process in process_list:
		process.join()

# Run with three treatments and 10 runs of each
treatments_list = ["Static", "Changing", "Plastic"]
n_runs = 10
tasks_list = [1, 1, 0, 0, 0, 0, 0, 0, 0] # Tasks that are included
exploration = {1: 100, 2: 100, 3: 100} # Max number of mutants to explore at each Hamming distance from base organism, 0 for exhaustive
main(treatments_list, n_runs, tasks_list, exploration)
