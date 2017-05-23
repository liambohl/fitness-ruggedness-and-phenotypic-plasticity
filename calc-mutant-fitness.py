###############################################################################
# Iterate over all treatments and all replicates to find specified generations of offspring
# of the final dominant genotype. Evaluate the phenotypic match score of the
# final dominant genotypes and offspring across all environments they would
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
	f: file pointer object
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
	environments: all permutations of values for previous tasks - [[int]]
	tasks_remaining: which tasks are measured in this experiment - [int]
	returns: expanded list of environments - [[int]]
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

	# Call next generation and return
	if len(tasks_remaining) > 1:
		return permute_environments(tasks_remaining[1:], new_environments)
	else:
		return new_environments

def score_orgs(sensing, run_name, environments, generation, n):
	'''
	Using Avida analyze mode, record tasks performed by each member of the specified generation of offspring in each specified environment.
		Score how well each organism matches its behavior to the environment (phenotypic match score)
	sensing: True iff organisms evolved with sensing. Determines which instruction set file to specify in Avida options. - bool
	run_name: treatment and replicate number - string
	environments: list of environments; each environment is a list of reaction values for each task - [[int]]
	generation: number of generations into the future from the base genome - int
	n: number of sample offspring that have been measured for this generation - int
	returns: list of phenotypic match scores for sample offspring in this generation - [int]
	'''
	# Set name of instruction set file
	if sensing:
		inst_set = "instset-heads-sense.cfg"
	else:
		inst_set = "instset-heads.cfg"
	
	# Initialize list of phenotypic match scores
	p_match_scores = [0 for i in range(n)]

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
		offspring_genomes_filename = "../offspring-fitness/{}/gen-step-{}.spop".format(run_name, generation)
		output_filename = "{}/data/step-{}{}.dat".format(run_name, generation, env_str)
		arg_tuple = (offspring_genomes_filename, output_filename) + tuple(env)

		sample_filename = "analyze-mutant-temp.cfg"
		analyze_filename = "data/{}/analyze/step-{}{}.cfg".format(run_name, generation, env_str)
		with open(sample_filename, "r") as sample_file, open(analyze_filename, "w") as analyze_file:
			i = 0
			for line in sample_file:
				if "%" in line:
					analyze_file.write(line.replace("%", str(arg_tuple[i])))
					i += 1
				else:
					analyze_file.write(line)
		
		# Run Avida in analyze mode
		subprocess.call("./avida -a -set ANALYZE_FILE {} -def INST_SET {} -set EVENT_FILE events/template.cfg -set VERBOSITY 0".format(analyze_filename, inst_set), shell = True)
		
		# Score phenotypic match in this environment
		with open("data/{}".format(output_filename), "r") as dat_file:
			skip_legend(dat_file)
			for i, line in enumerate(dat_file):
				tasks_str = line.split()[-1]
				tasks = [bool(int(t)) for t in tasks_str] # List of bools for whether the organism performed each respective task
				for task, performed in enumerate(tasks):
					if env[task] == 0: # Task is not measured in this environment.
						continue
					if (env[task] > 0 and performed) or (env[task] < 0 and not performed): # Task was correctly regulated for this environment
						p_match_scores[i] += 1
					else: # Task was misregulated. This is redundant, but keeps the range of possible scores centered at 0.
						p_match_scores[i] -= 1

	return p_match_scores

def summarize_landscape(treatment, run, max_score, base_score, offspring_scores):
	'''
	Summarize this mutational landscape in summary.txt.
		More detail needed!
	treatment: name of this treatment - string
	run: index of this replicate - int
	max_score: Highest possible phenotypic match score in this treatment - int
	base_score: Phenotypic match score of base organism - int
	offspring_scores: - Phenotypic match score of each sample offspring from each generation - [[int]]
	'''
	# Output to file
	out_filename = "../offspring-fitness/{}_{}/summary.txt".format(treatment, run)
	with open(out_filename, "w") as summary_file:
		# Max score and base score
		summary_file.write("{:12s}: {:6.3f}\n".format("Max score", float(max_score)))
		summary_file.write("{:12s}: {:6.3f}\n".format("Base org", float(base_score)))
		
		for i, generation in enumerate(offspring_scores):
			# Header
			summary_file.write("\n---- Step {} -----\n".format(i + 1))
			
			# Summary Stats
			total_mutations = len(generation)
			mean_score = numpy.mean(generation)
			stdev_score = numpy.std(generation)
			summary_file.write("{:12s}: {:6d}\n".format("Sample Size", total_mutations))
			summary_file.write("{:12s}: {:6.3f}\n".format("Mean score", mean_score))
			summary_file.write("{:12s}: {:6.3f}\n\n".format("Standard dev", stdev_score))

			# Count number of beneficial, deleterious mutations
			del_mutations, ben_mutations, neu_mutations = 0, 0, 0
			for score in generation:
				if score > base_score:
					ben_mutations += 1
				elif score < base_score:
					del_mutations += 1
				else:
					neu_mutations += 1
			p_del = round(del_mutations / total_mutations, 4)
			p_neu = round(neu_mutations / total_mutations, 4)
			p_ben = round(ben_mutations / total_mutations, 4)

			summary_file.write("{:12}: {:6.2%}\n".format("Deleterious", p_del))
			summary_file.write("{:12}: {:6.2%}\n".format("Neutral", p_neu))
			summary_file.write("{:12}: {:6.2%}\n".format("Beneficial", p_ben))
			
		# Write all scores with header for each generation
		summary_file.write("\n")
		score_lines = '\n'.join(["generation {}:\n{}".format(i + 1, generation) for i, generation in enumerate(offspring_scores)])
		summary_file.write(score_lines)

def dict_add(d, x):
	'''
	Add x to dict d or increase d[x]
	d: dict of unspecified type - {?: int}
	x: item to be counted - same type as dictionary keys
	'''
	if x in d:
		d[x] += 1
	else:
		d[x] = 1

def generate_offspring(genome_str, instruction_set, generation, n, mutation_rates):
	'''
	Generate a random sample of possible descendents at the given generation from the given genome
	genome_str: base organism genome with one character per instruction - string
	instruction_set: set of characters used to represent instructions - string
	generation: number of generations into the future from the base genome - int
	n: number of offspring to generate from this generation - int
	mutation_rates: probabilities of specific mutation types - see end of script for detailed description - {string: float}
	return: dict of possible offspring genomes and their frequencies - {string: int}
	'''
	a = len(instruction_set)

	offspring_dict = {} # Offspring and frequencies for this generation
	for i in range(n):
		
		new_genome = genome_str
		# Apply mutation operations once for each generation
		for s in range(generation):

			l = len(new_genome)

			# Substitution
			for point in range(len(new_genome)):
				if random.random() < mutation_rates["copy_mut"]:
					new_inst = instruction_set[int(random.random() * a)]
					new_genome = new_genome[:point] + new_inst + new_genome[point + 1:]
			
			# Insertion
			if random.random() < mutation_rates["divide_ins"]:
				point = int(random.random() * (l + 1))
				new_inst = instruction_set[int(random.random() * a)]
				new_genome = new_genome[:point] + new_inst + new_genome[point:]
			
			# Deletion
			if random.random() < mutation_rates["divide_del"]:
				point = int(random.random() * l)
				new_genome = new_genome[:point] + new_genome[point + 1:]

		dict_add(offspring_dict, new_genome) # Add the new genome
	return offspring_dict

def save_genomes(treatment, run, exploration, mutation_rates):
	'''
	For base organism and each future generation, generate a .spop file with nothing but genome sequences.
		These are used for evaluating each organism's phenotype.
	treatment: name of this treatment - string
	run: index of this replicate - int
	exploration: Number of offspring that have been measured at each generation after the base organism - {int: int}
	mutation_rates: probabilities of specific mutation types - see end of script for detailed description - {string: float}
	'''
	# Get genome
	filename = "../analysis/{}_{}/final_dom.dat".format(treatment, run)
	with open(filename, "r") as genome_file:
		skip_legend(genome_file)
		genome_str = genome_file.readline().strip()
	
	# Formatting info for Avida parser
	formatting_lines = "#filetype genotype_data\n#format sequence\n\n# Legend:\n# 1: Genome Sequence\n\n"
	
	# Save base genome
	base_filename = "../offspring-fitness/{}_{}/gen-step-0.spop".format(treatment, run)
	with open(base_filename, "w") as genome_file:
		genome_file.write(formatting_lines)
		genome_file.write(genome_str)
	
	# Generate and save all specified offspring
	for generation in exploration:
		offspring_dict = generate_offspring(genome_str, instruction_alphabet, generation, exploration[generation], mutation_rates)
		# Output offspring genomes to file
		offspring_filename = "../offspring-fitness/{}_{}/gen-step-{}.spop".format(treatment, run, generation)
		with open(offspring_filename, "w") as offspring_file:
			offspring_file.write(formatting_lines)
			for gen in offspring_dict:
				# Write the offspring genome as many times as it was generated.
				for i in range(offspring_dict[gen]):
					offspring_file.write("{}\n".format(gen))
			
def explore_landscape(treatment, run, tasks_list, exploration, mutation_rates):
	'''
	Find final dominant organism for this run, explore its potential offspring through specified generations,
		collect and summarize data on their phenotypes.
	treatment: name of this treatment - string
	run: index of this replicate - int
	tasks_list: which tasks are measured in this treatment - [int]
	exploration: Number of offspring that have been measured at each generation after the base organism - {int: int}
	mutation_rates: probabilities of specific mutation types - see end of script for detailed description - {string: float}
	'''
	# Generate .spop files with base org and sample offspring genomes
	save_genomes(treatment, run, exploration, mutation_rates)
	
	# Build list of environments, find max possible score
	n_tasks = len(tasks_list) - tasks_list.count(0)
	if treatment == "Static":
		environments = [tasks_list]
		max_score = n_tasks
	else:
		environments = permute_environments(tasks_list)
		max_score = n_tasks * 2 ** n_tasks # Number of measured tasks x number of environments
	
	# Extra scoring variables
	sensing = (treatment == "Plastic") # If sensing is true, organisms are given sensing instructions.
	run_name = "{}_{}".format(treatment, run) # Name of this run

	# Get phenotypic match scores for base organism and all offspring.
	base_score = score_orgs(sensing, run_name, environments, 0, 1)[0]
	offspring_scores = [score_orgs(sensing, run_name, environments, generation, exploration[generation]) for generation in exploration]
	
	# Summarize base org and offspring phenotypes for this run
	summarize_landscape(treatment, run, max_score, base_score, offspring_scores)
	
def main(treatments_list, n_runs, tasks_list, exploration, mutation_rates):
	'''
	Cycle through all replicates of all treatments to find final dominant org, generate near offspring from these orgs,
		and calculate phenotypic match score for base org and offspring
	treatments_list: names of experimental and control treatments - [string]
	tasks_list: tasks rewarded or punished - [int]
	exploration: max number of offspring to explore at each generation - {int: int}
	'''
	# One process to evaluate landscape of each run
	process_list = []
		
	for start, treatment in enumerate(treatments_list):
		for run in range(100 * start, 100 * start + n_runs):
			# Begin processing
			process_list.append(multiprocessing.Process(target = explore_landscape, args = (treatment, run, tasks_list, exploration, mutation_rates)))
			process_list[-1].start()
	
	# Finish processing
	for process in process_list:
		process.join()

# Mutation rates used in Avida trials and in generating sample offspring
# copy_mut: Probability, per site, of substituting the instruction for a random instruction
# divide_ins: Probability of inserting a random instruction at a random location in the genome
# divide_del: Probability of deleting the instruction at a random location in the genome
mutation_rates = {"copy_mut": 0.0075, "divide_ins": 0.05, "divide_del": 0.05}

# Run with three treatments and 10 runs of each
treatments_list = ["Static", "Changing", "Plastic"]
n_runs = 10
tasks_list = [1, 1, 0, 0, 0, 0, 0, 0, 0] # Tasks that are included
exploration = {1: 10000, 2: 10000, 3: 10000} # Max number of offspring to explore at each generation
main(treatments_list, n_runs, tasks_list, exploration, mutation_rates)
