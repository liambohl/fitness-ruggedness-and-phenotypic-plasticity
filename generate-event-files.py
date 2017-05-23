###############################################################################
# Generate events.cfg files for static environment and for each replicate of
# the changing environment treatments.
###############################################################################

import random

def generate_changing_file(task_dict, template_file, out_filename, env_freq, run_time):
	'''
	Make a complete event file for one run of a changing environment treatment.
	tasks_dict: absolute value of reaction value for each task - {str: int}
	template_file: .cfg file which contains basic info - file pointer
	out_filename: name of .cfg file to create - str
	env_freq: frequency of environment switching, measured in updates - int
	run_time: total length of Avida run in updates - int
	'''
	with open(out_filename, "w") as out_file:

		# Copy basic info
		for line in template_file:
			out_file.write(line)
		
		# Set reaction value for each task
		for update in range(0, run_time, env_freq):
			out_file.write("# Update {}:\n".format(update))
			for task in task_dict:
				value = task_dict[task] * random.choice([1, -1])
				line = "u {:6d} setReactionValue {:4s} {:3.1f}\n".format(0, task, value)
				out_file.write(line)
			out_file.write("\n")

def generate_event_files(task_dict, n_runs, env_freq = 50, run_time = 100000):
	'''
	Make complete event file for static environment. Use random seed of each
		changing run to make event file with random task value switching.
	tasks_dict: absolute value of reaction value for each task - {str: int}
	n_runs: number of replicates for each treatment - int
	env_freq: frequency of environment switching, measured in updates - int
	run_time: total length of Avida run in updates - int
	'''
	with open("template.cfg", "r") as template_file:

		# Static environment
		with open("events/static.cfg", "w") as out_file:

			# Copy basic info
			for line in template_file:
				out_file.write(line)

			# Set reaction value for each task
			for task in task_dict:
				value = task_dict[task]
				line = "u {:6d} setReactionValue {:4s} {:3.1f}\n".format(0, task, value)
				out_file.write(line)

		# Changing environments
		for run in range(n_runs):
			
			# Changing treatment
			seed = 100 + run
			random.seed(seed)
			filename = "events/changing_{}.cfg".format(seed)
			generate_changing_file(task_dict, template_file, filename, env_freq, run_time)
			
			
			# Plastic treatment
			seed = 200 + run
			random.seed(seed)
			filename = "events/plastic_{}.cfg".format(seed)
			generate_changing_file(task_dict, template_file, filename, env_freq, run_time)
		

generate_event_files({"NOT": 1, "NAND": 1, "AND": 1, "ORN": 1, "OR": 0, "ANDN": 0, "NOR": 0, "XOR": 0, "EQU": 0}, 10)
