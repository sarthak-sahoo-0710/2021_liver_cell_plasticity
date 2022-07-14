## this code counts the number of parameter sets that converge to a given stability regime

import os
import numpy as np

def reading_the_cfg_file(path_to_folder):
	for filename in os.listdir(path_to_folder+"r1/"):
		#print(filename[-4:])
		if filename[-4:] == ".cfg":
			name = filename[:-4]
	with open(path_to_folder+"r1/"+name+".cfg") as f:
		for line in f:
			a = line[:-1].split("\t")
			if a[0] == "NumberOfRACIPEModels":
				num_models = float(a[1])
	return name,num_models

def calculate_plasticity_score(path_to_folder,name):
	name,num_models = reading_the_cfg_file(path_to_folder)
	filename_monostable = name+"_solution_1.dat"
	mean_mono,std_mono = function_run_over_replicates(path_to_folder,filename_monostable)
	plasticity_score_mean = 1 - mean_mono/num_models
	plasticity_score_std = std_mono/num_models
	return plasticity_score_mean,plasticity_score_std

def function_count_number_of_parameter_sets(path_to_filename):
	c = 0
	with open(path_to_filename) as f:
		for line in f:
			c += 1
	return c

def function_run_over_replicates(path_to_folder,filename):
	a = []
	for rep in ["r1","r2","r3"]:
		path_to_filename = path_to_folder+rep+"/"+filename
		count = function_count_number_of_parameter_sets(path_to_filename)
		a.append(count)
	mean = np.mean(a)
	std = np.std(a, ddof = 1)
	return mean,std

def function_run_over_all_stability_regimes(path_to_folder):
	name, num_models = reading_the_cfg_file(path_to_folder)
	quality_of_model_data = 0
	d_freq = {}
	for i in range(1,5):
		filename = name+"_solution_"+str(i)+".dat"
		mean,std = function_run_over_replicates(path_to_folder,filename)
		quality_of_model_data += mean/num_models
		d_freq[i] = [mean,std]

	quality_of_model_data = 1 - quality_of_model_data
	if quality_of_model_data > 0.05:
		print("WARNING! Too many multistable solutions.", quality_of_model_data)

	plasticity_score_mean,plasticity_score_std = calculate_plasticity_score(path_to_folder,name)
	print(path_to_folder,"\t", plasticity_score_mean,"\t",plasticity_score_std)
	
def run_for_all_analysis():
	core_path = "./../../"
	## running for core
	function_run_over_all_stability_regimes(core_path+"core/")
	## running for modified core
	modified_core = ["/","_OE/","_DE/"]
	for i in modified_core:
		function_run_over_all_stability_regimes(core_path+"Modified_core/core_modified"+i)
	## running for self_activation perturbation
	self_act_per = ["core_dys/","c_wo_sa/","t_wo_sa/","s_wo_sa/"]
	for i in self_act_per:
		function_run_over_all_stability_regimes(core_path+"self_activation_perturbation/"+i)
	## running for over/under expression
	over_under = ["core_dys/", "DE/c_DE/", "DE/t_DE/", "DE/s_DE/", "OE/c_OE/", "OE/t_OE/", "OE/s_OE/"]
	for i in over_under:
		function_run_over_all_stability_regimes(core_path+"over_under_expression/"+i)
	

run_for_all_analysis()