import os
import itertools
import numpy as np
import pandas as pd
import seaborn as sns
from textwrap import wrap
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.signal import argrelextrema
from sklearn.mixture import GaussianMixture
from sklearn.cluster import AgglomerativeClustering

def reading_the_cfg_file(path_to_folder):
	for filename in os.listdir(path_to_folder):
		#print(filename[-4:])
		if filename[-4:] == ".cfg":
			name = filename[:-4]
	flag = -1
	d_genes = {}
	with open(path_to_folder+name+".cfg") as f:
		for line in f:
			a = line[:-1].split("\t")
			if a[0] == "NumberOfRACIPEModels":
				num_models = float(a[1])
			if a[0] == "NumberOfGenes":
				nodes = int(a[1])
				flag = 0
				continue
			if flag >= 0 and flag < nodes:
				flag += 1
				d_genes[a[0]] = a[1]

	return name,num_models,nodes,d_genes

def collating_all_runs_for_z_score_calculation(path_to_dat_files,network_name,genes,num_stability_to_consider):

    state_dataframe = pd.DataFrame(columns = np.sort(genes))

    # iterating over all the files with mono, bi, tri and tetra stable cases
    for i in range(1,num_stability_to_consider+1):
        
        # reading each file separately, getting the sub data structures and appending it to the state dataframe
        data = pd.read_csv(path_to_dat_files+network_name+"_solution_"+str(i)+"_z_norm.dat",delimiter="\t",header=None)
        
        for j in range(0,i):
            sub_dataframe = data[data.columns[len(genes)*j+2:len(genes)*j+(len(genes)+2)]]
            sub_dataframe.columns = np.sort(genes)
            state_dataframe = state_dataframe.append(sub_dataframe,ignore_index = True)

    return state_dataframe

def plotting_comparision_histograms(dict_state_dataframe,genes,network_name,path_to_plots,perturbed_gene):
	for g in genes:
		print("Making histogram for: ",g)
		for i in dict_state_dataframe:
			ax = sns.distplot(dict_state_dataframe[i][g],hist=False,kde_kws={"label": perturbed_gene+i[:-1]})
		plt.title(network_name+"_Gene_"+g+"_"+perturbed_gene)
		plt.xlim(-3.5,3)
		plt.xlabel("Z-normalised Expression")
		plt.ylabel("Frequency")

		plt.savefig(path_to_plots+"Gene_"+g+"_"+perturbed_gene+"_comp_histogram.png")
		plt.close()


num_stability_to_consider = 4

def run_for_all_analysis(replicate):
	core_path = "./../../"
	
	## running for modified core
	print("Running for modified core ...")
	network_name = "core"
	path_to_dat_files = core_path+"Modified_core/core_modified/"+replicate+"/"
	name,num_models,nodes,d_genes = reading_the_cfg_file(path_to_dat_files)
	genes = list(d_genes.values())
	modified_core = ["/","_OE/","_DE/"]
	dict_state_dataframes = {}
	for i in modified_core:
		network_name = "core"+i[:-1]
		if "E" in i:
			network_name += "_4"
		path_to_dat_files = core_path+"Modified_core/core_modified"+i+replicate+"/"
		path_to_output_z_norm = path_to_dat_files+"Z_normed/"
		
		state_dataframe = collating_all_runs_for_z_score_calculation(path_to_output_z_norm,network_name,genes,num_stability_to_consider)
		dict_state_dataframes[i] = state_dataframe
	path_to_plots = core_path+"Modified_core/core_modified/"+replicate+"/"+"plots/comp_histograms/"
	if not os.path.exists(path_to_plots):
		os.makedirs(path_to_plots)
	plotting_comparision_histograms(dict_state_dataframes,genes,network_name,path_to_plots,"tgfb")
		
	
	## running for self_activation perturbation
	print("Running for self-activation perturbations: ")
	
	over_under = ["_DE","_OE"]

	network_name = "core"
	path_to_dat_files = core_path+"over_under_expression/core/"+replicate+"/"
	name,num_models,nodes,d_genes = reading_the_cfg_file(path_to_dat_files)
	genes = list(d_genes.values())

	for g in genes:
		network_name = "core"
		path_to_dat_files = core_path+"over_under_expression/core/"+replicate+"/"
		name,num_models,nodes,d_genes = reading_the_cfg_file(path_to_dat_files)
		genes = list(d_genes.values())
		dict_state_dataframes = {}
		path_to_output_z_norm = path_to_dat_files+"Z_normed/"
		path_to_plots = path_to_dat_files+"plots/comp_histograms/"
		if not os.path.exists(path_to_plots):
			os.makedirs(path_to_plots)
		state_dataframe = collating_all_runs_for_z_score_calculation(path_to_output_z_norm,network_name,genes,num_stability_to_consider)
		dict_state_dataframes["_ref_"] = state_dataframe
		for i in over_under:
			print(i+"   "+g)
			idx = str(genes.index(g) + 1)
			network_name = "core"+i+"_"+idx
			path_to_dat_files = core_path+"over_under_expression/core"+i+"/"+g+"/"+replicate+"/"
			path_to_output_z_norm = path_to_dat_files+"Z_normed/"
			
			state_dataframe = collating_all_runs_for_z_score_calculation(path_to_output_z_norm,network_name,genes,num_stability_to_consider)
			dict_state_dataframes[i+"_"] = state_dataframe
		plotting_comparision_histograms(dict_state_dataframes,genes,network_name,path_to_plots,g)	
			
	## running for self activation perturbation experiments
	print("Running for self activation perturbation experiments: ")
	network_name = "core"
	path_to_dat_files = core_path+"self_activation_perturbation/core/"+replicate+"/"
	name,num_models,nodes,d_genes = reading_the_cfg_file(path_to_dat_files)
	genes = list(d_genes.values())
	path_to_plots = core_path+"self_activation_perturbation/core/"+replicate+"/"+"plots/comp_histograms/"
	if not os.path.exists(path_to_plots):
		os.makedirs(path_to_plots)
	
	per = ["core_c_wo_sa","core_t_wo_sa","core_s_wo_sa"]
	for i in per:
		dict_state_dataframes = {}
		print("core")
		path_to_dat_files = core_path+"self_activation_perturbation/core/"+replicate+"/"
		name,num_models,nodes,d_genes = reading_the_cfg_file(path_to_dat_files)
		genes = list(d_genes.values())
		path_to_output_z_norm = path_to_dat_files+"Z_normed/"
		
		state_dataframe = collating_all_runs_for_z_score_calculation(path_to_output_z_norm,network_name,genes,num_stability_to_consider)
		dict_state_dataframes["core_"] = state_dataframe

		print(i)
		path_to_dat_files = core_path+"self_activation_perturbation/"+i+"/"+replicate+"/"
		name,num_models,nodes,d_genes = reading_the_cfg_file(path_to_dat_files)
		genes = list(d_genes.values())
		path_to_output_z_norm = path_to_dat_files+"Z_normed/"
		
		state_dataframe = collating_all_runs_for_z_score_calculation(path_to_output_z_norm,network_name,genes,num_stability_to_consider)
		dict_state_dataframes[i+"_"] = state_dataframe
		plotting_comparision_histograms(dict_state_dataframes,genes,network_name,path_to_plots,i)
	
run_for_all_analysis("r1")