import os
import itertools
import numpy as np
import pandas as pd
import seaborn as sns
from textwrap import wrap
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.signal import argrelextrema
from sklearn.cluster import AgglomerativeClustering

from sklearn.mixture import GaussianMixture

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

def multiple_gaussian_fitting(array_to_be_plotted,gene_name,path_to_plots,num_gaussians,sample_size,max_clusts):

	AIC_all = {}
	BIC_all = {}
	for i in range(1,max_clusts+1):
		AIC_all[i] = []
		BIC_all[i] = []
	for sample in range(1,sample_size):

		random_state = np.random.RandomState(seed=1)

		X = np.array(array_to_be_plotted).reshape(-1, 1)

		# fit models with 1-10 components
		N = np.arange(1, max_clusts+1)
		models = [None for i in range(len(N))]

		for i in range(len(N)):
			models[i] = GaussianMixture(N[i]).fit(X)

		# compute the AIC and the BIC
		AIC = [m.aic(X) for m in models]
		BIC = [m.bic(X) for m in models]

		for idx,val in enumerate(AIC):
			AIC_all[idx+1].append(AIC[idx])
			BIC_all[idx+1].append(BIC[idx])

		M_best = models[num_gaussians-1]               # change number of gaussians here
		s = ""
		a = [M_best.means_[0][0],M_best.means_[1][0],M_best.means_[2][0]]
		arr = sorted(range(len(a)), key=lambda k: a[k])
		for x in arr:
			s += str(M_best.means_[x][0])+"\t"+str(M_best.covariances_[x][0][0])+"\t"+str(M_best.weights_[x])+"\t"
			#print(M_best.means_[0][0],"\t",M_best.means_[1][0],"\t",M_best.means_[2][0],"\t",M_best.covariances_[0][0][0],"\t",M_best.covariances_[1][0][0],"\t",M_best.covariances_[2][0][0],"\t",M_best.weights_[0],"\t",M_best.weights_[1],"\t",M_best.weights_[2])
		print(s)
	fig = plt.figure(figsize=(45, 10))
	fig.subplots_adjust(left=0.12, right=0.97,bottom=0.21, top=0.9, wspace=0.5)


	# plot 1: data + best-fit Mixture
	ax = fig.add_subplot(131)

	x = np.linspace(-6, 6, 1000)
	logprob = M_best.score_samples(x.reshape(-1, 1))
	responsibilities = M_best.predict_proba(x.reshape(-1, 1))
	pdf = np.exp(logprob)
	pdf_individual = responsibilities * pdf[:, np.newaxis]

	ax.hist(X, 30, density=True, histtype='stepfilled', alpha=0.4)
	ax.plot(x, pdf, '-k')
	ax.plot(x, pdf_individual, '--k')
	ax.text(0.04, 0.96, "Best-fit Mixture",ha='left', va='top', transform=ax.transAxes,fontsize=18)
	ax.set_xlim(-3.1,3.1)
	ax.set_xlabel('Z-normalised Expression',fontsize=18)
	ax.set_ylabel('Frequency',fontsize=18)
	plt.xticks(fontsize=18)
	plt.yticks(fontsize=18)

	# removing the first occurance 
	#AIC_all.pop(1)
	#BIC_all.pop(1)

	# plot 2: AIC plots
	ax = fig.add_subplot(132)
	ax.boxplot(AIC_all.values())
	ax.set_xticklabels(AIC_all.keys())
	plt.xlabel('Cluster Count',fontsize=18)
	plt.ylabel('AIC metric',fontsize=18)
	plt.title('AIC Values by Cluster Count',fontsize=18)

	plt.xticks(fontsize=18)
	plt.yticks(fontsize=18)
	

	# plot 3: BIC plots
	ax = fig.add_subplot(133)
	ax.boxplot(BIC_all.values())
	ax.set_xticklabels(BIC_all.keys())
	plt.xlabel('Cluster Count',fontsize=18)
	plt.ylabel('BIC metric',fontsize=18)
	plt.title('BIC Values by Cluster Count',fontsize=18)

	plt.xticks(fontsize=18)
	plt.yticks(fontsize=18)
	plt.savefig(path_to_plots+gene_name+".png")
	plt.close()

def plotting_overall_histograms(state_dataframe,genes,network_name,path_to_plots):
	for g in genes:
		print("Making histogram for: ",g)
		if g in ["Cebpa","Sox9"]:
			gene,threshold_value = find_minima_in_distribution(state_dataframe[g],g,[-3,3],[-0.5,1.5],1000)
			print("high-low minima for "+g+" is: "+str(threshold_value))
		multiple_gaussian_fitting(state_dataframe[g],g,path_to_plots,3,10,7)
		ax = sns.distplot(state_dataframe[g],hist=True,kde_kws={"color": "black"})#, "label": g})
		plt.title(network_name+"_Gene_"+g)
		plt.xlim(-3.5,3)
		plt.xlabel("Z-normalised Expression")
		plt.ylabel("Frequency")
		plt.savefig(path_to_plots+network_name+"_Gene_"+g+"_histogram.png")
		plt.close()

def find_minima_in_distribution(gene_values,gene,overall_range,minima_range,accuracy):
    
    density = gaussian_kde(gene_values)
    
    data = []
    
    for i in range(overall_range[0]*accuracy,(overall_range[1]*accuracy)+1):
        data.append(density(i/accuracy)[0])
        
    lb_index = int(abs((overall_range[0]*accuracy)-(minima_range[0]*accuracy)))
    len_minima_range = int((minima_range[1]-minima_range[0])*accuracy)
    
    minima_index = argrelextrema(np.array(data[lb_index:lb_index+len_minima_range]), np.less)
    
    threshold_value = ((minima_index[0]/accuracy)+minima_range[0])[0]
    #print(overall_range,minima_range,overall_range[0]*accuracy,(overall_range[1]*accuracy)+1,lb_index,lb_index+len_minima_range,data[lb_index],data[lb_index+len_minima_range])
    
    return gene,threshold_value
    
    pass


num_stability_to_consider = 4

def run_for_all_analysis(replicate):
	core_path = "./../../"
	
	## running for core
	"""print("Running for core ...")
				network_name = "core"
				path_to_dat_files = core_path+"core/"+replicate+"/"
				path_to_output_z_norm = path_to_dat_files+"Z_normed/"
				path_to_plots = path_to_dat_files+"plots/histograms/"
				if not os.path.exists(path_to_plots):
					os.makedirs(path_to_plots)
				name,num_models,nodes,d_genes = reading_the_cfg_file(path_to_dat_files)
				genes = list(d_genes.values())
				state_dataframe = collating_all_runs_for_z_score_calculation(path_to_output_z_norm,network_name,genes,num_stability_to_consider)
				plotting_overall_histograms(state_dataframe,genes,network_name,path_to_plots)"""

	## running for modified core
	print("Running for modified core ...")
	"""network_name = "core"
	path_to_dat_files = core_path+"Modified_core/core_modified/"+replicate+"/"
	name,num_models,nodes,d_genes = reading_the_cfg_file(path_to_dat_files)
	genes = list(d_genes.values())
	modified_core = ["/","_OE/","_DE/"]
	for i in modified_core:
		network_name = "core"+i[:-1]
		if "E" in i:
			network_name += "_4"
		path_to_dat_files = core_path+"Modified_core/core_modified"+i+replicate+"/"
		path_to_output_z_norm = path_to_dat_files+"Z_normed/"
		path_to_plots = path_to_dat_files+"plots/histograms/"
		if not os.path.exists(path_to_plots):
			os.makedirs(path_to_plots)
		state_dataframe = collating_all_runs_for_z_score_calculation(path_to_output_z_norm,network_name,genes,num_stability_to_consider)
		plotting_overall_histograms(state_dataframe,genes,network_name,path_to_plots)
	"""
	
	## running for self_activation perturbation
	print("Running for self-activation perturbations: ")
	network_name = "core"
	path_to_dat_files = core_path+"over_under_expression/core/"+replicate+"/"
	name,num_models,nodes,d_genes = reading_the_cfg_file(path_to_dat_files)
	genes = list(d_genes.values())
	path_to_output_z_norm = path_to_dat_files+"Z_normed/"
	path_to_plots = path_to_dat_files+"plots/histograms/"
	if not os.path.exists(path_to_plots):
		os.makedirs(path_to_plots)
	state_dataframe = collating_all_runs_for_z_score_calculation(path_to_output_z_norm,network_name,genes,num_stability_to_consider)
	plotting_overall_histograms(state_dataframe,genes,network_name,path_to_plots)
					
	"""over_under = ["_DE","_OE"]
	for i in over_under:
		for g in genes:
			print(i+"   "+g)
			idx = str(genes.index(g) + 1)
			network_name = "core"+i+"_"+idx
			path_to_dat_files = core_path+"over_under_expression/core"+i+"/"+g+"/"+replicate+"/"
			path_to_output_z_norm = path_to_dat_files+"Z_normed/"
			path_to_plots = path_to_dat_files+"plots/histograms/"
			if not os.path.exists(path_to_plots):
				os.makedirs(path_to_plots)
			state_dataframe = collating_all_runs_for_z_score_calculation(path_to_output_z_norm,network_name,genes,num_stability_to_consider)
			plotting_overall_histograms(state_dataframe,genes,network_name,path_to_plots)"""
							
	## running for over/under expression
	print("Running for over/under expression: ")
	"""network_name = "core"
	path_to_dat_files = core_path+"self_activation_perturbation/core/"+replicate+"/"
	name,num_models,nodes,d_genes = reading_the_cfg_file(path_to_dat_files)
	genes = list(d_genes.values())
					
	per = ["core","core_c_wo_sa","core_t_wo_sa","core_s_wo_sa"]
	for i in per:
		print(i)
		path_to_dat_files = core_path+"self_activation_perturbation/"+i+"/"+replicate+"/"
		name,num_models,nodes,d_genes = reading_the_cfg_file(path_to_dat_files)
		genes = list(d_genes.values())
		path_to_output_z_norm = path_to_dat_files+"Z_normed/"
		path_to_plots = path_to_dat_files+"plots/histograms/"
		if not os.path.exists(path_to_plots):
			os.makedirs(path_to_plots)
		state_dataframe = collating_all_runs_for_z_score_calculation(path_to_output_z_norm,network_name,genes,num_stability_to_consider)
		plotting_overall_histograms(state_dataframe,genes,network_name,path_to_plots)"""
				
run_for_all_analysis("r1")