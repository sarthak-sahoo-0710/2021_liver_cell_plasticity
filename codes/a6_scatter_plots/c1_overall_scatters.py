import os
import itertools
import numpy as np
import pandas as pd
import seaborn as sns
from textwrap import wrap
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
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

    state_dataframe = pd.DataFrame(columns = genes)

    # iterating over all the files with mono, bi, tri and tetra stable cases
    for i in range(1,num_stability_to_consider+1):
        
        # reading each file separately, getting the sub data structures and appending it to the state dataframe
        data = pd.read_csv(path_to_dat_files+network_name+"_solution_"+str(i)+"_z_norm.dat",delimiter="\t",header=None)
        
        for j in range(0,i):
            sub_dataframe = data[data.columns[len(genes)*j+2:len(genes)*j+(len(genes)+2)]]
            sub_dataframe.columns = genes
            state_dataframe = state_dataframe.append(sub_dataframe,ignore_index = True)

    return state_dataframe

def plotting_overall_scatter_plots(state_dataframe,genes,path_to_plots):
	for i in range(len(genes)):
		for j in range(i+1,len(genes)):
			gene1 = genes[i]
			gene2 = genes[j]
			sc = sns.regplot(state_dataframe[gene1],state_dataframe[gene2],marker="o",color='black', scatter_kws={'s':1})
			#plt.title()
			plt.xlabel(gene1)
			plt.ylabel(gene2)
			plt.ylim(-3, 2)
			plt.xlim(-3, 2)
			plt.title(spearmanr(state_dataframe[gene1],state_dataframe[gene2]))
			plt.savefig(path_to_plots+gene1+"_"+gene2+"_scatter.png",dpi=1000)
			plt.close()

def plotting_overall_scatter_plots_with_clustering(state_dataframe,genes,path_to_plots,num_clusters):
	#cluster = AgglomerativeClustering(n_clusters=num_clusters, affinity='euclidean', linkage='ward')
	#cluster.fit_predict(state_dataframe)
	for i in range(len(genes)):
		for j in range(i+1,len(genes)):
			gene1 = genes[i]
			gene2 = genes[j]
			sub_dataframe_scatter = state_dataframe[[gene1, gene2]]
			cluster = AgglomerativeClustering(n_clusters=num_clusters, affinity='euclidean', linkage='ward')
			cluster.fit_predict(sub_dataframe_scatter)
			plt.scatter(state_dataframe[gene1],state_dataframe[gene2], c=cluster.labels_, cmap='rainbow',marker="o",s=0.5)
			
			plt.title("Clusters: "+str(num_clusters))
			plt.xlabel(gene1)
			plt.ylabel(gene2)
			plt.ylim(-3, 2)
			plt.xlim(-3, 2)
			plt.savefig(path_to_plots+gene1+"_"+gene2+"_scatter_"+str(num_clusters)+"clusters.png",dpi=1000)
			plt.close()

num_stability_to_consider = 4
replicate = "r1"
core_path = "./../../"
network_name = "core"
path_to_dat_files = core_path+"over_under_expression/core/"+replicate+"/"
path_to_output_z_norm = path_to_dat_files+"Z_normed/"
name,num_models,nodes,d_genes = reading_the_cfg_file(path_to_dat_files)
genes = list(d_genes.values())
state_dataframe = collating_all_runs_for_z_score_calculation(path_to_output_z_norm,network_name,genes,num_stability_to_consider)
path_to_plots = path_to_dat_files+"plots/scatters/"
if not os.path.exists(path_to_plots):
	os.makedirs(path_to_plots)
plotting_overall_scatter_plots(state_dataframe,genes,path_to_plots)
#plotting_overall_scatter_plots_with_clustering(state_dataframe,genes,path_to_plots,3)