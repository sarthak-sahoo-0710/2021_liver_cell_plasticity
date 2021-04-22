

import os
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

def reading_the_cfg_file(path_to_folder):
	for filename in os.listdir(path_to_folder+"r1/"):
		#print(filename[-4:])
		if filename[-4:] == ".cfg":
			name = filename[:-4]
	flag = -1
	d_genes = {}
	with open(path_to_folder+"r1/"+name+".cfg") as f:
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

# function to read all the data for all the state solutions combined and then plot the expression density for each gene, scatter plots for HNF4A and PPARG with or without clustering results
def collating_all_runs_for_z_score_calculation(path_to_dat_files,name,nodes,number_of_states):

    # initializing the state dataframe that stores the expression data for all(1 to 4 state solutions) irrespective of them being in a mono, bi, tri or tetra stable cases
    state_dataframe = pd.DataFrame(columns = nodes)

    # iterating over all the files with mono, bi, tri and tetra stable cases
    for i in range(1,number_of_states+1):
        
        # reading each file separately, getting the sub data structures and appending it to the state dataframe
        data = pd.read_csv(path_to_dat_files+name+"_solution_"+str(i)+".dat",delimiter="\t",header=None)
        
        for j in range(0,i):
            sub_dataframe = data[data.columns[len(nodes)*j+2:len(nodes)*j+(len(nodes)+2)]]
            sub_dataframe.columns = nodes
            state_dataframe = state_dataframe.append(sub_dataframe,ignore_index = True)
    
    mean_of_the_columns = state_dataframe.mean(axis = 0)
    std_of_the_columns = state_dataframe.std(axis = 0)

    return mean_of_the_columns, std_of_the_columns

# function that creates a z_normalised dat file
def creating_z_normalised_dat_files(path_to_dat_files,name,genes,number_of_states,mean_of_the_columns,std_of_the_columns,path_to_output_z_norm):

    # iterating over all the files with mono, bi, tri and tetra stable cases
    for i in range(1,number_of_states+1):
        
        # reading each file separately, getting the sub data structures and appending it to the state dataframe
        data = pd.read_csv(path_to_dat_files+network_name+"_solution_"+str(i)+".dat",delimiter="\t",header=None)
        
        z_normalised_dataframe = pd.DataFrame(columns = sorted(genes))
        z_normalised_dataframe = data[data.columns[0:2]]
        z_normalised_dataframe.columns = ["id","states"]
        
        for j in range(0,i):
            sub_dataframe = data[data.columns[len(genes)*j+2:len(genes)*j+(len(genes)+2)]]
            sub_dataframe.columns = genes
            for gene_name in genes:
                sub_dataframe[gene_name] = (sub_dataframe[gene_name] - mean_of_the_columns[gene_name])/std_of_the_columns[gene_name]
            sub_dataframe = sub_dataframe[sorted(genes)]
            
            z_normalised_dataframe = pd.concat([z_normalised_dataframe, sub_dataframe],axis=1)

        z_normalised_dataframe.to_csv(path_to_output_z_norm+network_name+"_solution_"+str(i)+"_z_norm.dat",sep="\t", header=False, index=False)


def data_array_of_arrays(path_to_folder,name,nodes,d_genes):
	for i in os.listdir(path_to_folder+"r1/"):
		if "solution" in i:
			number = int(i.split("solution_")[1][:-4])
			if number > 4:
				continue
			d = []
			with open(path_to_folder+"r1/"+i) as f:
				for line in f:
					a = line[:-1].split("\t")
					c = 0
					b = []
					for x in a[2:]:
						b.append(float(x))
						c += 1
						if c == nodes:
							d.append(b)
							c = 0
							b =[]
			pca = PCA(n_components=2)
			pca.fit(d)
			X = pca.transform(d)
			print(pca.explained_variance_ratio_)
			plt.scatter(X[:, 0], X[:, 1],s=1)
			plt.show()

path_to_folder = "./../../Modified_core/core_modified/"
name,num_models,nodes,d_genes = reading_the_cfg_file(path_to_folder)
#get_z_stats(path_to_folder,name,nodes,d_genes)
data_array_of_arrays(path_to_folder,name,nodes,d_genes)