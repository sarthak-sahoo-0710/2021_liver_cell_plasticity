import os
import itertools
import numpy as np
import collections
import pandas as pd
import seaborn as sns
from gap_statistic import OptimalK
from textwrap import wrap
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.signal import argrelextrema
from sklearn.cluster import AgglomerativeClustering
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn import metrics
from scipy.spatial.distance import cdist
from sklearn.metrics import silhouette_samples, silhouette_score
from sklearn.mixture import GaussianMixture
from scipy.stats import ttest_ind

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
            sub_dataframe.columns = genes
            state_dataframe = state_dataframe.append(sub_dataframe,ignore_index = True)

    return state_dataframe

def kmeans_silhouette(principalDf,num_clusters,path_to_plots):
    score_array = [num_clusters]
    for N in range(15):
        range_n_clusters = num_clusters
        x = []
        for n_clusters in range_n_clusters:
            clusterer = KMeans(n_clusters=n_clusters)
            preds = clusterer.fit_predict(principalDf)
            centers = clusterer.cluster_centers_
            score = silhouette_score(principalDf, preds)
            x.append(score)
        score_array = np.vstack((score_array,x))
    #print(score_array)
    silhouette_score_data_frame = pd.DataFrame(score_array)
    silhouette_score_data_frame = silhouette_score_data_frame.drop([0])
    
    column_name=[]
    for i, j in zip(num_clusters, range(len(num_clusters))):
        name = str(i) + ' clusters'
        silhouette_score_data_frame = silhouette_score_data_frame.rename(columns={j:name})
        column_name.append(name)
    
    #plotting the boxplot to find the ideal number of clusters for k means
    fig, ax = plt.subplots()
    plot = silhouette_score_data_frame.boxplot(column=column_name,ax=ax)
    ax.set(xlabel='Clusters', ylabel='mean silhouette coeff', title='Boxplot grouped by cluster pca_clustering\n'+network_name)
    plt.savefig(path_to_plots+network_name+'_' +"boxplot for silhouette analysis_pca_clustsering.png",dpi=1000)
    plt.close()

def information_criterion(principalDf,path_to_plots):
    n_components = np.arange(1, 10)
    models = [GaussianMixture(n, covariance_type='full', random_state=0).fit(principalDf)
              for n in n_components]

    plt.plot(n_components, [m.bic(principalDf) for m in models], label='BIC')
    plt.plot(n_components, [m.aic(principalDf) for m in models], label='AIC')
    plt.legend(loc='best')
    plt.xlabel('n_components')
    plt.savefig(path_to_plots+'_' +"info_crit_pca_clustsering.png",dpi=1000)
    plt.close()

def gene_levels(state_dataframe_PCA,genes):
    cluster = AgglomerativeClustering(n_clusters=3, affinity='euclidean', linkage='ward')
    cluster.fit_predict(state_dataframe_PCA)
    state_dataframe_new = state_dataframe_PCA.copy(deep=True)
    state_dataframe_new["cluster_label"] = cluster.labels_
    cluster_0 = state_dataframe_new[state_dataframe_new["cluster_label"] == 0]
    cluster_1 = state_dataframe_new[state_dataframe_new["cluster_label"] == 1]
    cluster_2 = state_dataframe_new[state_dataframe_new["cluster_label"] == 2]

    plt.scatter(cluster_0["Cebpa"],cluster_0["Tgfbr2"])
    plt.show()
    plt.close()

    for g in genes:
        print("cluster_0", g, np.mean(cluster_0[g]), np.std(cluster_0[g]))
        print("cluster_1", g, np.mean(cluster_1[g]), np.std(cluster_1[g]))
        print("cluster_2", g, np.mean(cluster_2[g]), np.std(cluster_2[g]))

        ttest_0_1 = ttest_ind(cluster_0[g],cluster_1[g],equal_var = True)
        ttest_1_2 = ttest_ind(cluster_1[g],cluster_2[g],equal_var = True)
        ttest_0_2 = ttest_ind(cluster_0[g],cluster_2[g],equal_var = True)
        print(ttest_0_1[0],ttest_0_1[1],ttest_1_2[0],ttest_1_2[1],ttest_0_2[0],ttest_0_2[1])

def PCA_analysis(state_dataframe_PCA,genes,network_name,path_to_plots):
    num_clusters = [2,3,4,5,6,7,8]
    #scaling the normalized z dataframe to reduce mean to 0 and std to 1
    scaled_state_dataframe = StandardScaler().fit_transform(state_dataframe_PCA)
    
    #performing pca
    pca = PCA(n_components = 2)
    principalComponents = pca.fit_transform(scaled_state_dataframe)
    principalDf = pd.DataFrame(data = principalComponents
             , columns = ['principal component 1', 'principal component 2'])

    cluster = AgglomerativeClustering(n_clusters=3, affinity='euclidean', linkage='ward')
    cluster.fit_predict(state_dataframe_PCA)

    gene_levels(state_dataframe_PCA,genes)

    elements_count = collections.Counter(cluster.labels_)
    # printing the element and the frequency
    for key, value in elements_count.items():
       print(f"{key}   {value}")
    
    PrincipalComponents = principalDf
    
    #to find the variance of each prinicpal component and the contribuition of each gene in the PCA
    explained_variance=pca.explained_variance_
    explained_variance_ratio=pca.explained_variance_ratio_
    
    pca_components = pd.DataFrame(pca.components_,columns=state_dataframe_PCA.columns,index=['PC-1','PC-2'])
    
    fig, ax = plt.subplots()
    fig = principalDf.plot(kind='scatter',x='principal component 1', y='principal component 2', c='black',s=0.5,ax=ax)
    ax.set_xlabel('PC-1('+str(round(explained_variance_ratio[0],4)*100)+'% variance)')
    ax.set_ylabel('PC-2('+str(round(explained_variance_ratio[1],4)*100)+'% variance)')
    #plt.show()
    plt.savefig(path_to_plots+network_name+"PCA_scatter_black.png")
    plt.close()

    fig, ax = plt.subplots()
    fig = principalDf.plot(kind='scatter',x='principal component 1', y='principal component 2', c=cluster.labels_, cmap='rainbow',s=0.5,ax=ax)
    ax.set_xlabel('PC-1('+str(round(explained_variance_ratio[0],4)*100)+'% variance)')
    ax.set_ylabel('PC-2('+str(round(explained_variance_ratio[1],4)*100)+'% variance)')
    #plt.show()
    plt.savefig(path_to_plots+network_name+"PCA_scatter.png")
    plt.close()

    print(explained_variance_ratio)
    print(pca_components)

    #scatter plots of PCA analysis for each gene expression
    """for gene in genes:
        fig, ax = plt.subplots()
        fig = principalDf.plot(kind='scatter',x='principal component 1', y='principal component 2',c=state_dataframe_PCA[gene],cmap="cool",s=0.5,ax=ax)
        ax.set_xlabel('PC-1('+str(round(explained_variance_ratio[0],4)*100)+'% variance)')
        ax.set_ylabel('PC-2('+str(round(explained_variance_ratio[1],4)*100)+'% variance)')
        plt.title(gene)
        #plt.show()
        plt.savefig(path_to_plots+network_name+gene+"_PCA_scatter.png",dpi=1000)
        plt.close()"""
        
    return principalDf,pca_components

def gap_stat(df_gap_stat,path_to_plots):
    d = {}
    for i in range(10):
        d[i+1] = []
    for i in range(1,100):
        optimalK = OptimalK(parallel_backend='none')  #you can change this to rust
        n_clusters = optimalK(df_gap_stat, cluster_array=np.arange(1, 11))
        for idx,j in enumerate(optimalK.gap_df.gap_value):
            d[idx+1].append(j)
    fig, ax = plt.subplots()
    ax.boxplot(d.values())
    ax.set_xticklabels(d.keys())
    plt.xlabel('Cluster Count')
    plt.ylabel('Gap Value')
    plt.title('Gap Values by Cluster Count')
    #plt.show()
    plt.savefig(path_to_plots+"gap_stat.png",dpi=1000)
    plt.close()
    

num_stability_to_consider = 4

def run_for_all_analysis(replicate):
    core_path = "./../../"
    
    ## running for core
    """print("Running for core ...")
                network_name = "core"
                path_to_dat_files = core_path+"core/"+replicate+"/"
                path_to_output_z_norm = path_to_dat_files+"Z_normed/"
                path_to_plots = path_to_dat_files+"plots/PCA/"
                if not os.path.exists(path_to_plots):
                    os.makedirs(path_to_plots)
                name,num_models,nodes,d_genes = reading_the_cfg_file(path_to_dat_files)
                genes = list(d_genes.values())
                state_dataframe = collating_all_runs_for_z_score_calculation(path_to_output_z_norm,network_name,genes,num_stability_to_consider)
                PCA_analysis(state_dataframe,genes,network_name,path_to_plots)"""

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
        path_to_plots = path_to_dat_files+"plots/PCA/"
        if not os.path.exists(path_to_plots):
            os.makedirs(path_to_plots)
        state_dataframe = collating_all_runs_for_z_score_calculation(path_to_output_z_norm,network_name,genes,num_stability_to_consider)
        PCA_analysis(state_dataframe,genes,network_name,path_to_plots)"""
        
    
    ## running for self_activation perturbation
    print("Running for self-activation perturbations: ")
    network_name = "core"
    path_to_dat_files = core_path+"over_under_expression/core/"+replicate+"/"
    name,num_models,nodes,d_genes = reading_the_cfg_file(path_to_dat_files)
    genes = list(d_genes.values())
    path_to_output_z_norm = path_to_dat_files+"Z_normed/"
    path_to_plots = path_to_dat_files+"plots/PCA/"
    if not os.path.exists(path_to_plots):
        os.makedirs(path_to_plots)
    state_dataframe = collating_all_runs_for_z_score_calculation(path_to_output_z_norm,network_name,genes,num_stability_to_consider)
    PCA_analysis(state_dataframe,genes,network_name,path_to_plots)
                    
    """over_under = ["_DE","_OE"]
    for i in over_under:
        for g in genes:
            print(i+"   "+g)
            idx = str(genes.index(g) + 1)
            network_name = "core"+i+"_"+idx
            path_to_dat_files = core_path+"over_under_expression/core"+i+"/"+g+"/"+replicate+"/"
            path_to_output_z_norm = path_to_dat_files+"Z_normed/"
            path_to_plots = path_to_dat_files+"plots/PCA/"
            if not os.path.exists(path_to_plots):
                os.makedirs(path_to_plots)
            state_dataframe = collating_all_runs_for_z_score_calculation(path_to_output_z_norm,network_name,genes,num_stability_to_consider)
            PCA_analysis(state_dataframe,genes,network_name,path_to_plots)"""
                            
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
        path_to_plots = path_to_dat_files+"plots/PCA/"
        if not os.path.exists(path_to_plots):
            os.makedirs(path_to_plots)
        state_dataframe = collating_all_runs_for_z_score_calculation(path_to_output_z_norm,network_name,genes,num_stability_to_consider)
        PCA_analysis(state_dataframe,genes,network_name,path_to_plots)"""
                
run_for_all_analysis("r3")