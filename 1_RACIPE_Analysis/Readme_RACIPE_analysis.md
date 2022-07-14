# RACIPE Analysis (Figure 1B, 1C, 1H, 4B, 4C, 4D)
Link to RACIPE package - https://github.com/simonhb1990/RACIPE-1.0

TOPO file for core circuit (Figure 1A): `1_RACIPE_Analysis/TOPO files/Core_circuit_fig1A.topo`
TOPO file for core circuit (Figure S1C): `1_RACIPE_Analysis/TOPO files/Core_circuit_figS1C.topo`

Command to run RACIPE for a particular topo file (No overexpression or down expression):
`RACIPE path_to_topo_file -num_paras 10000 -threads 10`

Command to run RACIPE for a particular topo file (Overexpression):
`RACIPE path_to_topo_file -num_paras 10000 -threads 10 -OEID gene_number -OEFD 10`

- `gene_number` can be obtained from `.cfg` file from basic run of RACIPE without overexpression or down expression. 
- `OEFD` represents an overexpression of 10 fold of the chosen gene.


Command to run RACIPE for a particular topo file (Down expression):
`RACIPE path_to_topo_file -num_paras 10000 -threads 10 -DEID gene_number -DEFD 10`

- `gene_number` can be obtained from `.cfg` file from basic run of RACIPE without overexpression or down expression. 
- `DEFD` represents an down expression of 10 fold of the chosen gene.


## Z Normalization of RACIPE data
The steady state data of RACIPE output is Z normalized within each gene in the data set (as described in methods). The python script to perform Z Normalization is `a2_z_normalisation/c1_performing_z_normalisation.py`

Z Normalization for OE/DE datasets is done with respect to mean and std. dev values for each gene from the core dataset.


## Histogram (Figure 1H, S1D)
Gene expression distribution (Z Normalized) are plotted using the code present in `a3_histograms/`


## Heatmap (Figure 1B, S1E)
Gene expression heatmap (Z Normalized) are plotted using the code present in `a4_heatmaps/`


## PC Analysis and clustering (Figure 1C)
PC analysis and clustering was done using python script in `a5_PCA_plots`

## Scatter plots (Figure S1B)
The Z normalized values were plotted as a scatter plot for pairs of genes using the code present in `a6_scatter_plots`


## TGFB OE/DE Analysis and plotting (Figure S3B)
TGFB over expression/down expression analysis codes are present in `a7_mixed_codes/tgfb_OE_DE`.



