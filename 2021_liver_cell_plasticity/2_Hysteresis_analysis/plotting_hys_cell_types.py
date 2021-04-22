import numpy as np
import matplotlib.pyplot as plt

def data(path_to_file,threshold1,threshold2):
    hep_array = []
    blast_array = []
    chol_array = []
    c = 0
    hep = 0
    blast = 0
    chol = 0
    num = 400
    with open(path_to_file) as f:
        for line in f:
            a = line[:-1].split("\t")
            c += 1
            if c%num == 0:
                hep_array.append(hep)
                blast_array.append(blast)
                chol_array.append(chol)
                hep = 0
                blast = 0
                chol = 0
            if float(a[0])>threshold1:
                hep += 1/num
            elif float(a[1])>threshold1 and float(a[1])<threshold2:
                blast += 1/num
            else:
                chol += 1/num
    return hep_array,blast_array,chol_array

fig, ax = plt.subplots()    
for j in ["Hepatocyte","Hepatoblast"]:
    if j == "not_hep":
        lab = "intially cholangiocytes"
    else:
        lab = "intially hepatocytes"
    y = []
    yerr = []
    for i in range(0,40002,2500):
        path_to_file = "out_"+j+"_"+str(i)+".txt"
        hep_array,blast_array,chol_array = data(path_to_file,6200,17500)
        y.append(np.mean(chol_array))
        yerr.append(np.std(chol_array)) 
    ax.errorbar(list(range(0,40002,2500)), y,
            yerr=yerr,
            fmt='-o',label = lab)
plt.xlabel("Tgf_beta_levels",fontsize=20)
plt.legend(fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.ylabel("Fraction of Cholangiocytes",fontsize=20)
plt.legend([],[], frameon=False)
plt.savefig("cholangiocytes.png",dpi=800)
plt.show()
