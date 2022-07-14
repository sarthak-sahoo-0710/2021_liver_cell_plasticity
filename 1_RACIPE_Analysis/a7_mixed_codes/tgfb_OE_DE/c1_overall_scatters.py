import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# gene expression levels of various clusters


def grouped_barplot(df, cat,subcat, val , err):
    u = df[cat].unique()
    x = np.arange(len(u))
    subx = df[subcat].unique()
    offsets = (np.arange(len(subx))-np.arange(len(subx)).mean())/(len(subx)+1.)
    width= np.diff(offsets).mean()
    for i,gr in enumerate(subx):
        dfg = df[df[subcat] == gr]
        plt.bar(x+offsets[i], dfg[val].values, width=width, 
                label="{} {}".format("",gr), yerr=dfg[err].values)
    plt.xlabel(cat, fontsize=20)
    plt.ylabel(val, fontsize=20)
    plt.xticks(x, u, fontsize=20)
    plt.yticks(fontsize=20)
    #plt.ylim([0,0.5])
    #plt.legend(fontsize=20)
    #plt.tight_layout()
    #plt.figure(figsize=(20,10))
    plt.savefig("sox9_hep.png")
    plt.show()

df = pd.read_csv("sox9_hep.txt")
print(df)

cat = "Candidate"
subcat = "Sample_Set"
val = "Fraction_of_cases"
err = "Error"
plt.bar(df[subcat],df[val],yerr = df[err],color = ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854', '#ffd92f', '#e5c494', '#b3b3b3'])
plt.ylim([0,0.3])
plt.savefig("sox9_hep.png",dpi = 800)
plt.show()
#grouped_barplot(df, cat, subcat, val, err )
