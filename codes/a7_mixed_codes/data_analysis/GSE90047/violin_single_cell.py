import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

#df = pd.read_csv("data_single_cell_hepato.txt",sep="\t")
#print(df)
#ax = sns.violinplot(x="Expression", y="Gene", data=df,bw=0.2,split=True)
#plt.savefig("hep_data_pts.png")
#plt.show()

gene = "Tgfbr2"

df = pd.read_csv("hepato_pt.txt",sep="\t")
print(df)

ax = sns.boxplot(x="time", y=gene, data=df)
#plt.legend([],[], frameon=False)
plt.savefig("hepato_pt_"+gene+".png")
plt.show()