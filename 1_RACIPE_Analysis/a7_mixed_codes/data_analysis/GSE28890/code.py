import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

gene = "TGFB1"

df = pd.read_csv("negative.txt",sep="\t")
print(df)

ax = sns.boxplot(x="time", y=gene, data=df)
#plt.legend([],[], frameon=False)
plt.savefig("GSE28890_negative_"+gene+".png")
plt.show()

df_1 = pd.read_csv("positive.txt",sep="\t")
print(df_1)

ax = sns.boxplot(x="time", y=gene, data=df_1)
#plt.legend([],[], frameon=False)
plt.savefig("GSE28890_positive_"+gene+".png")
plt.show()

df_2 = pd.read_csv("all.txt",sep="\t")
print(df_2)

ax = sns.boxplot(x="case",hue="time", y=gene, data=df_2)
plt.legend([],[], frameon=False)
plt.savefig("GSE28890_all_"+gene+".png")
plt.show()