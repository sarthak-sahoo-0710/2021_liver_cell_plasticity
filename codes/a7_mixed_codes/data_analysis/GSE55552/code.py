import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

df = pd.read_csv("alt.txt",sep="\t")
print(df)

ax = sns.boxplot(x="Gene", y="expression",hue="cell_type", data=df)
plt.legend([],[], frameon=False)
plt.savefig("GSE55552.png")
plt.show()