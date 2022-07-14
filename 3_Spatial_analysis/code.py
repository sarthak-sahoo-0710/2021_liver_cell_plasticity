import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
threshold1,threshold2 = 6200,17500
a = []
for i in [0,2500,5000,7500,10000,12500,15000,17500,20000,22500,25000,27500,30000,32500,35000,37500,40000]:
    print(i)
    b = []
    with open("./data/out_Hepatocyte_"+str(i)+".txt") as f:
        c = 0
        for line in f:
            c += 1
            if c > 15:
                break
            arr = line[:-1].split("\t")
            if float(arr[0])>threshold1:
                b.append(float(arr[1]))
            elif float(arr[1])>threshold1 and float(arr[1])<threshold2:
                b.append(float(arr[1]))
            else:
                b.append(float(arr[1]))
    a.append(b)

a = np.transpose(a)
cmap = matplotlib.colors.ListedColormap(['blue','orange','green'])
plt.pcolor(a, edgecolors='k', linewidths=4,cmap='Oranges')
plt.savefig("Tgfbr2.png",dpi=800)
plt.show()
