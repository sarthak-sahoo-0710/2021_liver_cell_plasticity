import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

bf = []
sox9 = []
tgfbr2 = []
cebpa = []
with open("branch_1.txt") as f:
    for line in f:
        a = line[:-1].split("\t")
        bf.append(float(a[0]))
        sox9.append(float(a[3]))
        cebpa.append(float(a[1]))
        tgfbr2.append(float(a[2]))

plt.plot(bf,cebpa,"b-",linewidth=2.5)

bf = []
sox9 = []
tgfbr2 = []
cebpa = []
with open("branch_2.txt") as f:
    for line in f:
        a = line[:-1].split("\t")
        bf.append(float(a[0]))
        sox9.append(float(a[3]))
        cebpa.append(float(a[1]))
        tgfbr2.append(float(a[2]))

plt.plot(bf,cebpa,"r--",linewidth=2.5)

bf = []
sox9 = []
tgfbr2 = []
cebpa = []
with open("branch_3.txt") as f:
    for line in f:
        a = line[:-1].split("\t")
        bf.append(float(a[0]))
        sox9.append(float(a[3]))
        cebpa.append(float(a[1]))
        tgfbr2.append(float(a[2]))

plt.plot(bf,cebpa,"b-",linewidth=2.5)

bf = []
sox9 = []
tgfbr2 = []
cebpa = []
with open("branch_1_hep.txt") as f:
    for line in f:
        a = line[:-1].split("\t")
        bf.append(float(a[0]))
        sox9.append(float(a[3]))
        cebpa.append(float(a[1]))
        tgfbr2.append(float(a[2]))

plt.plot(bf,cebpa,color = "#008000",linewidth=2.5)

bf = []
sox9 = []
tgfbr2 = []
cebpa = []
with open("branch_3_hep.txt") as f:
    for line in f:
        a = line[:-1].split("\t")
        bf.append(float(a[0]))
        sox9.append(float(a[3]))
        cebpa.append(float(a[1]))
        tgfbr2.append(float(a[2]))

plt.plot(bf,cebpa,c="#FF8C00",linestyle="dashed",linewidth=2.5)

bf = []
sox9 = []
tgfbr2 = []
cebpa = []
with open("branch_2_hep.txt") as f:
    for line in f:
        a = line[:-1].split("\t")
        bf.append(float(a[0]))
        sox9.append(float(a[3]))
        cebpa.append(float(a[1]))
        tgfbr2.append(float(a[2]))

plt.plot(bf,cebpa,c="#008000",linewidth=2.5)


plt.ylabel("cebpa")
plt.xlim([-0.1,50100])
plt.savefig("cebpa.png")
plt.show()