import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

bf = []
sox9 = []
tgfbr2 = []
cebpa = []
with open("branch_1_core.txt") as f:
    for line in f:
        a = line[:-1].split("\t")
        bf.append(float(a[0]))
        sox9.append(float(a[3]))
        cebpa.append(float(a[1]))
        tgfbr2.append(float(a[2]))

plt.plot(bf,sox9,c="#0000FF",linewidth=2.5)

bf = []
sox9 = []
tgfbr2 = []
cebpa = []
with open("branch_2_core.txt") as f:
    for line in f:
        a = line[:-1].split("\t")
        bf.append(float(a[0]))
        sox9.append(float(a[3]))
        cebpa.append(float(a[1]))
        tgfbr2.append(float(a[2]))

plt.plot(bf,sox9,"r--",linewidth=2.5)

bf = []
sox9 = []
tgfbr2 = []
cebpa = []
with open("branch_3_core.txt") as f:
    for line in f:
        a = line[:-1].split("\t")
        bf.append(float(a[0]))
        sox9.append(float(a[3]))
        cebpa.append(float(a[1]))
        tgfbr2.append(float(a[2]))

plt.plot(bf,sox9,c="#0000FF",linewidth=2.5)


bf = []
sox9 = []
tgfbr2 = []
cebpa = []
with open("low_sox9_to_cebpa_link.txt") as f:
    for line in f:
        a = line[:-1].split("\t")
        bf.append(float(a[0]))
        sox9.append(float(a[3]))
        cebpa.append(float(a[1]))
        tgfbr2.append(float(a[2]))

plt.plot(bf,sox9,c="green",linewidth=2.5)

bf = []
sox9 = []
tgfbr2 = []
cebpa = []
with open("high_sox9_to_cebpa_link.txt") as f:
    for line in f:
        a = line[:-1].split("\t")
        bf.append(float(a[0]))
        sox9.append(float(a[3]))
        cebpa.append(float(a[1]))
        tgfbr2.append(float(a[2]))

plt.plot(bf,sox9,c="orange",linewidth=2.5)

plt.ylabel("sox9")
plt.xlim([-0.1,50100])
plt.savefig("sox9.png")
plt.show()