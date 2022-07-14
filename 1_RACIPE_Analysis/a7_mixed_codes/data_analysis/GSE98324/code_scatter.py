import scipy.stats
import seaborn as sns
import matplotlib.pyplot as plt

cebpa = []
tgfbr2 = []
sox9 = []

sns.regplot(tgfbr2,sox9)
plt.xlim([2.1,6.8])
plt.xlabel("tgfbr2")
plt.ylabel("sox9")
#plt.savefig("t_s.png")
plt.show()

x,y = scipy.stats.spearmanr(sox9,tgfbr2)
print("sox9","tgfbr2",x,y)
x,y = scipy.stats.spearmanr(sox9,cebpa)
print("sox9","cebpa",x,y)
x,y = scipy.stats.spearmanr(cebpa,tgfbr2)
print("cebpa","tgfbr2",x,y)
