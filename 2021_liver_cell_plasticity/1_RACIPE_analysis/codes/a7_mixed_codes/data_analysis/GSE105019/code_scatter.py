import scipy.stats
import seaborn as sns
import matplotlib.pyplot as plt

cebpa = [3.514452396,
6.102356371,
3.459519076,
3.539586865,
0.49596141,
0.104080545,
5.600338678,
4.376641277,
2.660454943,
0.419740242,
0.148002261,
1.699024629]
tgfbr2 = [4.890445138,
2.336883716,
5.401322466,
4.606394297,
5.994753954,
6.710816817,
2.468598638,
5.37165164,
4.568183351,
6.151971281,
6.22591548,
5.311058861]
sox9 = [3.044870513,
1.58832308,
2.751756223,
1.340136391,
2.666489367,
3.044552291,
0.56153584,
2.954490214,
1.624370136,
3.311877041,
3.246344615,
4.441880061]

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
