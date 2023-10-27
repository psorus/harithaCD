import numpy as np
import pandas as pd
from plt import *

from draw import draw_cd
from stats import wilcoxon_holm


data=pd.read_csv("example.csv",index_col=False)
pval, avg_rank,_=wilcoxon_holm(df_perf=data)


plt.figure(figsize=(10, 7))

draw_cd(pval,avg_rank)

plt.savefig("last.png")

plt.show()


