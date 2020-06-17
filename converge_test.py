"""
测试晶格势能的收敛情况
"""
from math import *
import matplotlib.pyplot as plt

t1 = 0.8
t2 = 0.8

def count(cut,t1,t2):
    s = 0
    for i in range(-cut,cut+1):
        for j in range(-cut,cut+1):
            s+=1.0/sqrt((i+t1)**2+(j+t2)**2+(i+t1)*(j+t2))
    return s

ns = []
cs = []
cuts = list(range(10))
for cut in cuts:
    c = count(cut,t1,t2)
    n = (2*cut+1)**2
    ns.append(n)
    cs.append(c)
    print(cut,n,c,end=' ')
    if cut:
        print((c-cs[0])/cut)
    else:
        print('')

plt.plot(ns,cs)
plt.show()

plt.plot(cuts,cs)
plt.show()
