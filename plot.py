import matplotlib.pyplot as plt

eig_file = r'EIGENVAL.txt'
path_file = r'KPATH.txt'

K_names = ['K','G','M','K']
#K_names = ['G','M','K','G']
excludes = [6,14]

excludes = []

e=1.602176487e-19
E0 = -5477  # 零点能量，按电子伏特

with open(eig_file,'r') as fp:
    data = list(map(lambda line:list(map(float,line.split())),
                    filter(bool,fp.readlines())))
trans = list(zip(*data))

ks = list(trans[0])
print(len(ks))
print("Excludes: ",excludes)

Es = list(map(list,trans[1:]))

for i in sorted(excludes,reverse=True):
    del ks[i]
    for E in Es:
        del E[i]




with open(path_file,'r') as fp:
    k_dots = list(map(lambda line:float(line.split()[0]),
                      filter(bool,fp.readlines())))

assert len(k_dots)==len(K_names)

for E in Es:
    for i in range(len(E)):
        E[i]=E[i]/e-E0

Emin=min(map(min,Es))
Emax=max(map(max,Es))

for k in k_dots:
    plt.plot([k,k],[Emin-100,Emax],'k')
plt.plot([k_dots[0],k_dots[-1]],[0,0],'--')
#plt.ylim((-100,200))

for E in Es:
    plt.plot(ks,[Ei for Ei in E],'b')
    plt.plot(ks,[Ei for Ei in E],'b.')

plt.xticks(ticks=k_dots,labels=K_names)
plt.ylabel('E (eV)')
plt.show()
