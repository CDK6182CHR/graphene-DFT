import matplotlib.pyplot as plt

eig_file = r'EIGENVAL.txt'
path_file = r'KPATH.txt'

K_names = ['K','G','M','K']

e=1.602e-19

with open(eig_file,'r') as fp:
    data = list(map(lambda line:list(map(float,line.split())),
                    filter(bool,fp.readlines())))
trans = list(zip(*data))

ks = trans[0]
Es = trans[1:]
Emin=min(map(min,Es))
Emax=max(map(max,Es))

with open(path_file,'r') as fp:
    k_dots = list(map(lambda line:float(line.split()[0]),
                      filter(bool,fp.readlines())))

assert len(k_dots)==len(K_names)

for E in Es:
    plt.plot(ks,[Ei/e for Ei in E],'b')

for k in k_dots:
    plt.plot([k,k],[Emin/e,Emax/e])
plt.plot([k_dots[0],k_dots[-1]],[0,0])

plt.xticks(ticks=k_dots,labels=K_names)
plt.show()
