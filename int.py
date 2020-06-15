"""
解析计算和赝势有关的积分
"""

from sympy import *

u,v = symbols('u v',real=True,postive=True)
Kx,Ky=symbols('Kx,Ky',real=True)
a,c,d = symbols('a c d',real=True)

x = a/sqrt(3)*(u+v)
y = a*(u-v)

fuv = exp(-I*(Kx*(x-c)+Ky*(y-d)))*exp(-sqrt((x-c)**2+(y-d)**2)/a)*2*a**2/sqrt(3)

print(1)
s1 = integrate(fuv,(u,0,1))
print(2)
s2 = simplify(integrate(s1,(v,0,1)))

print(s2)


