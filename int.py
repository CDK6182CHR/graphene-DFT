"""
解析计算和赝势有关的积分
"""

from sympy import *

"""
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
"""

t1,t2 = symbols('t1 t2',real=True,postive=True)
l1,l2 = symbols('l1 l2',real=True,integer=True)
t1=0
t2=0

f = 1/sqrt((t1+l1)**2+(t2+l2)**2+(t1+l1)*(t2+l2))

print(1)
s1 = summation(f,(l1,-oo,+oo))
print(2)
s2 = simplify(summation(s1,(l2,-oo,+oo)))

x,y = symbols('x,y',real=True)

ss = simplify(integrate(integrate(1/sqrt(x**2+y**2+x*y),(x,-oo,+oo)),(y,-oo,+oo)))



