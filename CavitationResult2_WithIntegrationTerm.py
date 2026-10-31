import matplotlib.pyplot as plt
import math
from scipy.interpolate import interp1d
from scipy import integrate
import math
import numpy as np
from sympy import *
import mpmath as mp
from decimal import *
Ro=10**-4
d=998.2
pV=0
pInf=10350
pGo=10**8
s=1
n=10
t=0.0
h=10**-7
i=0
k=[0,0,0,0]
l=[0,0,0,0]
Radi1=[]
time1=[]
R_dot=[]
Radi2=[]
time2=[]
Radi3=[]
time3=[]

lam=1000
def funcA(R,z):
    #f = (((pV-pInf)/(d*R))-((2*s/d)/R**2)-(4*n/d)*(z/(R**2))-((1.5*z*z)/R))*10**-6
    f = (((pV - pInf) / (d * R)) + (pGo * ((Ro / R) ** 3)) / (R * d) - ((2 * s / d) / R ** 2) - (4 * n / d) * (z / (R ** 2)) - ((1.5 * z * z) / R))
    return f

R = Ro
z = 0.0
Radi1.append(Ro)
time1.append(0.0)
R_dot.append(z)
for i in range(1,100000,1):
    l[0] = h * funcA(R, z)
    k[0] = h * z
    l[1] = h * funcA((R + (k[0] / 2.0)), (z + (l[0] / 2.0)))
    k[1] = h * (z + (l[0] / 2.0))
    l[2] = h * funcA((R + (k[1] / 2.0)), (z + (l[1] / 2.0)))
    k[2] = h * (z + (l[1] / 2.0))
    t = t + h
    l[3] = h * funcA(R + k[2], (z + l[2]))
    k[3] = h * (z + l[2])
    R = R + ((k[0] + 2.0 * k[1] + 2.0 * k[2] + k[3]) / 6.0)
    z = z + ((l[0] + 2.0 * l[1] + 2.0 * l[2] + l[3]) / 6.0)
    Radi1.append(R)
    time1.append(t)
    R_dot.append(z)
plt.plot(time1,Radi1,'b-')
plt.xlabel("time")
plt.ylabel("Radius")
plt.show()
r = interp1d(time1,Radi1,kind='cubic')
r1 = interp1d(time1,R_dot,kind='cubic')

def function_A(ti,T):
    if(float(r(ti)-r(T))==0.0):
        return 0.0
    else:
        expr=((r1(ti) * r(ti)**2) / ((r(ti) ** 3) - (r(T) ** 3)))*(mp.exp((-(T - ti) / lam)))*math.log(abs(r(ti)/r(T)))
        return expr

def Integrate(T):
    sum = 0.0
    y=[]
    if(T==0):
        return sum
    else:
        H = t / 1000
        accu_1 = 0
        for x in range(1,1000):
            f_v_1 = function_A(x * H,T)
            accu_1 = accu_1 + (2 * f_v_1) if (x % 2 == 0) else (4 * f_v_1)

        accu_1 = (accu_1 * H) / 3.0
        return accu_1

def funcB(R,z,t):
    #fi = (((pV - pInf) / (d * R)) + (pGo * ((Ro / R) ** 3)) / (R * d) - ((2 * s / d) / R ** 2) - (4 * n / d) * (z / (R ** 2)) - ((1.5 * z * z) / R))-(12.0/(R*d))*Integrate(t)
    fi=funcA(R,z)-(12.0/(R*d))*Integrate(t)
    return fi

R=Ro
z=0.0
t=0.0
Radi2.append(R)
time2.append(t)
for i in range(1,100000,1):
    l[0] = h * funcB(R,z,t)
    k[0] = h * z
    l[1] = h * funcB((R + (k[0] / 2.0)),(z + (l[0] / 2.0)),t+(h/2.0))
    k[1] = h * (z + (l[0] / 2.0))
    l[2] = h * funcB((R + (k[1] / 2.0)), (z + (l[1] / 2.0)),(t+h/2.0))
    k[2] = h * (z + (l[1] / 2.0))
    t = t + h
    l[3] = h * funcB(R + k[2], (z + l[2]),(t))
    k[3] = h * (z + l[2])
    R = R + ((k[0] + 2.0 * k[1] + 2.0 * k[2] + k[3]) / 6.0)
    z = z + (l[0] + 2.0 * l[1] + 2.0 * l[2] + l[3]) / 6.0
    Radi2.append(R)
    time2.append(t)



plt.plot(time2,Radi2,'g-')
plt.xlabel("time")
plt.ylabel("Radius")
plt.show()

sum=0.0
newRadi=[]
newRadi.append(0.0)
for x in range(1,10000,1):
    sum=sum+(Radi2[x]-Radi1[x])**2
    newRadi.append(math.sqrt(sum))

plt.plot(time1,newRadi,'r-')
plt.xlabel("time in seconds")
plt.ylabel("Radius in meters")
plt.show()