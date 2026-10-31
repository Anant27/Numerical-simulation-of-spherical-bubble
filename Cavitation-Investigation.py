import matplotlib.pyplot as plt
import math
from scipy import integrate
from scipy.interpolate import interp1d
Ro=10**-4
d=999
pV=0
pInf=10350
pGo=10**8
s=1
n=0.1
t=0.0
h=10**-7
i=0
k=[0,0,0,0]
l=[0,0,0,0]
lm=1000
Radi1=[]
time1=[]
R_dot=[]
Radi2=[]
time2=[]
R_dot2=[]

ran = 1000
def funcA(R,z):
    #f = (((pV-pInf)/(d*R))-((2*s/d)/R**2)-(4*n/d)*(z/(R**2))-((1.5*z*z)/R))*10**-6
    f = (((pV - pInf) / (d * R)) + (pGo * ((Ro / R) ** 3)) / (R * d) - ((2 * s / d) / R ** 2) - (4 * n / d) * (z / (R ** 2)) - ((1.5 * z * z) / R))
    return f

def funcB(R,z,t):
    #fi = (((pV - pInf) / (d * R)) + (pGo * ((Ro / R) ** 3)) / (R * d) - ((2 * s / d) / R ** 2) - (4 * n / d) * (z / (R ** 2)) - ((1.5 * z * z) / R))-(12.0/(R*d))*Integrate(t)
    fi=funcA(R,z)-(12.0/(R*d))*Integrate(t)
    return fi

R = Ro
z = 0.0
Radi1.append(Ro)
time1.append(0.0)
R_dot.append(z)
for i in range(0,ran+100000,1):
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

plt.plot(time1,Radi1)
plt.title('Without Integration')
plt.show()

Rad=interp1d(time1,Radi1,kind='cubic',fill_value='extrapolate')
RDOT=interp1d(time1,R_dot,kind='cubic',fill_value='extrapolate')

#def f_1(x,b):
#    if(x==b):
#        return 0.0
#    value = ((R_dot[x] * Radi1[x] / (math.pow(Radi1[x], 3) - math.pow(Radi1[b], 3))) * math.exp(-(time1[x] - time1[b]) / lm) * (math.log(Radi1[x]) - math.log(Radi1[b])))
#    return value

def f_1(t,T):
    if(t==T):
        return 0.0
    value = ((RDOT(t) * Rad(t) / (math.pow(Rad(t), 3) - math.pow(Rad(T), 3))) * math.exp(-(t - T) / lm) * (math.log(Rad(t)) - math.log(Rad(T))))
    return value

def Integrate(t):
    if(t==0.0):
        return 0.0
    H=t/1000.0
    accu_1=0
    for x in range(1,1000):
        f_v_1=f_1(x*H,t)
        accu_1 = accu_1 + (2 * f_v_1) if (x % 2 == 0) else (4 * f_v_1)

    accu_1 = (accu_1*H)/3.0
    return accu_1

R = Ro
z = 0.0
Radi2.append(Ro)
time2.append(0.0)
R_dot2.append(z)
for i in range(0,ran-2,1):
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
    R_dot2.append(t)

sum=0.0
newRadi=[]
newRadi.append(Ro)
newTime=[]
newTime.append(0.0)
for x in range(1,ran-2,1):
    sum=sum+(Radi2[x]-Radi1[x])**2
    newRadi.append(math.sqrt(sum))
    newTime.append(x*h)

plt.plot(time2,Radi2)
plt.title('With Integration')
plt.show()
plt.plot(newTime,newRadi)
plt.title('with L2 Norm')
plt.show()

