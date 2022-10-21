#initial equations that are needed: 
#7.14 - phi** + 3Hphi* + v = 0
#7.8 - rho = 1/2 phi*^2 + v
#Friedman - H^2 = 8piG rho/3
import math 
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.integrate import cumtrapz

#def second_deg(t, y): return([-3*6.67e-11**0.5*(8*3.14 * (6.67e-11)/3*(1/2 *y[0]+ y[1]**2/(2*6.67e-11)))**(0.5) *y[0] - y[1]/6.67e-11,y[0] ])

def second_deg(t, y): return([-3*(8*3.14 * (6.67e-11)/3*(1/2 *y[0]*6.67e-11+ y[1]**2/(2)))**(0.5) *y[0]/6.67e-11**0.5 - y[1]/6.67e-11 , y[0]])
#def second_deg(t, y): return([-3*(8*3.14 * (6.67e-11)/3*(1/2 *y[0]*6.67e-11+ y[1]**2/(2)))**(0.5) *y[0]/6.67e-10**0.5 - y[1]/6.67e-10 , y[0]])
#def second_deg(t, y): return([- y[1]/(6.67e-11),y[0]/6.67e-11**0.5])

#def second_deg(t, y): return([-3/6.67e-11*(8*3.14 * (6.67e-11)/3*(1/2 *y[0]+ y[1]**2 /(2*6.67e-11)))**(0.5) *y[0] - y[1]/6.67e-11**0.5 ,y[0]])

sol = solve_ivp(second_deg,[0,50*(6.67e-11)**(1/2)], [0,5/(6.67e-11)**(1/2)])
#print(sol.t)
#print(sol.y)
#print(len(sol.t)-1)
f = np.zeros(len(sol.t))
for i in range (len(sol.t)-1):
    f[i] += 0.5 * (sol.t[i+1] - sol.t[i]) * ((sol.y[1][i]**2 + sol.y[1][i+1]**2))
   
#print(f) 

h =  (8*3.14 * (6.67e-11)/3*(1/2 *sol.y[0] + sol.y[1]**2 /(2* 6.67e-11)))**0.5
#h = 6.67e-11**(-1)*(8*3.14 * (6.67e-11)/3*(1/2 *sol.y[0] * 6.67e-11 + 6.67e-11 **2 * sol.y[1]**2 /(2* 6.67e-11)))**0.5
#print(h)
v = 1/2 * sol.y[1]**2 / (6.67e-11)
#v = sol.y[1]**2
#print(v)

tf = 0
final = 200
for i in range(len(sol.t)):
    if (sol.y[1][i] < 0 and final > i):
       tf = sol.y[1][i]
       final = i

derH = np.gradient(h*6.67e-11)

#hInt = (0.5)*(8 * 3.14 * 6.67e-11 /3 *  (1/2 * sol.y[1] + f /(2*6.67e-11))) ** (-0.5)
hInt = cumtrapz(h,sol.t)
#print(hInt[160])
a = np.zeros(len(hInt))
for i in range(len(hInt)):
    a[i] =  math.exp(hInt[i]-240) 
t = np.zeros(183)
derHf = np.zeros(183)
hf = np.zeros(183)
#print(derH)
#print(a)
#print(h)
for i in range(183):
   t[i] += sol.t[i]
   derHf[i] += derH[i] 
   hf[i] += h[i]

eSlow = -derHf*6.67e-11**(5/2)/((a*((hf**2))))
pR = 2*math.pi *1/(20000*6.67e-11**0.5)* (hf/3e8)**2/(0.05**3*eSlow)
#2.1e-9
aS = 0.05**3/(2*math.pi**2)*pR
#print(aS)
#0.9682
ns = 1-4*eSlow
#print(ns)
aInfl = np.zeros(161)
hInfl = np.zeros(161)
eSlowInfl = np.zeros(161)
nsInfl = np.zeros(161)
asInfl = np.zeros(161)
for i in range(161):
    aInfl[i] += a[i]
    hInfl[i] += h[i]
    eSlowInfl[i] += eSlow[i]
    nsInfl[i] += ns[i]
    asInfl[i] += aS[i]
r = 16*eSlowInfl
#t = np.zeros(183)
#print(eSlowInfl[-1])
#print(nsInfl[-1])
#print(asInfl[-1])
#print(r[-1])
#plt.plot(a,hf)
#plt.plot(sol.t, sol.y[1])
#plt.plot(aInfl,hInfl)
#plt.plot(aInfl*hInfl, eSlowInfl)
#plt.plot(a,a*hf)
#plt.yscale("log")
#plt.xscale("log")
#plt.plot(sol.y[1] * 6.67e-11 **0.5,v * 6.67e-11 **2)
#plt.show()

#def second_deg(t,y): return
