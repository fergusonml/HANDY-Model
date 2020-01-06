
# coding: utf-8

# In[3]:

#get_ipython().magic('matplotlib inline')
# %load ../newman/euler.py
import numpy as np
import matplotlib.pyplot as plt

apred = 3e-5   #predator's birth rate
bpred = 2e-2   #predator's death rate
cprey = 3e-2   #prey's birth rate
dprey = 2e-4   #prey's death rate

def x_prime(x):
    return (apred*y)*x - bpred*x   #predator population rate of change

def y_prime(y):
    return cprey*y-(dprey*x)*y   #prey population rate of change

#def f(x,t):
#    return -x**3 + np.sin(t)

a = 0.0           # Start of the interval
b = 1000.0          # End of the interval
N = 10000.0          # Number of steps
h = (b-a)/N       # Size of a single step
x = 100.0           # Initial condition
y = 1000.0          # Initial condition

#create the array of time values and (empty) x-values, and loop over 
#all times supplementing x each time
tpoints = np.arange(a,b,h)
xpoints = []
ypoints = []

for t in tpoints:
    xpoints.append(x)
    ypoints.append(y)
    x += h*x_prime(x)
    y += h*y_prime(y)

plt.plot(tpoints,xpoints, color = 'red' , label='Predator')
plt.plot(tpoints,ypoints, color = 'green' , label='Prey')

plt.legend(loc='upper right')

plt.title("Predator Prey Euler Solution")
plt.xlabel("Years")
plt.ylabel("Wolves or Rabbits")
plt.show()


# In[4]:

small_r = []
big_r = []

for x in xpoints:
    big_r.append(x*5)

#plt.plot(tpoints,big_r, color = 'red' , label='PredatorX5')
#plt.plot(tpoints,ypoints, color = 'green' , label='Prey')

plt.legend(loc='upper right')

plt.title("Predator Prey Euler Solution")
plt.xlabel("Years")
plt.ylabel("Wolves or Rabbits")
#plt.show()


# In[ ]:

#parameters

alph_m = 1.0*10**(-2)
alph_max = 7.0*10**(-2)
beta_c = 3.0*10**(-2)
beta_e = 3.0*10**(-2)
s = 5.0*10**(-4)
rho = 5.0*10**(-3)
gam = 1.0*10**(-2)
lam = 1.0*10**(2)
kap = 10 # 1, 10, 100
del_ = 1#No typical value

#initial values

chi_cnot = 1.0*10**2
chi_enot = 1 # 0, 1, 25
chi_c = chi_cnot  #I'm not quite sure how to deal with this specific condition
chi_e = chi_enot
y = lam
w = 0 

#HANDY EQs

def xp_c(chi_c):
    return (beta_c - alph_c)*chi_c
def xp_e(chi_e):
    return (beta_e - alph_e)*chi_e
def yp(y):
    return (gam*(lam-y)-del_)*y
def wp(w):
    return del_*chi_c*y-c_c-c_e  #this function doesn't have the variable "w" in it explicitly: FIX

#variables and equations

w_th = rho*(chi_c + kap*chi_e)
omega = w/(w_th)
c_c = min(1, omega)*s*chi_c  #TO DO: figure out what min() means...
c_e = min(1, omega)*kap*s*chi_e  #TO DO: ^same
alph_c = alph_m + max(0, 1-c_c/(s*chi_c))*(alph_max - alph_m)  #TO DO: ^same
alph_e = alph_m + max(0, 1-c_e/(s*chi_e))*(alph_max - alph_m)  #TO DO: ^same
eta = (alph_max - beta_c)/(alph_max - alph_m)
chi = gam/(del_)*(lam-eta*s/del_)
egal_opt = 2*eta*s/(lam)
chi_m = gam*lam/(2*egal_opt) #or (gam*(lam/2)^2)/(eta*s)
phi = chi_enot/(chi_cnot)
equi_opt = 2*eta*s*(1 + phi)/(lam)
psy = chi_e/(chi_c)
uneq_opt = 2*eta*s*(1 + kap*psy)/(lam)

#Euler Solver

a = 0.0           # Start of the interval
b = 1000.0          # End of the interval
N = 10000.0          # Number of steps
h = (b-a)/N       # Size of a single step
x = 100.0           # Initial condition
y = 1000.0          # Initial condition

#create the array of time values and (empty) x-values, and loop over 
#all times supplementing x each time
tpoints = np.arange(a,b,h)
chi_cpts = []
chi_epts = []
ypts = []
wpts = []

for t in tpoints:
    chi_cpts.append(chi_c)
    chi_epts.append(chi_e)
    ypts.append(y)
    wpts.append(w)
    chi_c += h*xp_c(chi_c)
    chi_e += h*xp_e(chi_e)
    y += h*yp(y)
    w += h*wp(w)
    

#Graphing    

plt.plot(tpoints,chi_cpts, color = 'brown' , label='Commoners')
plt.plot(tpoints,chi_epts, color = 'blue' , label='Elite')
plt.plot(tpoints,ypts, color = 'green' , label='Nature')
plt.plot(tpoints,wpts, color = 'yellow' , label='Wealth')

plt.legend(loc='upper right')

plt.title("HANDY Model Simulation")
plt.xlabel("Years")
plt.ylabel("Chi & Lambda")
plt.show()

