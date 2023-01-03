import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def func(x, a, b):
    return a * np.exp(-b * x) + (1-a)

#x = np.linspace(0,4,50)
#y = func(x, 2.5, 1.3, 0.5)
#yn = y + 0.2*np.random.normal(size=len(x))

data = np.loadtxt("avg_diab_pop_200000K_0.1VER_withmetal1.out")
x=data[:,0]
yn=1-data[:,2]
p=np.array([0.8,2.e-5])

popt, pcov = curve_fit(func, x, yn,p)

print("rate (ps-1)=",popt[1],(pcov[1,1])**0.5)
print("thermal pop =",popt[0],pcov[0,0]**0.5)

plt.figure()
plt.plot(x, yn, 'ko', label="Original Noised Data")
plt.plot(x, func(x, *popt), 'r-', label="Fitted Curve")
plt.legend()
plt.show()

