#   diffraction.py
#    Single Slit Diffraction

from pylab import *
from math import *
from numpy import swapaxes, transpose, concatenate, sinc, around
import matplotlib.pyplot as plt

#plt.font('serif')

lamb = 0.7 * 10**(-6)   # wavelength (m)  500 nm
k=2.0*pi/lamb
d = 0.05*10**(-3)        # slit width  (m)
dmm = d * 1000.0        # slit width  (mm)
lambnm=lamb*10**9       # wavelength  (nm)

lamb2 = 0.5 * 10**(-6)   # wavelength (m)  500 nm
k2=2.0*pi/lamb2
d2 = 0.05*10**(-3)        # slit width  (m)
dmm2 = d2 * 1000.0        # slit width  (mm)
lambnm2=lamb2*10**9       # wavelength  (nm)


num=1000		                 # number of points	
n = 500

intenrojo=zeros(num)
thetxrojo=zeros(num)

for ithetrojo in arange(-499, 501):
    thet = ithetrojo/(50.0*n)
    beta = k*d*sin(thet)/2.0
    intenrojo[ithetrojo+499] = sinc(beta/pi)**2
    thetxrojo[ithetrojo+499] = thet

intenverde=zeros(num)
thetxverde=zeros(num)
	
for ithetverde in arange(-499, 501):
    thet = ithetverde/(50.0*n)
    beta = k2*d2*sin(thet)/2.0
    intenverde[ithetverde+499] = sinc(beta/pi)**2
    thetxverde[ithetverde+499] = thet
			
plt.plot(thetxrojo,intenrojo,'r-', label=str(lambnm) + ' nm')
plt.plot(thetxverde,intenverde,'g-', label=str(lambnm2) + ' nm')

plt.legend()
axes = plt.gca()
axes.set_xlim([around(thetxverde[0],decimals=3), thetxverde[999]])
plt.title('Difraccion para una rendija larga y estrecha de ancho ' + str(dmm) + ' mm')
plt.xlabel('Theta')
plt.ylabel('Intensidad')
plt.grid()
plt.show()

