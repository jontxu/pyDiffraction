# -*- coding: utf-8 -*- 
# doubleslitdiffraction.py
#     Diffraction + interference pattern for the double-slit experiment

from pylab import *
from math import *
from numpy import swapaxes, transpose, concatenate, sinc, around
import matplotlib.pyplot as plt

lamb = 0.7 * 10**(-6)   # wavelength (m)  500 nm
k=2.0*pi/lamb
a = 0.05*10**(-3)        # slit width  (m)
d = a * 5.0
dmm = d * 1000.0        # slit width  (mm)
amm = a * 1000.0
lambnm=lamb*10**9       # wavelength  (nm)

lamb2 = 0.5 * 10**(-6)   # wavelength (m)  500 nm
k2=2.0*pi/lamb2
a2 = 0.05*10**(-3)        # slit width  (m)
d2 = a * 5.0
dmm2 = d2 * 1000.0        # slit width  (mm)
amm2 = a2 * 1000.0
lambnm2=lamb2*10**9       # wavelength  (nm)


num=1000		                 # number of points	
n = 500
N= 5
intenrojo=zeros(num)
thetxrojo=zeros(num)

for ithetrojo in arange(-499, 501):
    thet = ithetrojo/(50.0*n)
    beta = k*a*sin(thet)/2.0
    delta = k*d*sin(thet)/2.0
    valor = sin(N * delta)**2
    valorb = sin(delta)**2
    inten = divide(valor, valorb)
    intenrojo[ithetrojo+499] = (inten * sinc(beta/pi)**2) / N**2
    thetxrojo[ithetrojo+499] = thet
	
intenrojo2=zeros(num)
thetxrojo2=zeros(num)

for ithetrojo2 in arange(-499, 501):
    thet = ithetrojo2/(50.0*n)
    beta = k*a*sin(thet)/2.0
    valorc = sin(beta)
    valord = beta
    inten2 = divide(valorc, valord)**2
    intenrojo2[ithetrojo2+499] = sinc(beta/pi)**2
    thetxrojo2[ithetrojo2+499] = thet

intenverde=zeros(num)
thetxverde=zeros(num)
	
for ithetverde in arange(-499, 501):
    thet = ithetverde/(50.0*n)
    beta = k2*a2*sin(thet)/2.0
    delta = k2*d2*sin(thet)/2.0
    valor = sin(N * delta)**2
    valorb = sin(delta)**2
    inten = divide(valor, valorb)
    intenverde[ithetverde+499] = (inten * sinc(beta/pi)**2) / N**2
    thetxverde[ithetverde+499] = thet

intenverde2=zeros(num)
thetxverde2=zeros(num)

for ithetverde2 in arange(-499, 501):
    thet = ithetverde2/(50.0*n)
    beta = k2*a2*sin(thet)/2.0
    valorc = sin(beta)
    valord = beta
    inten2 = divide(valorc, valord)**2
    intenverde2[ithetverde2+499] = sinc(beta/pi)**2
    thetxverde2[ithetverde2+499] = thet
	
plt.plot(thetxrojo,intenrojo,'r-', label='I ' + str(lambnm) + ' nm')
plt.plot(thetxrojo2,intenrojo2,'r--', label='D ' + str(lambnm) + ' nm')
#plt.plot(thetxverde,intenverde,'g-', label='I ' + str(lambnm2) + ' nm')
#plt.plot(thetxverde2,intenverde2,'g--', label='D ' + str(lambnm2) + ' nm')

plt.legend()
axes = plt.gca()
axes.set_xlim([around(thetxverde[0],decimals=3), thetxverde[999]])
plt.title('Red de difraccion para dos rendijas de ancho ' + str(amm) + ' mm y distancia ' + str(dmm) + ' mm')
plt.xlabel('Theta')
plt.ylabel('Intensidad')
plt.grid()
plt.show()

