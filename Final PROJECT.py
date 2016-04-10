# -*- coding: utf-8 -*-
"""
Created on Tue Mar 01 15:36:28 2016
@author: Hp
"""
import scipy
from scipy import array
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy import linspace
'Fluids    HOT FLUID       COLD FLUID'                                               
'Fluids    Waste Water     Cooling Water'
Qh=11712.4; Qc=11698.4;   #kW
delTlm=23;  #oC
Fh=140;muc=0.000966;muh=0.000509;Prh=3.31;Prc=5.19;Kh=0.645;Kc=0.617;Kp=17.5;
Lvp=1.55;Dp=0.2;Lhp=0.43;Lw=0.63;phi=1.25;Lc=0.38;t=0.0006;Np=1;

e=100;m=0;n=100;rhoh=985;rhoc=995;Cph=4183;qmax=1.3959;k=1.002;n=5;
Q=array([1.0,2.0,3.0,4.0,5.0])
print Fh
while abs(e)>10:
    Ua=(5000+m*1)
    A=(Qh*1000/(Ua*delTlm))
    Lp=(Lvp-Dp)
    Ap=Lp*Lw
    Asp=phi*Ap
    Ne=(A/Asp)
    Nt=Ne+2
    p=(Lc/Nt)
    b=(p-t)
    Ach=b*Lw
    Dh=(2*b/phi)
    Ncp=((Nt-1)/(2*Np))
    mch=(Fh/Ncp)
    Gch=(mch/Ach)
    Reh=(Gch*Dh/muh)
    Rec=(Gch*Dh/muc)
    Nuh=0.3*(Reh)**(0.663)*(Prh)**(0.333)
    Hh=(Nuh*Kh/Dh)
    Nuc=0.3*(Rec)**(0.663)*(Prc)**(0.333)
    Hc=(Nuc*Kc/Dh)
    Uc=((1/Hc)+(1/Hh)+(t/Kp))**(-1)
    e=abs(Uc-Ua)
    print e
    m=m+1
print e
print Ua
print Uc
Gp=(Fh)/(3.142*0.25*Dp**2)
fh=(1.441/Reh**0.206)
fc=(1.441/Rec**0.206)
####
'Assumption=First all the water will flow in the X-direction for simplify the compicated eq of vector + PDE'
"CASE 1"
' Water Flow in the X-direction.' 
'Component balace n Darcy-Forchheimer eq had solved for getting Pressure and Velocity Profile along X-direction'
def epid(y,x):
    rho=985
    mu=0.000509
    k=10**-2
    B=(mu/k)
    dy1dx=rho*y[0]*y[1]+B*y[0]
    dy0dx=y[1]
    return (dy1dx,dy0dx)
    
x=linspace(0,0.63,5)
yinitial=[1.3959,11]
solx=odeint(epid,yinitial,x)
print solx 
a=[solx[0,0],solx[2,0],solx[3,0],solx[4,0]]
print a
def rock(y,x):
    rho=985
    mu=0.000509
    k=5.14
    f=rho*9.81+y[0]*((mu/k))
    return f
    
x=linspace(0,0.63,5)
y0=1.3959
ansx=odeint(rock,y0,x)
print ansx
print "Graph of DelP vs x"
fig=plt.figure()
ax=fig.add_subplot(111)
plt.plot(x,ansx,'r')
ax.title.set_text('DelP Vs x')
ax.xaxis.label.set_text('x')
ax.yaxis.label.set_text('DelP')
plt.show()
'From graph of DelP vs x we can see that Pressure drop increase as we move forward in X-direction that was as expected..'
print "Graph of DelP vs Inlet Velocity(Q)"
fig=plt.figure()
ax=fig.add_subplot(111)
plt.plot(Q,ansx,'r')
ax.title.set_text('DelP Vs Q')
ax.xaxis.label.set_text('Q')
ax.yaxis.label.set_text('DelP')
plt.show()
'From graph of DelP vs Q(Inlet Velocity), we can see that as Velocity increase DelP increase that was as expected.'
"CASE 2"
' Water Flow in the Y-direction.' 
'Component balace n Darcy-Forchheimer eq had solved for getting Pressure and Velocity Profile along Y-direction'
def epid(z,y):
    rho=985
    mu=0.000509
    k=10**-2
    B=(mu/k)
    dz1dy=rho*z[0]*z[1]+B*z[0]
    dz0dy=z[1]
    return (dz1dy,dz0dy)
    
y=linspace(0,1.55,5)
zinitial=[1.3959,11]
soly=odeint(epid,zinitial,y)
print soly 
a=[soly[0,0],soly[2,0],soly[3,0],soly[4,0]]
print a
def rock(z,y):
    rho=985
    mu=0.000509
    k=5.14
    f=rho*9.81+z[0]*((mu/k))
    return f
    
y=linspace(0,1.55,5)
z0=1.3959
ansy=odeint(rock,z0,y)
print ansy
print "Graph of DelP vs y"
fig=plt.figure()
ax=fig.add_subplot(111)
plt.plot(y,ansy,'r')
ax.title.set_text('DelP Vs y')
ax.xaxis.label.set_text('y')
ax.yaxis.label.set_text('DelP')
plt.show()
'From graph of DelP vs x we can see that Pressure drop increase as we move forward in Y-direction that was as expected..'
print "Graph of DelP vs Inlet Velocity(Q)"
fig=plt.figure()
ax=fig.add_subplot(111)
plt.plot(Q,ansy,'r')
ax.title.set_text('DelP Vs Q')
ax.xaxis.label.set_text('Q')
ax.yaxis.label.set_text('DelP')
plt.show()

'''
From the above calculations and plots,
we can see that result that I am getting was as expected.
DelP increases with velocity in both direction,and also with length.
'''
