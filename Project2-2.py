# -*- coding: utf-8 -*-
"""
Created on Tue Mar 01 15:36:28 2016
@author: Hp
"""
import scipy
'Fluids    HOT FLUID       COLD FLUID'                                               
'Fluids    Waste Water     Cooling Water'
Qh=11712.4; Qc=11698.4;   #kW
delTlm=23;  #oC
Fh=140;muc=0.000966;muh=0.000509;Prh=3.31;Prc=5.19;Kh=0.645;Kc=0.617;Kp=17.5;
Lvp=1.55;Dp=0.2;Lhp=0.43;Lw=0.63;phi=1.25;Lc=0.38;t=0.0006;Np=1;
#Ua=5000 # Assumption (W/m2.K)
e=100;m=0;n=100;
print Fh
while abs(e)>1:
    Ua=(5000+m*0.1)
    
    
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
#print A
#print Nt
#print Ncp
#print Hh
#print Hc
#print Uc


      