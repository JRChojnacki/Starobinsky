# -*- coding: utf-8 -*-
"""
Created on Sat Feb 29 21:34:31 2020
@author: Janek
"""
import matplotlib.pyplot as plt
import numpy as np
import Gauss
from openpyxl import load_workbook
t0=0                #chwila początkowa                  #wartosc poczatkowa fi
M=1# 2.435 * 10**18
V0 = 8*10**(-11)*M**4
f0=6 *M
N=10500   #liczba iteracji w równaniu
              #Stała częsc potencjału
X=[]
              #Masa plancka
def V(f):
    return(V0*(1-np.exp(-np.sqrt(2/3) *f/M))**2)
def epsilon(dV,V):
    return(0.5*M**2*(dV/V)**2)
def Hubble(V0,M,a,p,f):         #ewoucja parametru hubble'a
    return(1/(np.sqrt(3)*M)*np.sqrt(np.abs(V0*(1-np.exp(-np.sqrt(2/3)*f/M))**2)))
def step(H):                    #ewolucja korku czasowego
    return(1/(100*H))
def df(f):                #pochodna potencjalu  
    return(2*V0*(1-np.exp(-np.sqrt(2/3)*f/M))*np.sqrt(2/3)/M*np.exp(-np.sqrt(2/3)*f/M))  
def solve(f0,n,H0,X):       #rozwiązanie numeryczne 
    global E,t,laststep
    f=np.zeros(n+1)
    t=np.zeros(n+1)
    H=np.zeros(n+1)
    fn=np.zeros(n+1)
    E=np.zeros(n+1)
    f[0]=f0
    fn[0]=f0
    H[0]=H0
    for i in range(n):
        if E[i]>1:
            f=f[:i]
            t=t[:i]
            H=H[:i]
            fn=fn[:i]
            E=E[:i]
            laststep=i
            break 
        deviation=np.sqrt(((H[i])*32)/((2*np.pi)**2)*step(H[i]))
        s=Gauss.fluctuation(0,deviation)
        f[i+1]=f[i]-1/(3*H[i])*df(f[i])*step(H[i])+s
        fn[i+1]=fn[i]-1/(3*H[i])*df(fn[i])*step(H[i])        
        t[i+1]=t[i]+step(H[i])
        H[i+1]=H0#Hubble(V0,M,a,p,fn[i])
        v=V(fn[i])
        dv=df(fn[i])
        E[i+1]=epsilon(dv,v)
        if H[i+1]>np.abs(f[i+1]):           
            X.append(H[i+1]*t[i+1])
    fig, axs=plt.subplots(2,sharex=True)
    axs[0].plot(t*H0,f, color='black', linestyle='-', linewidth=1,label='fluctuating')
    axs[0].plot(t*H0,fn, color='green', linestyle='-', linewidth=1,label='slow-roll')   
    axs[0].set_title("Exemplary field evolution, $\phi_0$={}".format(f0))
    axs[0].set(ylabel='$\phi$')
    #axs[0].set_ylim(0,2*f0)
    axs[1].plot(t*H0,E,color='blue')
    #axs[1].set_ylim(E[2],np.max(E))
    axs[1].set_title("Slow-roll $\epsilon$ parameter evolution")
    axs[1].set(xlabel='$Ht$')
    axs[1].hlines(1,0,laststep/100,colors='red',linestyle='--')
    axs[1].set(ylabel='$\epsilon$')
    fig.legend()
    for ax in axs.flat:
        ax.label_outer()
    #plt.xlabel('Ht')axs[0].set_xlim(0, 2)
    #plt.ylabel('$\phi$')
    #plt.grid(True)
    plt.savefig('Samplef0{}.png'.format(f0))
    return[X]
H0=np.sqrt(np.abs(V(f0)/3))/M
solve(f0,N,H0,X)