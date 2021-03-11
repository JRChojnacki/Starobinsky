# -*- coding: utf-8 -*-
"""
Created on Sat Feb 29 21:34:31 2020
@author: Janek
"""
import matplotlib.pyplot as plt    #Rysuneczki
import numpy as np   #Podstawowe operacje matematyczne, głównie wektory, macierze
#import Gauss          #krótki moduł losujący liczbę z rozkładu Gaussa   
#from openpyxl import load_workbook  #operacja na plikach w excelu
t0=0                #warunek początkowy czasu w równaniach różniczkowych
M= 1#2.435 * 10**18 #Masa Plancka, przyjęta za 1 oraz jej wartosć w GeV   
f0=15*M #warunek poczatkowy na wartosc pola w równaniach różniczkowych
V0 = 8*10**(-11)*M**4  #Wartoćc parametru V0 w modelu Starobinskiego, wzieta przez Rudeliusa z powerspectrum
N=10500   #liczba iteracji w równaniu, praktycznie dziala, jak czas przez jaki puszczamy ewolucje pola
X=[]          # baza histogramu, zmienna w którą wrzucam wartosć H*t, dla spełnionego warunku inflacji \phi>1
laststep=0
def V(f):           #potencjał
    return(V0*(1-np.exp(-np.sqrt(2/3) *f/M))**2)
def epsilon(dV,V):  #inflacyjny parametr epsilon<<1 => inflacja zachodzi
    return(0.5*(M*dV/V)**2)
def Hubble(V0,M,a,p,f):         #ewoucja parametru hubble'a, nie wykorzystana w programie
    return(1/(np.sqrt(3)*M)*np.sqrt(np.abs(V0*(1-np.exp(-np.sqrt(2/3)*f/M))**2)))
def fluctuation(mean,deviation):
    mu, sigma = 0, 0.1 # mean and standard deviation for quantum fluctuation
    s = np.random.normal(mu, deviation, 1)
    return(s)
def step(H):                    #ewolucja korku czasowego
    return(1/(100*H))
def df(f):                #pochodna potencjalu  
    return(2*V0*(1-np.exp(-np.sqrt(2/3)*f/M))*np.sqrt(2/3)/M * np.exp(-np.sqrt(2/3)*f/M))  
def solve(f0,n,H0,X):       #rozwiązanie numeryczne
    global En,laststep, t  # w E wsadzam wartosci epsilona, w t czas
    f=np.zeros(n+1)  #wartosć pola Langevin
    t=np.zeros(n+1) 
    H=np.zeros(n+1)  #Wartosc parametru Hubble'a
    fn=np.zeros(n+1) #wartosc pola slow-roll
    E=np.zeros(n+1)
    En=np.zeros(n+1)
    f[0]=f0
    fn[0]=f0
    H[0]=H0
    for i in range(n):
        if E[i]>= 1:
            f=f[:i]
            t=t[:i]
            H=H[:i]
            fn=fn[:i]
            E=E[:i]
            laststep+=i
            #print(laststep)
            break 
        deviation=np.sqrt(((H[i])*32)/((2*np.pi)**2)*step(H[i])) #parametr sigma w gaussie
        s=fluctuation(0,deviation) #fluktuacja- losowa liczba z rozkladu gaussa
        f[i+1]=f[i]-1/(3*H[i])*df(f[i])*step(H[i])+s #zdyskretyzowany Langevin
        fn[i+1]=fn[i]-1/(3*H[i])*df(fn[i])*step(H[i])        #zdyskretyzowany Slow-roll
        t[i+1]=t[i]+step(H[i])
        H[i+1]=H0#Hubble(V0,M,a,p,fn[i])
        v=V(f[i]) 
        dv=df(f[i])
        E[i+1]=epsilon(dv,v)
        X.append(H[i+1]*t[i+1])          
    return[X]
H0=np.sqrt(np.abs(V(f0)/3))/M #Wartosc parametru Hubble, wynika z drugiego rownania slow-rollu
for i in range(10000):      #10000 krotna ewolucja
    solve(f0,N,H0,X)
    #print(i/100)        #procent całosci
    #print(' ')
B=[]
C=[]
counts,bins=np.histogram(X,1000,density=True) # counts to liczba zliczeń w danym przedziale czasowym bins
#plt.hist(X,200,range=(0,6),density=1) #rysuje histogram
for i in range(1000):
    if counts[i]!=0: #bez tego logarytm wybucha
        C.append(-np.log(counts[i])) #
        B.append(bins[i])
lastbin=int(laststep/(N*10))
CL=C[:lastbin]
BL=B[:lastbin]
alpha,betha=np.polyfit(BL, CL, 1)
Y=alpha*np.array(B)+betha
#fig, axs=plt.subplots(2,sharex=True)
plt.plot(B,3*H0*np.array(B),linestyle='--',color='red',label='EI boundary')
plt.plot(B,Y,linewidth=2.0,color='green',label="Linear fit")
plt.plot(B, C, linewidth=1.0, color='black',label='Probability')   
plt.title(r'Probability that inflation occurs for the Starobinsky model, $\phi_0$={}'.format(f0))
plt.ylabel('-log(P)')
plt.xlabel('$Ht$')
plt.legend()

#axs[1].plot(H0*t,En/10000,color='blue')
#axs[1].hlines(1,0,laststep/100,colors='red',linestyle='--')
#axs[1].set_ylim(En[2],None)
#axs[1].set_title("Slow-roll $\epsilon$ parameter evolution")
#axs[1].set(xlabel='$Ht$')
#axs[1].set(ylabel='$\epsilon$')
#fig.legend()
plt.savefig('ProbabilityStaryf0={}.png'.format(f0))
print("f0= ",f0,"alpha= ",alpha,'3H0= ',3*H0,"avrg. inflation time= ", laststep/100)
#for ax in axs.flat:
#    ax.label_outer()
#new_row_data = [
#    ["f0= ",f0,"alpha= ",alpha,'3H0= ',3*H0,"avrg. inflation time= ", laststep/100],
#    ]
#wb = load_workbook("Stary.xlsx")
## Select First Worksheet
#ws = wb.worksheets[0]
## Append 2 new Rows - Columns A - D
#for row_data in new_row_data:
#    # Append Row Values
#    ws.append(row_data)
#wb.save("Stary.xlsx")