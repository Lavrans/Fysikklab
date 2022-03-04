# TFY41xx Fysikk vaaren 2021.
#
# Programmet tar utgangspunkt i hoeyden til de 8 festepunktene.
# Deretter beregnes baneformen y(x) ved hjelp av 7 tredjegradspolynomer, 
# et for hvert intervall mellom to festepunkter, slik at baade banen y, 
# dens stigningstall y' = dy/dx og dens andrederiverte
# y'' = d2y/dx2 er kontinuerlige i de 6 indre festepunktene.
# I tillegg velges null krumning (andrederivert) 
# i banens to ytterste festepunkter (med bc_type='natural' nedenfor).
# Dette gir i alt 28 ligninger som fastlegger de 28 koeffisientene
# i de i alt 7 tredjegradspolynomene.

# De ulike banene er satt opp med tanke paa at kula skal 
# (1) fullfoere hele banen selv om den taper noe mekanisk energi underveis;
# (2) rulle rent, uten aa gli ("slure").

# Vi importerer noedvendige biblioteker:
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import CubicSpline
import math

# Horisontal avstand mellom festepunktene er 0.200 m
h = 0.200
xfast=np.asarray([0,h,2*h,3*h,4*h,5*h,6*h,7*h])


# Vi begrenser starthøyden (og samtidig den maksimale høyden) til
# å ligge mellom 250 og 300 mm
ymax = 300
# yfast: tabell med 8 heltall mellom 50 og 300 (mm); representerer
# høyden i de 8 festepunktene
yfast=np.asarray(np.random.randint(50, ymax, size=8))
#konverter fra m til mm
yfast =yfast/1000
# inttan: tabell med 7 verdier for (yfast[n+1]-yfast[n])/h (n=0..7); dvs
# banens stigningstall beregnet med utgangspunkt i de 8 festepunktene.
inttan = np.diff(yfast)/h
attempts=1
# while-løkken sjekker om en eller flere av de 3 betingelsene ovenfor
# ikke er tilfredsstilt; i så fall velges nye festepunkter inntil
# de 3 betingelsene er oppfylt
# while (yfast[0] < yfast[1]*1.04 or
#        yfast[0] < yfast[2]*1.08 or
#        yfast[0] < yfast[3]*1.12 or
#        yfast[0] < yfast[4]*1.16 or
#        yfast[0] < yfast[5]*1.20 or
#        yfast[0] < yfast[6]*1.24 or
#        yfast[0] < yfast[7]*1.28 or
#        yfast[0] < 0.250 or
#        np.max(np.abs(inttan)) > 0.4 or
#        inttan[0] > -0.2):
#           yfast=np.asarray(np.random.randint(0, ymax, size=8))
#
#           #konverter fra m til mm
#           yfast =yfast/1000
#
#           inttan = np.diff(yfast)/h
#           attempts=attempts+1

yfast = np.asarray([0.350,0.315,0.299,0.306,0.281,0.236,0.168,0.200])

# Omregning fra mm til m:
# xfast = xfast/1000
# yfast = yfast/1000

# Når programmet her har avsluttet while-løkka, betyr det at
# tallverdiene i tabellen yfast vil resultere i en tilfredsstillende bane. 

#Programmet beregner deretter de 7 tredjegradspolynomene, et
#for hvert intervall mellom to nabofestepunkter.


#Med scipy.interpolate-funksjonen CubicSpline:
cs = CubicSpline(xfast, yfast, bc_type='natural')

xmin = 0.000
xmax = 1.401
dx = 0.001

x = np.arange(xmin, xmax, dx)   

#funksjonen arange returnerer verdier paa det "halvaapne" intervallet
#[xmin,xmax), dvs slik at xmin er med mens xmax ikke er med. Her blir
#dermed x[0]=xmin=0.000, x[1]=xmin+1*dx=0.001, ..., x[1400]=xmax-dx=1.400, 
#dvs x blir en tabell med 1401 elementer
Nx = len(x)
y = cs(x)       #y=tabell med 1401 verdier for y(x)
dy = cs(x,1)    #dy=tabell med 1401 verdier for y'(x)
d2y = cs(x,2)   #d2y=tabell med 1401 verdier for y''(x)

#Eksempel: Plotter banens form y(x)
baneform = plt.figure('y(x)',figsize=(12,6))
plt.plot(x,y,xfast,yfast,'*')
plt.title('Banens form')
plt.xlabel('$x$ (m)',fontsize=20)
plt.ylabel('$y(x)$ (m)',fontsize=20)
plt.ylim(0.0,0.40)
plt.grid()
plt.show()
#Figurer kan lagres i det formatet du foretrekker:
#baneform.savefig("baneform.pdf", bbox_inches='tight')
baneform.savefig("baneform.png", bbox_inches='tight')
#baneform.savefig("baneform.eps", bbox_inches='tight')

#konstanter
c = 2/5
g = 9.81

#v(x)
velocity = np.sqrt((2*g*(y[0]- y))/(1+c))

baneform = plt.figure('v(x)',figsize=(12,6))
plt.plot(x,velocity)
plt.title('Banens form')
plt.xlabel('$x$ (m)',fontsize=20)
plt.ylabel('$v(x)$ (m/s)',fontsize=20)
plt.grid()
plt.show()

#Beta (helningsvinkel)
beta = np.arctan(dy)

baneform = plt.figure('Beta',figsize=(12,6))
plt.plot(x,beta)
plt.title('Banens form')
plt.xlabel('$x$ (m)',fontsize=20)
plt.ylabel('$Beta$',fontsize=20)
plt.grid()
plt.show()

#x(t)
dt = np.zeros(Nx)
t = np.zeros(Nx)

vx = velocity*np.cos(beta)

for i in range(1, Nx):
    dt[i] = (2*dx)/(vx[i-1]+vx[i])
    t[i] = t[i-1] + dt[i]

baneform = plt.figure('X(t)',figsize=(12,6))
plt.plot(t,x)
plt.title('Banens form')
plt.xlabel('$t$ (s)',fontsize=20)
plt.ylabel('$x$ (m)',fontsize=20)
plt.grid()
plt.show()

#f/N
M = 0.031 #kg
a = ((velocity**2)*d2y)/np.sqrt((1+(dy**2))**3)
N = M*(g*np.cos(beta) + a)
f = (c*M*g*np.sin(beta))/(1+c)
fn = np.abs(f/N)

baneform = plt.figure('f/N',figsize=(12,6))
plt.plot(x,fn)
plt.title('Banens form')
plt.xlabel('$x$ (m)',fontsize=20)
plt.ylabel('$|f/N|$',fontsize=20)
plt.grid()
plt.show()



print('Antall forsøk',attempts)
print('Festepunkthøyder (m)',yfast)
print('Banens høyeste punkt (m)',np.max(y))

print('NB: SKRIV NED festepunkthøydene når du/dere er fornøyd med banen.')
print('Eller kjør programmet på nytt inntil en attraktiv baneform vises.')