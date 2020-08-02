from math import sin,cos,exp,sqrt,log
from numpy import arange
import matplotlib.pyplot as plt

# Parâmetros da exibição dos gráficos
plt.rcParams['xtick.labelsize'] = 28
plt.rcParams['ytick.labelsize'] = 28
plt.rcParams['axes.labelsize'] = 32

c = 1.0

def f(x,t):
    return -2.0*exp(-c*t)*(c*cos(t)+sin(t))

def sol_exata(t,a,xa):
    b = xa - 2.0*exp(-c*a)*cos(a)
    return 2.0*exp(-c*t)*cos(t) + b

a = 0.0           # Início do intervalo
b = 10.0          # Final do intervalo
xa = 0.0          # Condição inicial, ou seja, x(a)

N_exato = 1000    # Número de pontos para a sol. exata
h_exato = (b-a)/N_exato

x_exato = [] 
t_exato = arange(a,b,h_exato)
for t in t_exato:
    x_exato.append(sol_exata(t,a,xa))

N_lista = [8,16,32,64,128,256,512,1024]
erromedio_rk2 = []
erromedio_euler = []

for N in N_lista:     # Número de passos da sol. por Runge-Kutta
    h = (b-a)/N       # Tamanho de um passo dessa solução
    x = xa
    
    t_rk2 = arange(a,b,h)
    x_rk2 = []
    x_euler = []
    xe = x
    erro_rk2 = 0.0
    erro_euler = 0.0

    for t in t_rk2:
        
        erro_rk2 += (x - sol_exata(t,a,xa))**2     # Acumulando o erro no cálculo de Runge-Kutta
        erro_euler += (xe - sol_exata(t,a,xa))**2  # Acumulando o erro no cálculo de Euler

        x_rk2.append(x)
        k1 = h*f(x,t)
        k2 = h*f(x+0.5*k1,t+0.5*h)
        x += k2
        x_euler.append(xe)
        xe += k1
        
    if N == 128:
       plt.figure(figsize=(16,12))
       plt.plot(t_rk2,x_rk2,'gs',t_rk2,x_euler,'ro',t_exato,x_exato)
       plt.xlabel("t")
       plt.ylabel("x(t)")
       plt.show() 

    erro_rk2 /= len(t_rk2)
    erro_euler /= len(t_rk2)
    erromedio_rk2.append(sqrt(erro_rk2))
    erromedio_euler.append(sqrt(erro_euler))

plt.figure(figsize=(16,12))
plt.loglog(N_lista,erromedio_rk2,'gs',N_lista,erromedio_euler,'ro')
plt.xlabel("N")
plt.ylabel("Erro médio")
plt.show()

print((log(erromedio_euler[0])-log(erromedio_euler[len(N_lista)-1]))/(log(1024.0)-log(8.0)))
print((log(erromedio_rk2[0])-log(erromedio_rk2[len(N_lista)-1]))/(log(1024.0)-log(8.0)))
