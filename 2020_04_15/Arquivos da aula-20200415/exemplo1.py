from math import sin,cos,pi
from numpy import arange, array
import matplotlib.pyplot as plt

# Parâmetros da exibição dos gráficos
plt.rcParams['xtick.labelsize'] = 22
plt.rcParams['ytick.labelsize'] = 22
plt.rcParams['axes.labelsize'] = 26
plt.rcParams['axes.titlesize'] = 30

# Constantes
m, g, l = 1.0, 9.8, 1.0

def f(r,t):
    theta, omega = r[0], r[1]
    return array([omega,-g/l*sin(theta)],float)

def passo_rk2(f,r,t,h):            # Calcula um passo no método de RK2
    k1 = h*f(r,t)
    k2 = h*f(r+0.5*k1,t+0.5*h)
    return k2

def passo_rk4(f,r,t,h):            # Calcula um passo no método de RK4
    k1 = h*f(r,t)
    k2 = h*f(r+0.5*k1,t+0.5*h)
    k3 = h*f(r+0.5*k2,t+0.5*h)
    k4 = h*f(r+k3,t+h)
    return (k1+2.0*(k2+k3)+k4)/6.0

a = 0.0             # Início do intervalo
b = 10.0            # Final do intervalo

theta_a, omega_a = 2*pi/3, 0.0                      # Condições iniciais
E_a = m*(l*omega_a)**2/2 + m*g*l*(1-cos(theta_a)) # Energia inicial 


h = 3e-2            # Tamanho do passo de integração
theta_rk2, theta_rk4 = [], []
E_rk2, E_rk4 = [], []
t_lista = []
r_rk2 = array([theta_a,omega_a],float)
r_rk4 = array([theta_a,omega_a],float)
t = 0
while t<= b:   # Realizando as integrações numéricas
    t_lista.append(t)
    theta_rk2.append(r_rk2[0])
    theta_rk4.append(r_rk4[0])
    E_rk2.append(m*(l*r_rk2[1])**2/2 + m*g*l*(1-cos(r_rk2[0])) - E_a)
    E_rk4.append(m*(l*r_rk4[1])**2/2 + m*g*l*(1-cos(r_rk4[0])) - E_a)
    r_rk2 += passo_rk2(f,r_rk2,t,h)
    r_rk4 += passo_rk4(f,r_rk4,t,h)
    t += h
    
plt.figure(figsize=(12,9))
plt.plot(t_lista,E_rk2,t_lista,E_rk4)
plt.title("Runge-Kutta (2 e 4)")
plt.xlabel("t")
plt.ylabel("E(t) - E(0)")
plt.show()

plt.figure(figsize=(14,9))
plt.plot(t_lista,E_rk2,t_lista,E_rk4)
plt.ylim(-0.0001,0.0001)
plt.title("Runge-Kutta (2 e 4)")
plt.xlabel("t")
plt.ylabel("E(t) - E(0)",labelpad=-950)
plt.show()
