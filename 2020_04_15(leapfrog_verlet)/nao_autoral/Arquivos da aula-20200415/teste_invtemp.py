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

theta_0, omega_0 = 1.0, 0.0                      # Condições iniciais

h = 1e-1            # Tamanho do passo de integração
print("Tamanho do passo de integração: h = ",h) 
print()

print("Integrando com RK2")
r_0 = array([theta_0,omega_0],float)
print("  Partimos de theta(0)=",theta_0,"com omega(0)=",omega_0)

t = 0
r_1 = r_0 + passo_rk2(f,r_0,t,h)
print("  Em direção ao futuro obtemos theta(h)=",r_1[0],"e omega(h)=",r_1[1])
print()

t = h
print("  Partimos de theta(h)=",r_1[0],"com omega(h)=",r_1[1])
r_zero = r_1 + passo_rk2(f,r_1,t,-h)
print("  Em direção ao passado obtemos theta(0)=",r_zero[0],"e omega(0)=",r_zero[1])

print()
print("Integrando com RK4")
r_0 = array([theta_0,omega_0],float)
print("  Partimos de theta(0)=",theta_0,"com omega(0)=",omega_0)

t = 0
r_1 = r_0 + passo_rk4(f,r_0,t,h)
print("  Em direção ao futuro obtemos theta(h)=",r_1[0],"e omega(h)=",r_1[1])
print()

t = h
print("  Partimos de theta(h)=",r_1[0],"com omega(h)=",r_1[1])
r_zero = r_1 + passo_rk4(f,r_1,t,-h)
print("  Em direção ao passado obtemos theta(0)=",r_zero[0],"e omega(0)=",r_zero[1])

