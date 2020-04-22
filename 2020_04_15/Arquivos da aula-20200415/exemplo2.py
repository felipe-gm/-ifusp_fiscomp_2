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
    return k2                      # Retorna o INCREMENTO em r

def passo_lpf(f,r,rpm,t,h):        # Calcula um passo no método leapfrog
    r += h*f(rpm,t+0.5*h)          # Calculando valores nos pontos "inteiros"
    rpm += h*f(r,t+h)              # Calculando valores nos pontos "médios"
    return r, rpm                  # Retorna os NOVOS VALORES de r e rpm

a = 0.0             # Início do intervalo
b = 10.0            # Final do intervalo

theta_a, omega_a = 2*pi/3, 0.0                    # Condições iniciais
E_a = m*(l*omega_a)**2/2 + m*g*l*(1-cos(theta_a)) # Energia inicial 


h = 5e-2       # Tamanho do passo de integração
theta_rk2, theta_lpf = [], []
E_rk2, E_lpf = [], []
t_lista = []
r_rk2 = array([theta_a,omega_a],float)
r_lpf = array([theta_a,omega_a],float)
t = 0
rpm = r_lpf + passo_rk2(f,r_lpf,t,h/2)  # Inicializando valores no ponto "médio"
while t<= b:   # Realizando as integrações numéricas
    t_lista.append(t)
    theta_rk2.append(r_rk2[0])
    theta_lpf.append(r_lpf[0])
    E_rk2.append(m*(l*r_rk2[1])**2/2 + m*g*l*(1-cos(r_rk2[0])) - E_a)
    E_lpf.append(m*(l*r_lpf[1])**2/2 + m*g*l*(1-cos(r_lpf[0])) - E_a)
    r_rk2 += passo_rk2(f,r_rk2,t,h)
    r_lpf, rpm = passo_lpf(f,r_lpf,rpm,t,h)
    t += h
    
plt.figure(figsize=(12,9))
plt.plot(t_lista,E_rk2,t_lista,E_lpf)
plt.title("RK2 e leapfrog")
plt.xlabel("t")
plt.ylabel("E(t) - E(0)")
plt.show()
