from math import sin,cos,pi,sqrt
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

def passo_eul(f,x,t,h):            # Calcula um passo no método de Euler
    k1 = h*f(x,t)
    return k1                      # Retorna o INCREMENTO em r

def passo_rk2(f,r,t,h):            # Calcula um passo no método de RK2
    k1 = h*f(r,t)
    k2 = h*f(r+0.5*k1,t+0.5*h)
    return k2                      # Retorna o INCREMENTO em r

def passo_rk4(f,r,t,h):            # Calcula um passo no método de RK4
    k1 = h*f(r,t)
    k2 = h*f(r+0.5*k1,t+0.5*h)
    k3 = h*f(r+0.5*k2,t+0.5*h)
    k4 = h*f(r+k3,t+h)
    return (k1+2.0*(k2+k3)+k4)/6.0 # Retorna o INCREMENTO em r

def passo_lpf(f,r,rpm,t,h):        # Calcula um passo no método leapfrog
    r += h*f(rpm,t+0.5*h)          # Calculando valores nos pontos "inteiros"
    rpm += h*f(r,t+h)              # Calculando valores nos pontos "médios"
    return r, rpm                  # Retorna os NOVOS VALORES de r e rpm

a = 0.0             # Início do intervalo
b = 10.0            # Final do intervalo

theta_a, omega_a = 2*pi/3, 0.0                       # Partindo do repouso
#theta_a, omega_a = 0.0, sqrt(2*g/l*(1-cos(theta_a))) # Partindo do fundo         
E_a = m*(l*omega_a)**2/2 + m*g*l*(1-cos(theta_a)) # Energia inicial 


h = 5e-2       # Tamanho do passo de integração
theta_rk2, theta_eul = [], []
omega_rk2, omega_eul = [], []
theta_rpm_rk2, theta_rpm_eul = [], []
omega_rpm_rk2, omega_rpm_eul = [], []
E_rk2, E_eul = [], []
t_lista = []
r_rk2 = array([theta_a,omega_a],float)
r_eul = array([theta_a,omega_a],float)
t = 0
rpm_eul = r_eul + passo_eul(f,r_eul,t,h/2) # Inicializando valores no ponto "médio"
rpm_rk2 = r_rk2 + passo_rk2(f,r_rk2,t,h/2) # Inicializando valores no ponto "médio"
while abs(t)<= b:   # Realizando as integrações numéricas
    t_lista.append(t)
    theta_rk2.append(r_rk2[0])
    theta_eul.append(r_eul[0])
    omega_rk2.append(r_rk2[1])
    omega_eul.append(r_eul[1])
    theta_rpm_rk2.append(rpm_rk2[0])
    theta_rpm_eul.append(rpm_eul[0])
    omega_rpm_rk2.append(rpm_rk2[1])
    omega_rpm_eul.append(rpm_eul[1])
    E_rk2.append(m*(l*r_rk2[1])**2/2 + m*g*l*(1-cos(r_rk2[0])) - E_a)
    E_eul.append(m*(l*r_eul[1])**2/2 + m*g*l*(1-cos(r_eul[0])) - E_a)
    r_rk2, rpm_rk2 = passo_lpf(f,r_rk2,rpm_rk2,t,h)
    r_eul, rpm_eul = passo_lpf(f,r_eul,rpm_eul,t,h)
    t += h
    
plt.figure(figsize=(12,9))
plt.plot(t_lista,theta_eul,t_lista,theta_rk2)
plt.title("Leapfrog com passo inicial de Euler ou RK2")
plt.xlabel("t")
plt.ylabel("theta(t)")
plt.show()

plt.figure(figsize=(12,9))
plt.plot(t_lista,theta_rpm_eul,t_lista,theta_rpm_rk2)
plt.title("Leapfrog com passo inicial de Euler ou RK2")
plt.xlabel("t")
plt.ylabel("theta_rpm(t)")
plt.show()

plt.figure(figsize=(12,9))
plt.plot(t_lista,omega_eul,t_lista,omega_rk2)
plt.title("Leapfrog com passo inicial de Euler ou RK2")
plt.xlabel("t")
plt.ylabel("omega(t)")
plt.show()

plt.figure(figsize=(12,9))
plt.plot(t_lista,omega_rpm_eul,t_lista,omega_rpm_rk2)
plt.title("Leapfrog com passo inicial de Euler ou RK2")
plt.xlabel("t")
plt.ylabel("omega_rpm(t)")
plt.show()

plt.figure(figsize=(12,9))
plt.plot(t_lista,E_eul,t_lista,E_rk2)
plt.title("Leapfrog com passo inicial de Euler ou RK2")
plt.xlabel("t")
plt.ylabel("E(t) - E(0)")
plt.show()
