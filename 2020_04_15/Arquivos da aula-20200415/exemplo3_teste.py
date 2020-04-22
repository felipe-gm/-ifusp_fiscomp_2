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

def f(theta,t):
    return -g/l*sin(theta)

def passo_vrl(f,r,v,vpm,t,h):   # Calcula um passo no método de Verlet
    r += h*vpm                  # Calculando a posição no próximo t "inteiro"
    k = h*f(r,t+h)
    v = vpm + 0.5*k             # Calculando a velocidade no próximo t "inteiro"
    vpm += k                    # Calculando a velocidade no próximo t "médio"
    return r, v, vpm            # Retorna os NOVOS VALORES de r, v e vpm

a = 0.0             # Início do intervalo
b = 10.0            # Final do intervalo

theta_a, omega_a = 2*pi/3, 0.0                     # Condições iniciais
E_a = m*(l*omega_a)**2/2 + m*g*l*(1-cos(theta_a))  # Energia inicial 


h = 5e-2       # Tamanho do passo de integração
theta_vrl, omega_vrl = [], []
E_vrl = []
t_lista = []
r_vrl, v_vrl = theta_a, omega_a
t = 0
k = (h/2)*f(r_vrl,t) # Inicializando o ponto "médio" com RK2
vpm = v_vrl + (h/2)*f(r_vrl+k/2,t+(h/2)/2)
while t<= b:   # Realizando a integração numérica
    t_lista.append(t)
    theta_vrl.append(r_vrl)
    omega_vrl.append(v_vrl)
    E_vrl.append(m*(l*v_vrl)**2/2 + m*g*l*(1-cos(r_vrl)) - E_a) # Var. energia
    r_vrl, v_vrl, vpm = passo_vrl(f,r_vrl,v_vrl,vpm,t,h)
    t += h
    
plt.figure(figsize=(12,9))
plt.plot(t_lista,E_vrl)
plt.title("Verlet")
plt.xlabel("t")
plt.ylabel("E(t) - E(0)")
plt.show()

plt.figure(figsize=(12,9))
plt.plot(theta_vrl,omega_vrl)
plt.title("Verlet")
plt.xlabel("theta(t)")
plt.ylabel("omega(t)")
plt.show()
