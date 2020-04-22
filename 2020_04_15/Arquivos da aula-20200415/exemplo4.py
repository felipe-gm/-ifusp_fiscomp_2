# Órbita da Terra em torno do Sol (Newman exercício 8.12)
from math import sin,cos,exp,sqrt,log
from numpy import arange,array,linalg
import matplotlib.pyplot as plt
# Todas as grandezas estão em unidades do SI

# Constantes (constante gravitacional, massa do Sol, massa da Terra
G, M, m = 6.67380e-11, 1.9891e30, 5.9722e24 
GM = G*M
GMm = GM*m

# Limites e passo de integração
a = 0.                  # Instante inicial em segundos
b = 16*365*(60*60*24)   # Instante final: 16 anos (em segundos)
h = 30*24*60*60         # Passo de tempo: 30 dias (em segundos)

# Condições inicias (supondo o Sol na origem)
r_a = array([-1.4710e11,0.,0.]) # Posição da Terra no periélio
v_a = array([0.,-3.0287e4,0.0]) # Velocidade da Terra no periélio

def f(r,t):                     # Força/massa atuando sobre a esfera
    magr = linalg.norm(r)
    return -GM*r/magr**3

def passo_vrl(f,r,v,vpm,t,h):   # Calcula um passo no método de Verlet
    r += h*vpm                  # Calculando a posição no próximo t "inteiro"
    k = h*f(r,t+h)
    v = vpm + 0.5*k             # Calculando a velocidade no próximo t "inteiro"
    vpm += k                    # Calculando a velocidade no próximo t "médio"
    return r, v, vpm            # Retorna os NOVOS VALORES de r, v e vpm

t = 0
r_vrl, v_vrl = r_a, v_a
vpm = v_vrl + 0.5*h*f(r_vrl,t)  # Inicializando valores no ponto "médio"
t_lista, K_lista, Ug_lista, E_lista = [], [], [], []
x_lista, y_lista = [], []
while t<= b:   # Realizando a integração numérica
    K = 0.5*m*linalg.norm(v_vrl)**2
    U_g = -GMm/linalg.norm(r_vrl)
    t_lista.append(t)
    K_lista.append(K)
    Ug_lista.append(U_g)
    E_lista.append(K+U_g)
    x_lista.append(r_vrl[0])
    y_lista.append(r_vrl[1])
    r_vrl, v_vrl, vpm = passo_vrl(f,r_vrl,v_vrl,vpm,t,h)
    t += h

# Parâmetros da exibição dos gráficos
plt.rcParams['xtick.labelsize'] = 20
plt.rcParams['ytick.labelsize'] = 20
plt.rcParams['axes.labelsize'] = 24
plt.rcParams['axes.titlesize'] = 28

plt.figure(figsize=(9,9))
plt.plot(x_lista,y_lista)
plt.title("Verlet")
plt.xlabel("x")
plt.ylabel("y")
plt.show()

plt.figure(figsize=(12,9))
plt.plot(t_lista,K_lista)
plt.plot(t_lista,Ug_lista)
plt.plot(t_lista,E_lista)
plt.title("Verlet")
plt.xlabel("t")
plt.ylabel("Energia")
plt.show()
