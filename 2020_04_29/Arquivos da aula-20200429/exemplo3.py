# Problema de valores de contorno: que velocidade inicial faz com que uma
# bola lançada para cima no instante 'ta' retorne à E inicial no
# instante 'tb'?
from math import sin,cos,exp,sqrt,log
from numpy import arange,array
import matplotlib.pyplot as plt

# Parâmetros da exibição dos gráficos
plt.rcParams['xtick.labelsize'] = 24
plt.rcParams['ytick.labelsize'] = 24
plt.rcParams['axes.labelsize'] = 28
plt.rcParams['text.usetex'] = True

# Construímos uma equação semelhante à de Schrödinger, mas com unidades
# arbitrárias

xa = -4.          # Início do intervalo da variável independente
xb = 4.           # Final do intervalo da variável independente
psia = 0.0        # Valor da variável independente no início do intervalo
psib = 0.0        # Valor alvo da variável independente no final do intervalo
h = 0.01          # Tamanho do passo de integração
prec = 1e-14      # Precisão do resultado para a energia

beta=0.25
def V(x):
    return beta*x**4

def f(r,x,E):
    psi, dpsi = r[0], r[1]
    f0, f1 = dpsi, 2*(V(x)-E)*psi
    return array([f0,f1],float)

def passo_rk4(f,r,x,h,E):            # Calcula um passo no método de RK4
    k1 = h*f(r,x,E)
    k2 = h*f(r+0.5*k1,x+0.5*h,E)
    k3 = h*f(r+0.5*k2,x+0.5*h,E)
    k4 = h*f(r+k3,x+h,E)
    return (k1+2.0*(k2+k3)+k4)/6.0

# Função que calcula a função de onda associada a uma energia E.
# O resultado é normalizado para que a integral do quadrado da
# função de onda seja igual a 1. Há também uma variável de
# reescala ('esc') que permite acomodar várias curvas em um só gráfico.
def produz_lista_psi(x_lista,ra,E,esc):
    r = ra
    psi_lista = []
    norm = 0.0
    for x in x_lista:
        psi_lista.append(r[0])
        norm += h*r[0]**2 # Calculando a integral do quadrado da função de onda
        r += passo_rk4(f,r,x,h,E)
    norm = sqrt(norm)
    psi_lista[:] = [esc*psi/norm + E for psi in psi_lista] 
    return psi_lista

# Função que integra a equação diferencial definida por f entre 'xa' e 'xb'
# com condição inicial 'ra' e retorna o valor final do vetor 'r'
def r_final(f,xa,xb,ra,h,E):
    r = ra.copy() # Não utilizar 'copy' faz com que mudar 'r' afete 'ra'
    for x in arange(xa,xb,h):
        r += passo_rk4(f,r,x,h,E)
    return r

# Função que determina uma energia permitida, impondo que a função de onda
# seja nula nos extremos do intervalo de x
def determina_E(E1,E2,ra,prec): # Solução via método da secante
    psi2 = r_final(f,xa,xb,ra,h,E1)[0] # Valor de psib com energia E1
    while abs((E2-E1)/((E1+E2)/2)) > prec:
        psi1 = psi2
        psi2 = r_final(f,xa,xb,ra,h,E2)[0] # Valor de psib com energia 'E2' 
        E1, E2 = E2, E2 - psi2*(E2-E1)/(psi2-psi1)
    return (E1+E2)/2       # Resultado final do cálculo

dpsi = 1e-6
ra = array([psia,dpsi],float)
E_lista = []

#
# Determinando as energias 
#
nmax = 8 # Número de energias a determinar
for n in range(nmax):
    E1 = 1.5*(n+0.02*n**2)  # Dois palpites para iniciar a solução
    E2 = E1 + 0.1
    E_lista.append(determina_E(E1,E2,ra,prec))
    print("A energia do nível",n,"é",E_lista[n])

# Produzindo gráficos da energia e comparando com a do oscilador harmônico
plt.figure(figsize=(12,9))
n_lista = array(range(1,nmax))
DeltaEohs, DeltaE = [], []
DeltaEohs[:] = [1 for i in range(1,nmax)]
DeltaE[:] = [E_lista[i]-E_lista[i-1] for i in range(1,nmax)]
plt.xlabel("$n$")
plt.ylabel("$E_n - E_{n-1}$")
plt.plot(n_lista,DeltaEohs,'b-',n_lista,DeltaE,'r-')
plt.plot(n_lista,DeltaEohs,'bo',n_lista,DeltaE,'rs')
plt.show()

# Produzindo gráficos das funções de onda e do potencial
esc=0.8
plt.figure(figsize=(12,9))
x_lista = arange(xa,xb,h)
V_lista = []
for x in x_lista:
    V_lista.append(V(x))
plt.plot(x_lista,V_lista,'k--')
plt.xlabel("$x$")
plt.ylabel("$\psi(x)$")
plt.axvline(0,0,E_lista[-1]+1,color='k',linewidth=0.5)
plt.xlim(xa,xb)
plt.ylim(0.0,1+E_lista[-1])
for E in E_lista:
    plt.fill_between(x_lista,E,produz_lista_psi(x_lista,ra,E,esc),alpha=0.5)
plt.show()
