from math import pi,exp
from random import random,randrange
from numpy import zeros,ones
import matplotlib.pyplot as plt

# Este programa implementa a simulação de Monte Carlo de um gás de bósons
# não interagentes. Os valores permitidos para cada número quântico nx, ny
# ou nz são números inteiros quaisquer, o que corresponde a uma escolha
# de condições de contorno periódicas para o cálculo das funções de onda.
# Caso escolhêssemos condições de contorno fixas, os valores permitidos
# seriam apenas os inteiros positivos. (É essa a escolha feita por Newman.)
# No limite em que N e L tendem ao infinito, não há diferença nas propriedades
# físicas previstas por ambas as escolhas.

# Parâmetros e constantes da simulação
N = 10**3                   # Número de átomos
L = 1e-6                    # Lado do cubo
h = 1.055e-34               # Valor da constante de Planck dividida por 2*pi
m = 6.645e-27               # Massa de um átomo de hélio-4
epsilon = (pi*h/L)**2/(2*m) # Escala de energia
kB = 1.381e-23              # Constante de Boltzmann
T = 10*epsilon/kB            # Temperatura
Passos = 300                # Número de passos de Monte Carlo da simulação
kBT = kB*T

# Matriz que armazena os números quânticos
n = zeros([N,3],int)        # Partimos do estado fundamental

# Laço principal
E_lista = []                # Lista para armazenar a energia
E = 0.0
for i in range(N):
    E += n[i,0]**2 + n[i,1]**2 + n[i,2]**2
E *= epsilon
for k in range(Passos):     # Percorremos os passos de Monte Carlo
    # Escolhemos um número quântico que tentaremos alterar
    for l in range(N):          # Cada partícula em média tem 1 chance de mudar 
        i = randrange(N)               # Partícula
        j = randrange(3)               # Número quântico da partícula
        if random() < 0.5:             # Com probabilidade 1/2 propomos n->n+1
            dn = 1
            dE = (2*n[i,j]+1)*epsilon  # Variação da energia
        else:                          # Caso contrário propomos n->n-1
            dn = -1             
            dE = (-2*n[i,j]+1)*epsilon
        if random() < exp(-dE/kBT):    # Decidimos se aceitamos a alteração
            n[i,j] += dn
            E += dE
    # Registramos a energia a cada passo de Monte Carlo
    E_lista.append(E/N/kBT)        # Energia por partícula (medida em kB*T)
    
# Parâmetros da exibição dos gráficos
plt.rcParams['xtick.labelsize'] = 28
plt.rcParams['ytick.labelsize'] = 28
plt.rcParams['axes.labelsize'] = 32

plt.figure(figsize=(12,9))
plt.plot(E_lista)
plt.ylabel("Energia por partícula (em $k_B T$)")
plt.xlabel("Passo de Monte Carlo")
plt.ylim(0.05,1.7)
plt.show()

