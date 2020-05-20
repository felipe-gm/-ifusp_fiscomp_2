from random import random
from math import pi,exp,sqrt,cos,sin,log
from numpy import arange,zeros
import matplotlib.pyplot as plt

# Constantes
N = 10**6               # Número inicial de píons
L = 20.0                # Distância que os píons percorrem
tau = 2.6e-8            # Meia-vida do píon no repouso, em segundos
c = 3e8                 # Velocidade da luz, em metros por segundo
mc2 = 140               # Energia de repouso de um píon, em MeV
K_medio = 200           # Energia cinética média dos píons, em MeV
dK = 30                 # Desvio-padrão da energia cinética dos píons

# Função que gera pares de números aleatórios gaussianos
mu = K_medio
sigma = dK
doispi = 2.0*pi
doissigma2 = 2.0*sigma*sigma
def gaussian():
    r = sqrt(-doissigma2*log(1-random()))
    theta = doispi*random()
    return mu+r*cos(theta), mu+r*sin(theta)

# Vamos sortear as energias de cada píon e delas calcular a velocidade
# (em unidades de c).
v_lista = zeros(N) 
for n in range(N//2):
    K0, K1 = gaussian()
    v_lista[2*n] = sqrt(1-1/(1+K0/mc2)**2)*c
    v_lista[2*n+1] = sqrt(1-1/(1+K1/mc2)**2)*c

# Para cada píon, determinamos o tempo próprio de percurso delta_t
# e a probabilidade p de que o píon sobreviva ao percurso.
sobreviventes = 0
for n in range(N):
    delta_t = sqrt(1-(v_lista[n]/c)**2)*L/v_lista[n]
    p = 2**(-delta_t/tau)
    if random() < p:
        sobreviventes +=1

print("Fração de píons sobreviventes =",sobreviventes/N)

