from random import random
from numpy import arange
import matplotlib.pyplot as plt

# Constantes
NTl = 1000          # Número inicial de átomos de tálio
NPb = 0             # Número inicial de átomos de chumbo
tau = 3.053*60      # Meia-vida do tálio em segundos
h = 1.0             # Tamanho do passo de tempo em segundos
p = 1 - 2**(-h/tau) # Probabilidade de decaimento em um passo
tmax = 1000         # Tempo total da simulação

# Criando as listas para os gráficos
t_lista = arange(0.0,tmax,h)
Tl_lista, Pb_lista = [], []

# Laço principal
for t in t_lista:
    Tl_lista.append(NTl)
    Pb_lista.append(NPb)

    # Dando a cada átomo a chance de decair
    decaimento = 0
    for i in range(NTl):
        if random() < p:
            decaimento += 1
    NTl -= decaimento
    NPb += decaimento

# Parâmetros da exibição dos gráficos
plt.rcParams['xtick.labelsize'] = 28
plt.rcParams['ytick.labelsize'] = 28
plt.rcParams['axes.labelsize'] = 32

# Traçando o gráfico
plt.figure(figsize=(12,9))
plt.semilogy(t_lista,Tl_lista,"b.")
plt.semilogy(t_lista,Pb_lista,"r.")
plt.xlabel("Tempo")
plt.ylabel("Número de átomos")
plt.show()
            

