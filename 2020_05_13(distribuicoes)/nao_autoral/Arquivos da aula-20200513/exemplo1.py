from random import random
from math import log
from numpy import sort,arange
import matplotlib.pyplot as plt

# Constantes
NTl = 10**4         # Número inicial de átomos de tálio
tau = 3.053*60      # Meia-vida do tálio em segundos

mu = log(2)/tau     # Constante da distribuição exponencial

# Sorteando os tempos de decaimento de cada átomo
t_lista = []
for i in range(NTl):
    t_lista.append(-log(1-random())/mu)

# Vamos fazer um gráfico, em função de t, do número de átomos que não
# decaíram até um certo instante t. Para isso, é útil reordenar a lista
# de tempos do menor para o maior.
t_ordenado = sort(t_lista)
N_ordenado = arange(NTl,0,-1)

# Parâmetros da exibição dos gráficos
plt.rcParams['xtick.labelsize'] = 28
plt.rcParams['ytick.labelsize'] = 28
plt.rcParams['axes.labelsize'] = 32

# Traçando o gráfico
plt.figure(figsize=(12,9))
plt.xlabel("$t$")
plt.ylabel("$N(t)$")
plt.semilogy(t_ordenado,N_ordenado,"b.")
plt.show()
            

