from random import random
from math import sqrt,pi,exp
from numpy import arange
from scipy.special import erfinv
import matplotlib.pyplot as plt
import time

tempo_inicial = time.time()

# Constantes
mu = 2.0      # Valor médio desejado das variáveis aleatórias gaussianas
sigma2 = 3.0  # Valor desejado da variância da gaussiana
N = 10**5     # Número de variáveis aleatórias gaussianas a obter
Bins = 10**2  # Número de caixas do histograma das variáveis gaussianas

# O laço abaixo produz N números aleatórios gaussianos utilizando a
# função erro inversa do pacote 'scipy'.
X_lista = []
for i in range(N):
    X_lista.append(mu + sqrt(2*sigma2)*erfinv(2*random()-1))

print("Tempo de execução (em segundos):",(time.time() - tempo_inicial))

# Vamos também produzir o gráfico de uma distribuição gaussiana com
# média mu e variância sigma, para comparação.
XX_min, XX_max = min(X_lista), max(X_lista)
h = (XX_max - XX_min)/Bins
XX_lista, Gauss = arange(XX_min,XX_max,h), []
Gauss[:] = [1/sqrt(2*pi*sigma2)*exp(-(X-mu)**2/2/sigma2) for X in XX_lista]

# Parâmetros da exibição dos gráficos
plt.rcParams['xtick.labelsize'] = 28
plt.rcParams['ytick.labelsize'] = 28
plt.rcParams['axes.labelsize'] = 32

# Traçando os gráficos
plt.figure(figsize=(12,9))
plt.hist(X_lista,bins=Bins,density=True) # Produzindo um histograma normalizado
plt.plot(XX_lista,Gauss)
plt.xlabel("$X$")
plt.ylabel("$P(X)$")
plt.show()
            

