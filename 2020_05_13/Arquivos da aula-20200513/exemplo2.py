from random import random
from math import sqrt,pi,exp
from numpy import arange
import matplotlib.pyplot as plt
import time

tempo_inicial = time.time()

# O teorema central do limite garante que, sob condições bastante amplas,
# a soma S de variáveis aleatórias independentes z extraídas de uma
# distribuição de probabilidades com média mu e variância sigma2 é descrita
# por uma distribuição de probabilidades gaussiana.
# Concretamente, se S_n é a média de n variáveis z, a distribuição da
# quantidade X = mu + sqrt(n)*(S_n - mu) converge quando n vai a infinito
# para uma gaussiana de média mu e variância sigma2.

# Constantes
mu = 2.0      # Valor médio desejado das variáveis aleatórias gaussianas
sigma2 = 3.0  # Valor desejado da variância da gaussiana
n = 10**3     # Número de variáveis z a sortear para definir cada X
N = 10**5     # Número de variáveis X a obter
Bins = 10**2  # Número de caixas do histograma das variáveis X

# Parâmetros da distribuição uniforme correspondente
a, b = mu - sqrt(3*sigma2), mu + sqrt(3*sigma2)

# Vamos produzir as N variáveis X, cada uma envolvendo a soma de n
# variáveis z distribuídas uniformemente entre a e b correspondendo
# à média mu e à variância sigma2. Depois vamos produzir um gráfico em
# histograma, para visualizar a distribuição dos X obtidos.
X_lista = []
for i in range(N):
    S = 0.0
    for j in range(n):
        S += (b-a)*random() + a
    X_lista.append(mu + sqrt(n)*(S/n-mu))

print("Tempo de execução (em segundos):",(time.time() - tempo_inicial))

# Vamos também produzir o gráfico de uma distribuição gaussiana com
# média mu e variância sigma2, para comparação.
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
            

