# Related third party imports.
import numpy as np
import matplotlib.pyplot as plt

# Variaveis globais
mu, sigma = .0, 1. # mean and standard deviation

# Constantes
T = 10**3   # Numero de passos da simulacao

# Criando listas para armazenar a trajetoria
x_lista, y_lista = [], []

# A particula sai da origem
x, y = 0, 0

x_lista.append(x)
y_lista.append(y)

# Draw samples from the distribution
s = np.random.default_rng().normal(mu, sigma, (2,T))

caminhada = np.cumsum(s,axis=1)      # sum over columns for each of s

x_lista.extend(caminhada[0])
y_lista.extend(caminhada[1])

# Tracando o grafico
plt.figure(figsize=(12,9))
plt.plot(x_lista,y_lista,label='trajetoria')
plt.plot(x_lista[0],y_lista[0],'o',label='inicio')
plt.plot(x_lista[-1],y_lista[-1],'o',label='fim')
plt.xlabel("x(t)")
plt.ylabel("y(t)")
plt.title('caminhada com 1000 passos')
plt.legend(loc='upper right')
plt.show()