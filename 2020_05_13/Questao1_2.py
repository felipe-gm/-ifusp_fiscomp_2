# Related third party imports.
import numpy as np
import matplotlib.pyplot as plt

# Variaveis globais
mu, sigma = .0, 1. # mean and standard deviation

# Constantes
T = 10**3   # Numero de passos da simulacao
M = int(1e4)     # Numero de caminhadas independentes

# Criando listas para armazenar a trajetoria
y_lista = []

# A particula sai da origem
y = 0

y_lista.append(y)

# Draw samples from the distribution
s = np.random.default_rng().normal(mu, sigma, (2,T,M))

caminhada = np.cumsum(s,axis=1)      # sum over columns for each of s

sqrd_dist = caminhada[0]**2+caminhada[1]**2

mean_sqrd_dist = np.mean(sqrd_dist, axis=1)

y_lista.extend(mean_sqrd_dist**(1/2))

# Tracando o grafico

fig, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel('t')
ax1.set_ylabel('R(t)', color=color)
ax1.plot(y_lista, color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('R^2(t)', color=color)  # we already handled the x-label with ax1
ax2.plot(mean_sqrd_dist, color=color)
ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()