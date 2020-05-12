from random import randrange
from numpy import arange
import matplotlib.pyplot as plt

# Constantes
L = 100     # L+1, com L par, é o lado da rede quadrada
N = 10**3   # Número de passos da simulação

# Criando listas para armazenar a trajetória
x_lista, y_lista = [], []

# A partícula sai da origem
x, y = 0, 0
n = 0
while n<N:
    x_lista.append(x)
    y_lista.append(y)
    z = randrange(4)# Sorteamos uma de quatro direções
    if z == 0:
        dr = [1,0]  # Passo para a direita
    if z == 1:
        dr = [0,1]  # Passo para cima
    if z == 2:
        dr = [-1,0] # Passo para a esquerda
    if z == 3:
        dr = [0,-1] # Passo para baixo
    # Não permitimos que a partícula saia da rede
    if (abs(x+dr[0])<L/2) and (abs(y+dr[1])<L/2):
        x += dr[0]
        y += dr[1]
        n += 1

# Parâmetros da exibição dos gráficos
plt.rcParams['xtick.labelsize'] = 28
plt.rcParams['ytick.labelsize'] = 28
plt.rcParams['axes.labelsize'] = 32

# Traçando o gráfico
plt.figure(figsize=(12,9))
plt.plot(x_lista,y_lista,'r-')
plt.plot(x_lista[0],y_lista[0],'ko')
plt.plot(x_lista[-1],y_lista[-1],'bs')
plt.xlim(-L/2,L/2)
plt.ylim(-L/2,L/2)
plt.xlabel("$x(t)$")
plt.ylabel("$y(t)$")
plt.show()
            

