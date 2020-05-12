from random import randrange,seed
from numpy import arange,empty
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Constantes
L = 100     # L+1, com L par, é o lado da rede quadrada
N = 10**6   # Número de passos da simulação

# Criando listas para armazenar a trajetória
data = empty((2,N),dtype=int)

# A partícula sai da origem
x, y = 0, 0
n = 0
seed(52)
while n<N:
    data[0,n]=x
    data[1,n]=y
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

# Atualização do quadro da animação, que contém todos 
# os passos, representados como linhas, até o passo n.
def atualiza_linha(n, data, line):
    line.set_data(data[..., :n])
    return line,

# Traçando o gráfico
fig = plt.figure(figsize=(12,9))
l, = plt.plot([], [], 'r-')
plt.xlim(-L/2,L/2)
plt.ylim(-L/2,L/2)
plt.xlabel("$x(t)$")
plt.ylabel("$y(t)$")
linha_ani = animation.FuncAnimation(fig, \
                                    atualiza_linha, \
                                    fargs=(data, l), \
                                    interval=10,\
                                    blit=True)
plt.show()
            

