from math import exp
from random import random,randrange,seed
from numpy import full
from vpython import cylinder,sphere,vector,scene

# Implementa a solução aproximada do problema da cobertura por dímeros
# utilizando recozimento simulado. 

# Parâmetros do problema
L = 50          # Lado da rede quadrada
T0 = float(L)   # Temperatura inicial do recozimento
Tmin = 1e-5     # Temperatura final do recozimento
tau = 1e5       # Tempo característico da redução da temperatura

# Definindo uma matriz para armazenar a que dímero pertencem os
# sítios da rede. O estado -1 indica que o sítio não pertence a nenhum dímero.
d = full([L,L],-1)
# Definindo uma lista que armazena os sítios pertencentes a cada dímero
sitios = []  

N = 0           # Inicialmente não há dímeros

# Laço principal
t = 0    # Zerando o tempo inicial
T = T0   # Inicializando a temperatura
# Vamos criar listas para indicar os deslocamentos até os 4 vizinhos de qualquer sítio
dx, dy = [1,0,-1,0], [0,1,0,-1]  
#seed(10)  # Retire o comentário para gerar sempre a mesma solução
while T>Tmin:
    # Escolhendo dois sítios vizinhos ao acaso
    x, y = randrange(L), randrange(L)
    viz = randrange(4)
    xv, yv = x+dx[viz], y+dy[viz]
    if xv == -1 or xv == L or yv == -1 or yv == L:
        continue
    # Resfriamento
    t += 1
    T = T0*exp(-t/tau)
    if d[x,y] != d[xv,yv]:        # Sítios não pertencem ao mesmo dímero?
        continue
    if d[x,y] == -1 and d[xv,yv] == -1:  # Sítios estão ambos vazios?
        sitios.append([[x,y],[xv,yv]]) # Acrescentamos o par à lista de dímeros
        d[x,y] = d[xv,yv] = N
        N += 1
        continue
    if d[x,y] != -1 and d[x,y] == d[xv,yv]: # Sítios formam um dímero?
        if random()<exp(-1/T):
            # Se o movimento for aceito, e o dímero não ocupava a última
            # posição na lista, movemos na lista o último dímero para
            # a posição do dímero que será removido e trocamos o rótulo
            indice = d[x,y]
            d[x,y] = d[xv,yv] = -1
            if indice < N-1:
                sitios[indice] = sitios[N-1]
                x, y = sitios[indice][0][0], sitios[indice][0][1]
                xv, yv = sitios[indice][1][0], sitios[indice][1][1]
                d[x,y] = d[xv,yv] = indice
            sitios.pop()
            N -= 1
            continue

# Produzindo uma visualização do resultado
R = 0.2
scene.center = vector(L/2,L/2,0)
for i in range(N):
    sitio0 = vector(sitios[i][0][0],sitios[i][0][1],0)
    sitio1 = vector(sitios[i][1][0],sitios[i][1][1],0)
    sphere(pos=sitio0,radius=R)
    sphere(pos=sitio1,radius=R)
    cylinder(pos=sitio0,axis=sitio1-sitio0,radius=R/2)
                
print("A estimativa para o número de dímeros é",N)
print("A solução ótima é",int(L*L/2))
