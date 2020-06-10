from math import exp
from random import random,randrange,seed
from vpython import sphere,curve,canvas,vector,mag

# Implementa a solução aproximada do problema do caixeiro viajante
# utilizando recozimento simulado. 

# Parâmetros do problema
N = 25          # Número de cidades a percorrer
T0 = float(N)   # Temperatura inicial do recozimento
Tmin = N/25000  # Temperatura final do recozimento
tau = 1e3       # Tempo característico da redução da temperatura

# Parâmetros da visualização
R = 0.02        # Raio das esferas que representam as cidades

# Inicializando os gráficos
canvas(center=vector(0.5,0.5,0))

# Função para calcular a distância percorrida durante a viagem
def distancia(r):
    s = 0.0
    for i in range(N):
        s += mag(r[i+1]-r[i])
    return s

# Retire o comentário da linha abaixo caso queira que as cidades tenham
# sempre as mesmas posições
seed(11)

# Escolhendo N posições para as cidades e calculando a distância inicial
r = [0]*(N+1) # Lista que armazena as posições das cidades
for i in range(N):
    r[i] = vector(random(), random(), 0.0)
r[N] = r[0]
D = distancia(r)

print("A estimativa original é",D)

# Criando os objetos da visualização
for i in range(N):
    sphere(pos=r[i],radius=R)
l = curve(pos=r,radius=R/2)           # Criando a curva que liga as cidades

# Retire o comentário da linha abaixo caso queira que os movimentos
# aceitos sejam sempre os mesmos
seed(40)

# Laço principal
t = 0    # Zerando o tempo inicial
T = T0   # Inicializando a temperatura
while T>Tmin:
    # Resfriamento
    t += 1
    T = T0*exp(-t/tau)
    # Atualizando a visualização a cada 100 tentativas de movimento
    if t%100==0:
        l.clear()                     # Limpando a curva
        l = curve(pos=r,radius=R/2)   # Recriando a curva
    # Escolhendo duas cidades distintas e trocando a ordem de visitação
    i,j = randrange(1,N),randrange(1,N)
    while i==j:
        i,j = randrange(1,N),randrange(1,N)
    D_antigo = D
    r[i],r[j] = r[j],r[i]
    D = distancia(r)        # Nova distância total da viagem
    deltaD = D - D_antigo   # Variação na distância
    # Se o movimento é rejeitado, trocamos as cidades de volta
    if random()>exp(-deltaD/T):
        r[i],r[j] = r[j],r[i]
        D = D_antigo

print("A estimativa para a distância otimizada é",D)
l.clear()                     # Limpando a curva
l = curve(pos=r,radius=R/2)   # Recriando a curva

