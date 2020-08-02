from math import sin
from random import random
from numpy import linspace
import time

tempo_inicial = time.time()

# Parâmetros de integração
a = 1e-10       # Limite inferior de integração
b = 1.0         # Limite superior de integração
N = 10**6       # Número de pontos a utilizar para a integração

# Definição da função a ser integrada
def f(x):
    return sin(1/x)**2

print("O valor exato da integral entre x=0 e x=1 é",0.673457)
print("")

# Integrando pela regra do trapézio
h = (b-a)/N
integral = 0.5*f(a) + 0.5*f(b)
for k in range(1,N):
    integral += f(a+k*h)
integral *= h
tempo = time.time() - tempo_inicial
print("A estimativa para a integral pela regra do trapézio é ",integral)
print("O tempo de cálculo em segundos foi",tempo)
print("")

tempo_inicial = time.time()

# Laço de integração pelo método de "jogar pedras".
k = 0
h = 1.0
for n in range(1,N+1):
    x = a+(b-a)*random()   # Sorteamos a coordenada x do ponto
    y = random()           # Sorteamos a coordenada y do ponto
    if y < f(x):           # O ponto sorteado está abaixo da curva da função?
        k += 1             # Atualizamos o contador desses pontos
integral = h*(b-a)*k/N     # Estimativa para a integral
tempo = time.time() - tempo_inicial
print("A estimativa para a integral pelo método de jogar pedras é ",integral)
print("O tempo de cálculo em segundos foi",tempo)
print("")

tempo_inicial = time.time()

# Laço de integração pelo método do valor médio.
f_medio = 0.0
for n in range(1,N+1):
    x = a+(b-a)*random()   # Sorteamos a coordenada x do ponto
    f_medio += f(x)
f_medio /= N
integral = f_medio*(b-a)   # Estimativa para a integral
tempo = time.time() - tempo_inicial
print("A estimativa para a integral pelo método do valor médio é ",integral)
print("O tempo de cálculo em segundos foi",tempo)
print("")


