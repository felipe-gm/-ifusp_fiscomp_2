from math import sqrt
from random import random
from numpy import empty
import time

tempo_inicial = time.time()

# Parâmetros de integração
ndim = 10       # Número de dimensões da hiperesfera
N = 10**7       # Número de pontos a utilizar para a integração

# Definição da função a ser integrada
def f(r):
    soma = 0.0
    for i in range(ndim-1):
        soma += r[i]**2
    if soma >= 1.0:
        return 0.0          # A função é nula fora da "calota"
    else:
        return sqrt(1-soma)

# Laço de integração pelo método do valor médio.
# Note que a função a ser integrada retorna zero se o ponto sorteado
# em ndmin-1 dimensões está fora da "calota", por isso podemos utilizar
# como intervalo de integração todo o "hiperquadrante".
r = empty(ndim-1)
f_medio = 0.0
for n in range(1,N+1):
    for i in range(ndim-1):
        r[i] = random()   # Sorteamos cada coordenada do ponto
    f_medio += f(r)
f_medio /= N
integral = f_medio # Estimativa para a integral
hipervolume = (2**ndim)*integral # Multiplicamos pelo número de "hiperquadrantes"
print("A estimativa para o volume pelo método do valor médio é ",hipervolume)

print("Tempo de execução em segundos: ",time.time() - tempo_inicial)

