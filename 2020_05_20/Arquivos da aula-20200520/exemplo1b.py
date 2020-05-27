from math import sqrt,pi
from random import random
from numpy import linspace
import matplotlib.pyplot as plt

# Parâmetros de integração
a = -1.0        # Limite inferior de integração
b = 1.0         # Limite superior de integração
h = 1.0         # Altura do retângulo de integração
N = 10**8       # Número máximo de pontos a utilizar para a integração
dN = 10**3      # A cada quantos pontos estimar a integral 

# Definição da função a ser integrada
def f(x):
    return sqrt(1-x*x)
	
# Laço de integração. Note que contamos os pontos abaixo da curva da função.
k, erro_lista = 0, []
for n in range(1,N+1):
    x = a+(b-a)*random()   # Sorteamos a coordenada x do ponto
    y = h*random()         # Sorteamos a coordenada y do ponto
    if y < f(x):   # O ponto sorteado está abaixo da curva da função?
        k += 1	   # Atualizamos o contador desses pontos
    if n%dN == 0:
        erro_lista.append(abs(pi-2*(b-a)*h*k/n)/pi)

integral = k/N*(b-a)*h   # Estimativa final para a integral
print("A estimativa para pi é ",2*integral)

# Parâmetros da exibição dos gráficos
plt.rcParams['xtick.labelsize'] = 24
plt.rcParams['ytick.labelsize'] = 24
plt.rcParams['axes.labelsize'] = 26

# Traçando o gráfico do erro como função do número de sorteios
plt.figure(figsize=(12,9))
abscissa, funcao = linspace(1,N,100), []
funcao[:] = [1/sqrt(x) for x in abscissa]
plt.xlabel("$n$")
plt.ylabel("Erro relativo")
plt.loglog(abscissa,funcao)
plt.loglog(range(1,N+1,dN),erro_lista,'r.')
plt.show()
