from math import sqrt,pi
from random import random
from numpy import linspace
import matplotlib.pyplot as plt

# Parâmetros de integração
a = -1.0        # Limite inferior de integração
b = 1.0         # Limite superior de integração
h = 1.0         # Altura do retângulo de integração
N = 10**4       # Número de pontos a utilizar para a integração

# Definição da função a ser integrada
def f(x):
    return sqrt(1-x*x)
	
# Laço de integração. Note que contamos os pontos sorteados e aqueles
# abaixo da curva da função, produzindo listas dos valores correspondentes.
k = 0
x_acima, x_abaixo, y_acima, y_abaixo = [], [], [], []
for n in range(N):
    x = a+(b-a)*random()   # Sorteamos a coordenada x do ponto
    y = h*random()         # Sorteamos a coordenada y do ponto
    if y < f(x):   # O ponto sorteado está abaixo da curva da função?
        k += 1	   # Atualizamos o contador desses pontos
        x_abaixo.append(x)
        y_abaixo.append(y)
    else:          # O ponto sorteado está acima da curva da função?
        x_acima.append(x)
        y_acima.append(y)

integral = k/N*(b-a)*h   # Estimativa para a integral
print("A estimativa para pi é ",2*integral)

# Parâmetros da exibição dos gráficos
plt.rcParams['xtick.labelsize'] = 28
plt.rcParams['ytick.labelsize'] = 28
plt.rcParams['axes.labelsize'] = 32

# Traçando o gráfico, incluindo a função
abscissa, funcao = linspace(a,b,100), []
funcao[:] = [f(x) for x in abscissa]
plt.figure(figsize=(12,9))
plt.xlabel("$x$")
plt.ylabel("$f(x)$")
plt.plot(abscissa,funcao)
plt.plot(x_abaixo,y_abaixo,"b.")
plt.plot(x_acima,y_acima,"k.")
plt.show()
