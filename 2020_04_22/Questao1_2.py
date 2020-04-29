# -*- coding: utf-8 -*-
"""Related third party imports"""

import matplotlib.pyplot as plt

from numpy import arange, array, concatenate, empty
from numpy import linalg as LA


"""Variaveis globais"""

a = 0.0                      # Início do intervalo de integração
b = 60*60*24*370*250         # Final do intervalo de integração
prec = 1e3/60*60*24*365.2425 # Precisão requerida para o cálculo
h = 60*60*24*7               # Tamanho do passo de integração

# Parâmetros da exibição dos gráficos
plt.rcParams['xtick.labelsize'] = 22
plt.rcParams['ytick.labelsize'] = 22
plt.rcParams['axes.labelsize'] = 26
plt.rcParams['axes.titlesize'] = 30


"""Constantes"""

M = 1.9891e30                   # massa do Sol
G = 6.6738e-11                  # constante gravitacional
GM = G*M


def f(r,t):
    x, vx, y, vy = r[0], r[1], r[2], r[3]
    GMr3 = GM/(x**2+y**2)**1.5
    f0, f1, f2, f3 = vx, -GMr3*x, vy, -GMr3*y
    return array([f0,f1,f2,f3],float)

def passo_mbs(f,r,t,H,prec): # Calcula um passo no método de Bulirsch-Stoer
    # Inicializamos com um passo do método do ponto médio modificado
    # A matriz R1 armazena a primeira linha da tabela de extrapolação.
    # Por agora, essa linha contém apenas a estimativa do método do
    # ponto médio modificado para a solução no final do intervalo.
    n = 1
    y = r + 0.5*H*f(r,t)
    x = r + H*f(y,t+0.5*H)
    R1 = empty([1,r.shape[0]],float)
    R1[0] = 0.5*(y + x + 0.5*H*f(x,t+H))
    # Agora fazemos um laço aumentando o valor de n até que a precisão
    # seja atingida.
    erro = 2*H*prec # Garantindo que o laço seja executado ao menos 1 vez
    while erro > H*prec:
        n += 1
        h = H/n
        # Método do ponto médio modificado
        y = r + 0.5*h*f(r,t)
        x = r + h*f(y,t+0.5*h)
        for i in range(n-1):
            y += h*f(x,t+(i+1.0)*h)
            x += h*f(y,t+(i+1.5)*h)
        # Calculando as estimativas por extrapolação.
        # As matrizes R1 e R2 armazenam a penúltima e a última
        # linhas mais recentes da tabela
        R2 = empty([n,r.shape[0]],float)
        R2[0] = 0.5*(y + x + 0.5*h*f(x,t+h))
        for m in range(1,n):
            epsilon = (R2[m-1]-R1[m-1])/((n/(n-1))**(2*m)-1)
            R2[m] = R2[m-1] + epsilon
        erro = abs(epsilon[0])
        R1 = R2
    # Fazemos r igual à estimativa mais precisa de que dispomos
    r = R2[n-1]
    #print("t=",t," n='",n) # Imprime o tempo e o n para convergência
    return r        # Retornamos o NOVO VALOR de r


# Condicoes iniciais (e.g.: r(a))
ra = array([-4.4368e12,.0, .0,6.1218e3],float)
r_mbs = ra
t = a
h_atual = h

x_lista, y_lista = [], []

while t<= b:   # Realizando a integração numérica
    x_lista.append(r_mbs[0])
    y_lista.append(r_mbs[2])
    r_mbs = passo_mbs(f,r_mbs,t,h,prec)
    t += h


plt.figure(figsize=(12,9))
plt.title("Bulirsch-Stoer")
plt.plot([0], [0], 'y.', label=f'Posicao do Sol', markersize=30)
plt.plot(x_lista, y_lista, 'b.', label=f'Orbita Terra', markersize=.1)
plt.xlabel("x")
plt.ylabel("y")
plt.legend(loc='upper right')
plt.show()
