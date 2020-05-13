from random import random
from numpy import arange
import matplotlib.pyplot as plt

# Constantes
NBi213 = int(1e5)          # Numero inicial de atomos de Bi213 
NTl = 0             # Numero inicial de atomos de Tl
NPb = 0             # Numero inicial de atomos de Pb
NBi209 = 0             # Numero inicial de atomos de Bi209
tauBi = 46.*60      # Meia-vida do Bi213 em segundos
tauTl = 2.2*60      # Meia-vida do Tl em segundos
tauPb = 3.3*60      # Meia-vida do Pb em segundos
h = 1.             # Tamanho do passo de tempo em segundos
pBi = 1 - 2**(-h/tauBi) # Probabilidade de decaimento de Bi213 em um passo
pTl = 1 - 2**(-h/tauTl) # Probabilidade de decaimento de Tl em um passo
pPb = 1 - 2**(-h/tauPb) # Probabilidade de decaimento de Pb em um passo
tmax = 2e4         # Tempo total da simulacao

# Criando as listas para os graficos
t_lista = arange(0.0,tmax,h)
Bi213_lista, Tl_lista, Pb_lista, Bi209_lista = [], [], [], []

# Laco principal
for t in t_lista:
    Bi213_lista.append(NBi213)
    Tl_lista.append(NTl)
    Pb_lista.append(NPb)
    Bi209_lista.append(NBi209)

    # Dando a cada atomo a chance de decair
    # Item 1
    decaimento = 0
    for i in range(NPb):
        if random() < pPb:
            decaimento += 1
    NPb -= decaimento
    NBi209 += decaimento
    # Item 2
    decaimento = 0
    for i in range(NTl):
        if random() < pTl:
            decaimento += 1
    NTl -= decaimento
    NPb += decaimento
    # Item 3
    decaimentoTl = 0
    decaimentoPb = 0
    for i in range(NBi213):
        if random() < pBi:
            if random() < .9791:
                decaimentoPb += 1
            else:
                decaimentoTl +=1
    NBi213 -= decaimentoPb + decaimentoTl
    NPb += decaimentoPb
    NTl += decaimentoTl

# Tracando o grafico
plt.figure(figsize=(12,9))
plt.semilogy(t_lista,Bi213_lista,label='Bi213')
plt.semilogy(t_lista,Tl_lista,label='Tl')
plt.semilogy(t_lista,Pb_lista,label='Pb')
plt.semilogy(t_lista,Bi209_lista,label='Bi209')
plt.xlabel("Tempo")
plt.ylabel("Numero de atomos")
plt.legend(loc='upper right')
plt.show()
            

