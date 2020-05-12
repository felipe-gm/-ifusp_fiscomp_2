from random import random

import matplotlib.pyplot as plt

# Variaveis globais

k_lista = [i+1 for i in range(16)]

def C(k): return sum([sequencia[i]*sequencia[i+k] for i in range(N-k)])/(N-k)

# Primeiro ponto

N = int(1e6)
a = 57
c = 1
m = 256
x0 = 10

x = x0
sequencia = []
for i in range(N):
    sequencia.append(x)
    x = ((a*x*m+c)%m)/m

ck1_lista = [C(k) for k in k_lista]

# Segundo ponto

sequencia = [random() for _ in range(N)]

ck2_lista = [C(k) for k in k_lista]

plt.figure(figsize=(12,9))
plt.xlabel("k")
plt.ylabel("$C(k)$")
plt.plot(k_lista, ck1_lista,'b.', label="Gerador do exemplo 2")
plt.plot(k_lista, ck2_lista,'r.', label="Funcao random() do pacote random")
plt.legend(loc='upper right')
plt.show()