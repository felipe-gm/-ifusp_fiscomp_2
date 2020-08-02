from math import exp
from numpy import random,linspace,full
import matplotlib.pyplot as plt

# Implementa o cálculo da energia interna para uma coleção de
# sistemas de dois níveis (0 e epsilon) não interagentes.
# Vamos simular a coleção primeiro na temperatura mais alta,
# de modo a decidir quantos passos de Monte Carlo são necessários
# para a equilibração

# Parâmetros
N = 100          # Número de sistemas compondo a coleção
T = 1.0          # Temperatura 
Passos = 100     # Número total de passos de Monte Carlo a realizar
    
# Simulação. 
s = full(N,1)           # Partículas inicialmente no estado excitado
boltz = exp(-1.0/T)     # Probabilidade de transição exp(-dE/T)
E_sim = []              # Lista para armazenar os valores da energia
for passo in range(Passos):
    z = random.rand(N)  # Sorteamos N números aleatórios
    # Damos a cada partícula sequencialmente a chance de mudar de estado
    for i in range(N):
        dE = 1 - 2*s[i] # Variação da energia em caso de mudança
        if dE < 0 or z[i] < boltz:  # Mudança foi aceita?
            s[i] = 1 - s[i]         # Registramos a mudança
    energia = sum(s)            # Energia ao final do passo
    E_sim.append(energia/N)

# Traçamos o gráfico
plt.rcParams['xtick.labelsize'] = 28
plt.rcParams['ytick.labelsize'] = 28
plt.rcParams['axes.labelsize'] = 32

plt.figure(figsize=(12,9))
plt.plot(E_sim)
plt.ylabel("Energia interna por partícula")
plt.xlabel("Passo de Monte Carlo")
plt.show()

