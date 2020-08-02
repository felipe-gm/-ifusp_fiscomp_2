from math import exp
from numpy import random,linspace,full
import matplotlib.pyplot as plt

# Implementa o cálculo da energia interna por partícula e do
# calor específico em função da temperatura para uma coleção de
# sistemas de dois níveis (0 e epsilon) não interagentes.
# O cálculo é feito tanto a partir das expressões analíticas
# exatas quanto por meio de simulações de Monte Carlo com o
# algoritmo de Metropolis. Trabalhamos com unidades em que epsilon
# e a constante de Boltzmann são ambas iguais a 1. Vamos simular
# a coleção primeiro na temperatura mais alta, e vamos resfriá-la
# progressivamente, em cada etapa aguardando um certo tempo de
# equilibração antes de acumular valores para o cálculo das médias.

# Parâmetros
N = 100          # Número de sistemas compondo a coleção
Tmin = 0.1       # Temperatura mínima 
Tmax = 2         # Temperatura máxima
nT = 20          # Número de valores de temperatura entre Tmin e Tmax
Passos = 2000     # Número total de passos de Monte Carlo a realizar
Equilibra = 10   # Número de passo de Monte Carlo para equilibração
M = Passos-Equilibra # Número de amostras para estimar médias

# Criamos listas para armazenar os valores da temperatura (em ordem
# decrescente), da energia interna média e do calor específico.
dT = (Tmax-Tmin)/(nT-1)
T_lista = [(Tmax - i*dT) for i in range(nT)]
E_sim, C_sim = [], []   # Energia e calor específico segundo a simulação
# Listas para os resultados analíticos
T_exato = linspace(Tmin,Tmax,200)
E_exato, C_exato = [], []
for T in T_exato:
    boltz = exp(-1.0/T)                 # Fator de Boltzmann
    uboltz = 1.0 + boltz
    E_exato.append(boltz/uboltz)
    C_exato.append(boltz/(T*uboltz)**2)
    
# Simulação. Há um laço externo que percorre as temperaturas e outro
# laço interno que percorre os passos de Monte Carlo, coletando médias
# apenas após os primeiros 10 passos.
s = full(N,1)   # Partículas inicialmente no estado excitado
for n in range(nT):
    T = T_lista[n]
    boltz = exp(-1.0/T)     # Probabilidade de transição exp(-dE/T)
    # Primeiro aguardamos os passos de equilibração
    for passos in range(Equilibra):
        z = random.rand(N)  # Sorteamos N números aleatórios
        # Damos a cada partícula sequencialmente a chance de mudar de estado
        for i in range(N):
            dE = 1 - 2*s[i] # Variação da energia em caso de mudança
            if dE < 0 or z[i] < boltz:  # Mudança foi aceita?
                s[i] = 1 - s[i]         # Registramos a mudança
    # Vamos começar a calcular médias
    acumula_E, acumula_E2 = 0.0, 0.0
    # Os passos restantes servem para cálculo das médias
    for passos in range(Equilibra,Passos):
        z = random.rand(N)  # Sorteamos N números aleatórios
        # Damos a cada partícula sequencialmente a chance de mudar de estado
        for i in range(N):
            dE = 1 - 2*s[i] # Variação da energia em caso de mudança
            if dE < 0 or z[i] < boltz:  # Mudança foi aceita?
                s[i] = 1 - s[i]         # Registramos a mudança
        energia = sum(s)            # Energia ao final do passo
        acumula_E += energia        # Acumulamos a energia
        acumula_E2 += energia**2    # Acumulamos a energia quadrática
    # Agora registramos as médias
    E_medio = acumula_E/M
    print(T,E_medio)
    E2_medio = acumula_E2/M
    E_sim.append(E_medio/N)
    C_sim.append((E2_medio-E_medio**2)/N/T**2)

# Traçamos os gráficos
plt.rcParams['xtick.labelsize'] = 28
plt.rcParams['ytick.labelsize'] = 28
plt.rcParams['axes.labelsize'] = 32

plt.figure(figsize=(12,9))
plt.plot(T_exato,E_exato)
plt.plot(T_lista,E_sim,'ro')
plt.ylabel("Energia interna por partícula")
plt.xlabel("Temperatura")
plt.show()

plt.figure(figsize=(12,9))
plt.plot(T_exato,C_exato)
plt.plot(T_lista,C_sim,'ro')
plt.ylabel("Calor específico")
plt.xlabel("Temperatura")
plt.show()
