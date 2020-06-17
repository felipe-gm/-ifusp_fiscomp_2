from math import exp,sqrt
from numpy import random,linspace,full
import matplotlib.pyplot as plt

# Implementa o cálculo da energia interna por partícula e do
# calor específico (c.e.) em função da temperatura para uma coleção de
# sistemas de dois níveis (0 e epsilon) não interagentes.
# O cálculo é feito tanto a partir das expressões analíticas
# exatas quanto por meio de simulações de Monte Carlo com o
# algoritmo de Metropolis. Trabalhamos com unidades em que epsilon
# e a constante de Boltzmann são ambas iguais a 1. Vamos simular
# a coleção primeiro na temperatura mais alta, e vamos resfriá-la
# progressivamente, em cada etapa aguardando um certo tempo de
# equilibração antes de acumular valores para o cálculo das médias.

# Parâmetros
N = 100              # Número de sistemas compondo a coleção
Tmin = 0.1           # Temperatura mínima 
Tmax = 2             # Temperatura máxima
nT = 20              # Número de valores de temperatura entre Tmin e Tmax
Passos = 500         # Número total de passos de Monte Carlo a realizar
Equilibra = 100       # Número de passo de Monte Carlo para equilibração
M = Passos-Equilibra # Número de amostras para estimar médias
Bloco = M//40        # No. de passos de MC para estimar erro no c.e.

# Criamos listas para armazenar os valores da temperatura (em ordem
# decrescente), da energia interna média e do calor específico.
dT = (Tmax-Tmin)/(nT-1)
T_lista = [(Tmax - i*dT) for i in range(nT)]
E_sim, C_sim = [], []   # Energia e calor específico segundo a simulação
E_erro, C_erro = [], [] # Erros na energia e no calor específico
# Listas para os resultados analíticos
T_exato = linspace(Tmin,Tmax,200)
E_exato, C_exato = [], []
for T in T_exato:
    boltz = exp(-1.0/T)                 # Fator de Boltzmann
    uboltz = 1.0 + boltz
    E_exato.append(boltz/uboltz)
    C_exato.append(boltz/(T*uboltz)**2)
    
# Simulação. Definimos uma função para implementar um único passo
# de Monte Carlo. Buscamos coletar médias apenas após a equilibração.
s = full(N,1)   # Partículas inicialmente no estado excitado

def passo_MC(s,boltz):  # Função que implementa um passo de Monte Carlo
    z = random.rand(N)  # Sorteamos N números aleatórios
    # Damos a cada partícula sequencialmente a chance de mudar de estado
    for i in range(N):
        dE = 1 - 2*s[i] # Variação da energia em caso de mudança
        if dE < 0 or z[i] < boltz:  # Mudança foi aceita?
            s[i] = 1 - s[i]         # Registramos a mudança
    return(s)

for n in range(nT):
    T = T_lista[n]
    boltz = exp(-1.0/T)     # Probabilidade de transição exp(-dE/T)
    # Primeiro aguardamos os passos de equilibração
    for passo in range(Equilibra):
        s = passo_MC(s,boltz)       
    # Vamos começar a calcular médias
    acumula_E, acumula_E2, acumula_C, acumula_C2 = 0.0, 0.0, 0.0, 0.0
    # Os passos restantes servem para cálculo das médias
    for jext in range(M//Bloco):
        acumula_Ece, acumula_E2ce = 0.0, 0.0 # Para estimar c.e. e seu erro
        for jint in range(Bloco):
            s = passo_MC(s,boltz)
            energia = sum(s)
            acumula_Ece += energia        # Acumulamos a energia
            acumula_E2ce += energia**2    # Acumulamos a energia quadrática
        Ece_medio = acumula_Ece/Bloco
        E2ce_medio = acumula_E2ce/Bloco
        dc = (E2ce_medio-Ece_medio**2)/N/T**2
        acumula_C += dc
        acumula_C2 += dc**2
        acumula_E += acumula_Ece      # Acumulamos a energia 
        acumula_E2 += acumula_E2ce    # Acumulamos a energia quadrática
    # Agora registramos as médias
    E_medio = acumula_E/M
    E2_medio = acumula_E2/M
    E_sim.append(E_medio/N)
    E_erro.append(sqrt((E2_medio-E_medio**2)/N**2/(M-1)))
    C_medio = acumula_C/(M//Bloco)
    C2_medio = acumula_C2/(M//Bloco)
    C_sim.append((E2_medio-E_medio**2)/N/T**2)
    C_erro.append(sqrt((C2_medio-C_medio**2)/(M//Bloco-1)))
    

# Traçamos os gráficos
plt.rcParams['xtick.labelsize'] = 28
plt.rcParams['ytick.labelsize'] = 28
plt.rcParams['axes.labelsize'] = 32

plt.figure(figsize=(12,9))
plt.plot(T_exato,E_exato)
plt.errorbar(T_lista,E_sim,E_erro,fmt='ro')
plt.ylabel("Energia interna por partícula")
plt.xlabel("Temperatura")
plt.show()

plt.figure(figsize=(12,9))
plt.plot(T_exato,C_exato)
plt.errorbar(T_lista,C_sim,C_erro,fmt='ro')
plt.ylabel("Calor específico")
plt.xlabel("Temperatura")
plt.show()
