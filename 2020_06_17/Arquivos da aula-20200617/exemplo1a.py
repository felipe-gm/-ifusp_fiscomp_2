from math import exp,sqrt
from numpy import random,empty,full,sum
import matplotlib.pyplot as plt

# Implementa o cálculo da magnetização como função do campo externo,
# a uma temperatura fixa, para o modelo de Ising em 1D, através de
# uma simulação de Monte Carlo com o algoritmo de Metropolis.
# O estado de cada spin é representado por uma variável s[i] que pode
# assumir valores -1 ou +1. São utilizadas condições de contorno
# periódicas.

# Parâmetros
L = 2**7             # Número de spins
T = 1.0              # Temperatura da simulação
J = 1.0              # Constante da interação entre spins
Bmax = 1.00          # Valor máximo do campo
Bmin = -1.00         # Valor mínimo do campo
nB = 21              # Número de valores do campo entre Bmin e Bmax
Passos = 1*10**4     # Número total de passos de Monte Carlo a realizar
Equilibra = 10**3    # Número de passos de Monte Carlo para equilibração

# Criamos listas para armazenar os valores do campo (em ordem
# decrescente), da magnetização por spin e do erro em sua estimativa.
dB = (Bmax-Bmin)/(nB-1)
B_lista = [(Bmax - iB*dB) for iB in range(nB)]
m_sim, m_erro = [], []   # Magnetização por spin e seu erro

# Para implementar as condições de contorno periódicas, vamos criar
# vetores com os vizinhos de cada sítio.
ve, vd = empty(L,int), empty(L,int)
for i in range(L):
    ve[i] = i-1     # Vizinho à esquerda
    if i == 0:      # Correção para o primeiro spin
        ve[i] = L-1
    vd[i] = i+1     # Vizinho à direita
    if i == L-1:    # Correção para o último spin
        vd[i] = 0

# Simulação. Definimos uma função para implementar um único passo
# de Monte Carlo. Buscamos coletar médias apenas após a equilibração.
s = full(L,1)   # Spins inicialmente estão todos no estado +1

def passo_MC(s,boltz):  # Função que implementa um passo de Monte Carlo
    z = random.rand(L)  # Sorteamos L números aleatórios
    # Damos a cada spin a chance de mudar de estado, mas aqui fazemos
    # isso aleatoriamente para minimizar correlações entre mudanças.
    for i in random.permutation(L):
        si = s[i]                # Estado do i-ésimo spin
        sv = s[ve[i]] + s[vd[i]] # Soma dos estados dos vizinhos
        if z[i] < boltz[si,sv]:  # Mudança foi aceita?
            s[i] *= -1           # Registramos a mudança
    return s

medidas = Passos - Equilibra     # Número de medidas no equilíbrio
for iB in range(nB):
    B = B_lista[iB]
    # A variação de energia quando um spin é invertido depende apenas do campo
    # e da soma dos estados dos seus vizinhos. Vale a pena criar uma lista
    # para armazenar os valores possíveis do fator de Boltzmann correspondente.
    boltz = full([3,5],0.0)
    for si in range(-1,2):       # O estado do spin vai de -1 a 1
        for sv in range(-2,3):   # Soma dos estados dos vizinhos vai de -2 a 2
            dE = 2*(J*sv + B)*si # Variação da energia na inversão do spin
            boltz[si,sv] = exp(-dE/T)
    # Primeiro aguardamos os passos de equilibração
    for passo in range(Equilibra):
        s = passo_MC(s,boltz)       
    # Vamos começar a calcular médias
    acumula_M, acumula_M2 = 0.0, 0.0
    # Os passos restantes servem para cálculo das médias
    for passo in range(medidas):
        s = passo_MC(s,boltz)
        M = sum(s)
        acumula_M += M        # Acumulamos a magnetização 
        acumula_M2 += M**2    # Acumulamos a magnetização quadrática
    # Agora registramos as médias
    M_medio = acumula_M/medidas
    M2_medio = acumula_M2/medidas
    m_sim.append(M_medio/L)
    m_erro.append(sqrt((M2_medio-M_medio**2)/L**2/(medidas-1)))

# Traçamos os gráficos
plt.rcParams['xtick.labelsize'] = 24
plt.rcParams['ytick.labelsize'] = 24
plt.rcParams['axes.labelsize'] = 28
plt.rcParams['text.usetex'] = True

plt.figure(figsize=(12,9))
plt.errorbar(B_lista,m_sim,m_erro,fmt='g-')
plt.plot(B_lista,m_sim,'go')
plt.ylabel("$m$")
plt.xlabel("$B$")
plt.axvline(0,-1,1,color='k',linewidth=0.5)
plt.axhline(0,-1,1,color='k',linewidth=0.5)
plt.savefig('exemplo1a.png', transparent=True)
plt.show()

