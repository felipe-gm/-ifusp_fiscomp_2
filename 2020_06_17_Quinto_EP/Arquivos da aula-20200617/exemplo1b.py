from math import exp,sqrt
from numpy import random,empty,full,sum
import matplotlib.pyplot as plt

# Implementa o cálculo da energia interna, do calor específico e da
# magnetização em função da temperatura, a um campo magnético fixo,
# para o modelo de Ising em 1D, através de uma simulação de Monte Carlo
# com o algoritmo de Metropolis. O estado de cada spin é representado
# por uma variável s[i] que pode assumir valores -1 ou +1. São utilizadas
# condições de contorno periódicas.

# Parâmetros
L = 2**7             # Número de spins
B = 0.0              # Campo magnético da simulação
J = 1.0              # Constante da interação entre spins
Tmax = 2.00          # Valor máximo da temperatura
Tmin = 0.10          # Valor mínimo da temperatura
nT = 20              # Número de valores da temperatura entre Tmin e Tmax
Passos = 12*10**3    # Número total de passos de Monte Carlo a realizar
Equilibra = 4*10**3  # Número de passos de Monte Carlo para equilibração
medidas = Passos - Equilibra # Número de medidas no equilíbrio
Bloco = medidas//100         # No. de passos de MC para estimar erro no c.e.

# Criamos listas para armazenar os valores da temperatura (em ordem
# decrescente), da magnetização por spin e do erro em sua estimativa.
dT = (Tmax-Tmin)/(nT-1)
T_lista = [(Tmax - iT*dT) for iT in range(nT)]
m_sim, m_erro = [], []   # Magnetização por spin e seu erro
E_sim, C_sim = [], []    # Energia e calor específico segundo a simulação
E_erro, C_erro = [], []  # Erros na energia e no calor específico

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

# Definimos uma função para calcular a energia de uma configuração
def Energia(s):
    soma = 0.0
    for i in range(L):
        soma += -J*s[i]*s[vd[i]] -B*s[i]
    return soma

# Simulação. Definimos uma função para implementar um único passo
# de Monte Carlo. Buscamos coletar médias apenas após a equilibração.
s = 1-2*random.randint(2, size=L)   # Estados iniciais aleatórios dos spins

def passo_MC(s,boltz):  # Função que implementa um passo de Monte Carlo
    z = random.rand(L)  # Sorteamos L números aleatórios
    # Damos a cada spin a chance de mudar de estado, mas aqui fazemos
    # isso aleatoriamente para minimizar correlações entre mudanças.
    for i in random.permutation(L):
        bf = boltz[s[i],s[ve[i]] + s[vd[i]]]
        if z[i] < bf or (bf == boltz[s[i],0] and z[i]<0.5):# Mudança foi aceita?
            s[i] = -s[i]                     # Registramos a mudança
    return s

for iT in range(nT):
    E_plot = []
    T = T_lista[iT]
    print("Simulando à temperatura T=",T)
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
        energia = Energia(s)
        E_plot.append(energia)    
    # Vamos começar a calcular médias
    acumula_M, acumula_M2 = 0.0, 0.0
    acumula_E, acumula_E2, acumula_C, acumula_C2 = 0.0, 0.0, 0.0, 0.0
    for jext in range(medidas//Bloco):
        acumula_Ece, acumula_E2ce = 0.0, 0.0 # Para estimar c.e. e seu erro
        for jint in range(Bloco):
            s = passo_MC(s,boltz)
            M = abs(sum(s))       # Note o uso do valor absoluto!
            acumula_M += M        # Acumulamos a magnetização 
            acumula_M2 += M**2    # Acumulamos a magnetização quadrática
            energia = Energia(s)
            acumula_Ece += energia        # Acumulamos a energia
            acumula_E2ce += energia**2    # Acumulamos a energia quadrática
            E_plot.append(energia)
        Ece_medio = acumula_Ece/Bloco
        E2ce_medio = acumula_E2ce/Bloco
        dc = (E2ce_medio-Ece_medio**2)/L/T**2
        acumula_C += dc
        acumula_C2 += dc**2
        acumula_E += acumula_Ece      # Acumulamos a energia 
        acumula_E2 += acumula_E2ce    # Acumulamos a energia quadrática    
    # Retire os comentários das duas linhas logo abaixo para ver a curva da
    # energia em função dos passos de MC e avaliar a equilibração.
    #plt.plot(E_plot)
    #plt.show()
    # Agora registramos as médias
    M_medio = acumula_M/medidas
    M2_medio = acumula_M2/medidas
    m_sim.append(M_medio/L)
    m_erro.append(sqrt((M2_medio-M_medio**2)/L**2/(medidas-1)))
    E_medio = acumula_E/medidas
    E2_medio = acumula_E2/medidas
    E_sim.append(E_medio/L)
    E_erro.append(sqrt(abs((E2_medio-E_medio**2)/L**2/(medidas-1))))
    C_medio = acumula_C/(medidas//Bloco)
    C2_medio = acumula_C2/(medidas//Bloco)
    C_sim.append((E2_medio-E_medio**2)/L/T**2)
    C_erro.append(sqrt(abs((C2_medio-C_medio**2)/(medidas//Bloco-1))))

# Traçamos os gráficos
plt.rcParams['xtick.labelsize'] = 24
plt.rcParams['ytick.labelsize'] = 24
plt.rcParams['axes.labelsize'] = 28
plt.rcParams['text.usetex'] = True

plt.figure(figsize=(12,9))
plt.errorbar(T_lista,E_sim,E_erro,fmt='r-')
plt.plot(T_lista,E_sim,'ro')
plt.ylabel("$E/L$")
plt.xlabel("$T$")
plt.show()

plt.figure(figsize=(12,9))
plt.errorbar(T_lista,C_sim,C_erro,fmt='r-')
plt.plot(T_lista,C_sim,'ro')
plt.ylabel("$c$")
plt.xlabel("$T$")
plt.xlim(0,2)
plt.ylim(0,0.5)
plt.savefig('exemplo1b_ce.png', transparent=True)
plt.show()

plt.figure(figsize=(12,9))
plt.errorbar(T_lista,m_sim,m_erro,fmt='r-')
plt.plot(T_lista,m_sim,'ro')
plt.ylabel("$m$")
plt.xlabel("$T$")
plt.xlim(0,2)
plt.ylim(0,1.1)
plt.savefig('exemplo1b_mag.png', transparent=True)
plt.show()

