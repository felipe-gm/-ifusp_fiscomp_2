from math import exp
from numpy import random,empty,full,sum
import matplotlib.pyplot as plt

# Implementa a simulação do modelo de Ising em 2D, através do método de
# Monte Carlo com o algoritmo de Metropolis.
# Este código fixa a temperatura e o campo, exibindo a configuração
# dos spins ao final da simulação.

# Parâmetros
L = 2**7             # Lado da rede quadrada
B = 0.0              # Campo magnético da simulação
J = 1.0              # Constante da interação entre spins
T = 3.0              # Valor da temperatura
Passos = 2000        # Número total de passos de Monte Carlo a realizar

N = L*L              # Número de spins

# Para implementar as condições de contorno periódicas, vamos criar
# uma matriz com os vizinhos de cada sítio.
v = empty((4,N),int)
for i in range(N):
    v[0,i] = i+1    # Vizinho à direita
    if (i+1)%L == 0:      # Correção para os spins na última coluna
        v[0,i] -= L
    v[1,i] = i+L    # Vizinho acima
    if v[1,i] > N-1:      # Correção para os spins na última linha
        v[1,i] -= N       
    v[2,i] = i-1    # Vizinho à esquerda
    if i%L == 0:          # Correção para os spins na primeira coluna
        v[2,i] += L       
    v[3,i] = i-L    # Vizinho abaixo
    if v[3,i] < 0:        # Correção para os spins na primeira linha
        v[3,i] += N

# Definimos uma função para calcular a energia de uma configuração
def Energia(s):
    soma = 0.0
    for i in range(N):
        soma += -(J*(s[v[0,i]]+s[v[1,i]]) + B)*s[i]
    return soma

# Simulação. Definimos uma função para implementar um único passo
# de Monte Carlo. Buscamos coletar médias apenas após a equilibração.
s = 1-2*random.randint(2, size=N)   # Estados iniciais aleatórios dos spins

def passo_MC(s,boltz):  # Função que implementa um passo de Monte Carlo
    z = random.rand(N)  # Sorteamos L números aleatórios
    # Damos a cada spin a chance de mudar de estado, mas aqui fazemos
    # isso aleatoriamente para minimizar correlações entre mudanças.
    for i in random.permutation(N):
        bf = boltz[s[i],s[v[0,i]] + s[v[1,i]] + s[v[2,i]] + s[v[3,i]]]
        if z[i] < bf or (bf == boltz[s[i],0] and z[i]<0.5):# Mudança foi aceita?
            s[i] = -s[i]                     # Registramos a mudança
    return s

print("Simulando à temperatura T=",T)
# A variação de energia quando um spin é invertido depende apenas do campo
# e da soma dos estados dos seus vizinhos. Vale a pena criar uma lista
# para armazenar os valores possíveis do fator de Boltzmann correspondente.
boltz = full([3,9],0.0)
for si in range(-1,2):       # O estado do spin vai de -1 a 1
    for sv in range(-4,5):   # Soma dos estados dos vizinhos vai de -4 a 4
        dE = 2*(J*sv + B)*si # Variação da energia na inversão do spin
        boltz[si,sv] = exp(-dE/T)
# Percorremos os passos de MC. Ao final, traçamos a energia e a
# magnetizaçãoem função do passo, para verificar que a configuração
# final é representativa do equilíbrio térmico
E_plot, M_plot = [], []
for passo in range(Passos):
    s = passo_MC(s,boltz)
    energia = Energia(s)/N
    mag = abs(sum(s))/N
    E_plot.append(energia)
    M_plot.append(mag)

plt.plot(E_plot)
plt.show()

plt.plot(M_plot)
plt.show()

# Agora exibimos a configuração
xp, yp = [], []
xn, yn = [], []
for i in range(N):
    x, y = i%L, i//L
    if s[i] == 1:
        xp.append(x)
        yp.append(y)
    else:
        xn.append(x)
        yn.append(y)

# Parâmetros da exibição dos gráficos
plt.rcParams['xtick.labelsize'] = 28
plt.rcParams['ytick.labelsize'] = 28
plt.rcParams['axes.labelsize'] = 32

plt.figure(figsize=(10,10))
tamanho = (540/L)**2
plt.scatter(xp,yp,c='b',marker='s',s=tamanho)
plt.scatter(xn,yn,c='r',marker='s',s=tamanho)
plt.savefig('exemplo2b.png', transparent=True)
plt.show()  
