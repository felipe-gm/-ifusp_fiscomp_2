from math import sin,cos,pi
from numpy import arange, array, empty, shape
import matplotlib.pyplot as plt
import time

tempo_inicial = time.time()

# Parâmetros da exibição dos gráficos
plt.rcParams['xtick.labelsize'] = 22
plt.rcParams['ytick.labelsize'] = 22
plt.rcParams['axes.labelsize'] = 26
plt.rcParams['axes.titlesize'] = 30

# Constantes
M = 1.989e30                    # massa do Sol
G = 6.674e-11                   # constante gravitacional
prec = 1000./(365*24*3600)      # precisão da solução (1 km/ano)
nmax = 8        # Número máximo de iterações do algoritmo de Bulirsch-Stoer
a = 0.0                     # Início do intervalo de integração
b = 2*49*365*24*3600.       # Final do intervalo de integração
H = 1*365*24*3600.         # Tamanho sugerido do passo de integração

# Condições iniciais
x_0, vx_0, y_0, vy_0 = 4.0e12, 0.0, 0.0, 5.0e2

GM = G*M
def f(r,t):
    x, vx, y, vy = r[0], r[1], r[2], r[3]
    GMr3 = GM/(x**2+y**2)**1.5
    f0, f1, f2, f3 = vx, -GMr3*x, vy, -GMr3*y
    return array([f0,f1,f2,f3],float)

def passo_mbs_indiv(f,r,t,H,prec,nmax):
    # Calcula um passo no método de Bulirsch-Stoer. Caso se atinjam
    # 'nmax' iterações sem convergência, retorna a informação
    # Inicializamos com um passo do método do ponto médio modificado
    # A matriz R1 armazena a primeira linha da tabela de extrapolação.
    # Por agora, essa linha contém apenas a estimativa do método do
    # ponto médio modificado para a solução no final do intervalo.
    converge = False
    n = 1
    y = r + 0.5*H*f(r,t)
    x = r + H*f(y,t+0.5*H)
    R1 = empty([1,r.shape[0]],float)
    R1[0] = 0.5*(y + x + 0.5*H*f(x,t+H))
    # Agora fazemos um laço aumentando o valor de n até que a precisão
    # seja atingida.
    for n in range(2,nmax+1):
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
        erro = max(abs(epsilon[0]),abs(epsilon[2]))
        if erro <= H*prec:
            converge = True
            break
        R1 = R2
    # Fazemos r igual à estimativa mais precisa de que dispomos
    r = R2[n-1]
    return converge, r  # Retornamos o NOVO VALOR de r

def passo_mbs_adapt(f,H,prec,nmax,r_lista,t_lista):
    # Calcula um passo no método de Bulirsch-Stoer adaptativo.
    # Esta função não retorna nenhum valor, mas apenas atualiza as listas
    # de t e de r ao atingir convergência em até 'nmax' iterações
    r = r_lista[-1]
    t = t_lista[-1]
    converge, r = passo_mbs_indiv(f,r,t,H,prec,nmax)
    if converge == False: # Se não houve convergência, divida o passo por 2
        passo_mbs_adapt(f,H/2,prec,nmax,r_lista,t_lista)
    else:
        t_lista.append(t+H)
        r_lista.append(r)

def integ_mbs_adapt(f,r_a,a,b,H,prec,nmax,r_lista,t_lista):
    # Esta função percorre o intervalo de integração, determinando os
    # valores de r e t com passo máximo de tamanho H, que é subdividido
    # caso não se atinja a precisão requerida em até 'nmax' iterações
    # do algoritmo de Bulirsch-Stoer. A função não retorna um valor,
    # mas atualiza as listas de r e t.
    t = a
    t_lista.append(t)    # Registramos o valor inicial de t
    r_lista.append(r_a)  # Registramos o valor inicial de r
    while t < b:
        passo_mbs_adapt(f,H,prec,nmax,r_lista,t_lista)
        t = t_lista[-1]  # Atualizamos t para o último valor calculado

# Invocando a integração e criando as listas para traçar os gráficos
r_a = array([x_0, vx_0, y_0, vy_0],float)    # Condição inicial
r_lista, t_lista = [], []
integ_mbs_adapt(f,r_a,a,b,H,prec,nmax,r_lista,t_lista)

x_lista, y_lista = [], []
h_lista = []
for i in range(len(t_lista)):
    x_lista.append(r_lista[i][0])
    y_lista.append(r_lista[i][2])
    h_lista.append(t_lista[i]-t_lista[i-1])
h_lista[0] = h_lista[1]

print("Tempo de execução (em segundos):",(time.time() - tempo_inicial))

plt.figure(figsize=(12,9))
plt.title("Bulirsch-Stoer adaptativo")
plt.xlabel("x")
plt.ylabel("y")
plt.plot(x_lista,y_lista,"b.")
plt.show()

plt.figure(figsize=(12,9))
plt.title("Bulirsch-Stoer adaptativo")
plt.xlabel("t")
plt.ylabel("h(t)")
plt.semilogy(t_lista,h_lista)
plt.show()
