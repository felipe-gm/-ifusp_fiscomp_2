"""Related third party imports"""

import matplotlib.pyplot as plt

from numpy import arange, array, sin, cos, arccos, exp, empty


# Variaveis globais

a = .0        # Inicio do intervalo da variavel independente 
b = 100.      # Final do intervalo da variavel independente
H = 4         # Tamanho inicial de um passo de integracao
prec = 1e-8   # Precisao desejada do passo
nmax = 8      # Número máximo de iterações do algoritmo de Bulirsch-Stoer

# Parâmetros da exibição dos gráficos
plt.rcParams['xtick.labelsize'] = 22
plt.rcParams['ytick.labelsize'] = 22
plt.rcParams['axes.labelsize'] = 26
plt.rcParams['axes.titlesize'] = 30


# Constantes

K, M, RHO, F_EXT, OMEGA_EXT = 1., 1., .4, 1., .1


# # Equacoes de diferenca

def f(r,t):
    x, y = r[0], r[1]
    fx, fy = y, (-K*x -RHO*y +F_EXT*cos(OMEGA_EXT*t))/M
    return array([fx, fy], float)


def solucao_analitica(t):
    return x_trans(t) + x_est(t)

def x_est(t):
    omega_0 = (K/M)**(1/2)
    gamma = RHO/M

    A_ext = F_EXT/(
        M
        *(
            (omega_0**2 - OMEGA_EXT**2)**2
            + gamma**2 * OMEGA_EXT**2
        )**(1/2)
    )

    phi = -arccos(
        (omega_0**2 - OMEGA_EXT**2)
        /(
            (omega_0**2 - OMEGA_EXT**2)**2
            + gamma**2 * OMEGA_EXT**2
        )**(1/2)
    )

    return A_ext * cos(OMEGA_EXT*t + phi)

def x_trans(t):
    omega_0 = (K/M)**(1/2)
    gamma = RHO/M

    omega = (omega_0**2 - (gamma**2)/4)**(1/2)

    alpha, beta = (1.-x_est(a)), (-1. + gamma*(1.-x_est(a))/2.)/omega

    return exp(-gamma*t/2)*(alpha*cos(omega*t) + beta*sin(omega*t))


# # Integracao numerica

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
        erro = abs(epsilon[0])
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

# Condicoes iniciais (e.g.: r(a))
x_0, vx_0 = 1, -1

# Invocando a integração e criando as listas para traçar os gráficos
r_a = array([x_0, vx_0],float)    # Condição inicial
r_lista, t_lista = [], []
integ_mbs_adapt(f,r_a,a,b,H,prec,nmax,r_lista,t_lista)

x_lista, h_lista = [], []
 
for i in range(len(t_lista)):
    x_lista.append(r_lista[i][0])
    h_lista.append(t_lista[i]-t_lista[i-1])
h_lista[0] = h_lista[1]


x_exato = [] 
t_exato = arange(a,b,1e-2)
for t in t_exato:  # Criando a lista com a solução exata
    x_exato.append(solucao_analitica(t))


plt.figure(figsize=(12,9))
plt.title("Bulirsch-Stoer adaptativo")
plt.plot(t_exato, x_exato, label='solucao exata')
plt.plot(t_lista, x_lista, 'b.', label='solucao numerica')
plt.plot(t_lista, h_lista, 'g.', label='tamanho de cada passo')
plt.xlabel("t")
plt.ylabel("x(t), h(t)")
plt.legend(loc='upper right')
plt.show()