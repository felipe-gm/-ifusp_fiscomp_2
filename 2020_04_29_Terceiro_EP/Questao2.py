from numpy import arange, array, empty, shape
import matplotlib.pyplot as plt

# Variaveis globais

prec = 1e-10    # precisao da solucao
nmax = 8    # Numero maximo de iteracoes do algoritmo de Bulirsch-Stoer
a = .0    # Inicio do intervalo de integracao
b = 20.    # Final do intervalo de integracao
H = 20.    # Tamanho sugerido do passo de integracao

# Parametros da exibicao dos graficos
plt.rcParams['xtick.labelsize'] = 22
plt.rcParams['ytick.labelsize'] = 22
plt.rcParams['axes.labelsize'] = 26
plt.rcParams['axes.titlesize'] = 30

# Constantes

A, B = 1., 3.

# Condicoes iniciais
x_0, y_0 = .0, .0

def f(r,t):
    x, y = r[0], r[1]
    f0, f1 = 1-(B+1)*x+A*(x**2)*y, B*x-A*(x**2)*y
    return array([f0,f1],float)

def passo_mbs_indiv(f,r,t,H,prec,nmax):
    # Calcula um passo no metodo de Bulirsch-Stoer. Caso se atinjam
    # 'nmax' iteracoes sem convergencia, retorna a informacao
    # Inicializamos com um passo do metodo do ponto medio modificado
    # A matriz R1 armazena a primeira linha da tabela de extrapolacao.
    # Por agora, essa linha contem apenas a estimativa do metodo do
    # ponto medio modificado para a solucao no final do intervalo.
    converge = False
    n = 1
    y = r + 0.5*H*f(r,t)
    x = r + H*f(y,t+0.5*H)
    R1 = empty([1,r.shape[0]],float)
    R1[0] = 0.5*(y + x + 0.5*H*f(x,t+H))
    # Agora fazemos um laco aumentando o valor de n ate que a precisao
    # seja atingida.
    for n in range(2,nmax+1):
        h = H/n
        # Metodo do ponto medio modificado
        y = r + 0.5*h*f(r,t)
        x = r + h*f(y,t+0.5*h)
        for i in range(n-1):
            y += h*f(x,t+(i+1.0)*h)
            x += h*f(y,t+(i+1.5)*h)
        # Calculando as estimativas por extrapolacao.
        # As matrizes R1 e R2 armazenam a penultima e a ultima
        # linhas mais recentes da tabela
        R2 = empty([n,r.shape[0]],float)
        R2[0] = 0.5*(y + x + 0.5*h*f(x,t+h))
        for m in range(1,n):
            epsilon = (R2[m-1]-R1[m-1])/((n/(n-1))**(2*m)-1)
            R2[m] = R2[m-1] + epsilon
        erro = max(abs(epsilon[0]),abs(epsilon[1]))
        if erro <= H*prec:
            converge = True
            break
        R1 = R2
    # Fazemos r igual a estimativa mais precisa de que dispomos
    r = R2[n-1]
    return converge, r  # Retornamos o NOVO VALOR de r

def passo_mbs_adapt(f,H,prec,nmax,r_lista,t_lista):
    # Calcula um passo no metodo de Bulirsch-Stoer adaptativo.
    # Esta funcao nao retorna nenhum valor, mas apenas atualiza as listas
    # de t e de r ao atingir convergencia em ate 'nmax' iteracoes
    r = r_lista[-1]
    t = t_lista[-1]
    converge, r = passo_mbs_indiv(f,r,t,H,prec,nmax)
    if converge == False: # Se nao houve convergencia, divida o passo por 2
        passo_mbs_adapt(f,H/2,prec,nmax,r_lista,t_lista)
    else:
        t_lista.append(t+H)
        r_lista.append(r)

def integ_mbs_adapt(f,r_a,a,b,H,prec,nmax,r_lista,t_lista):
    # Esta funcao percorre o intervalo de integracao, determinando os
    # valores de r e t com passo maximo de tamanho H, que e subdividido
    # caso nao se atinja a precisao requerida em ate 'nmax' iteracoes
    # do algoritmo de Bulirsch-Stoer. A funcao nao retorna um valor,
    # mas atualiza as listas de r e t.
    t = a
    t_lista.append(t)    # Registramos o valor inicial de t
    r_lista.append(r_a)  # Registramos o valor inicial de r
    while t < b:
        passo_mbs_adapt(f,H,prec,nmax,r_lista,t_lista)
        t = t_lista[-1]  # Atualizamos t para o último valor calculado

# Invocando a integracao e criando as listas para tracar os graficos
r_a = array([x_0, y_0],float)    # Condicao inicial
r_lista, t_lista = [], []
integ_mbs_adapt(f,r_a,a,b,H,prec,nmax,r_lista,t_lista)

x_lista = [r[0] for r in r_lista]
y_lista = [r[1] for r in r_lista]

plt.figure(figsize=(12,9))
plt.title("Bulirsch-Stoer adaptativo")
plt.xlabel("t")
plt.ylabel("x(t), y(t)")
plt.plot(t_lista,x_lista,".", label='x')
plt.plot(t_lista,y_lista,".", label='y')
plt.legend(loc='upper right')
plt.show()