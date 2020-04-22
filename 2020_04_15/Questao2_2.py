#!/usr/bin/env python
# coding: utf-8

# Related third party imports

import matplotlib.pyplot as plt

from numpy import arange, array, dot
from numpy import linalg as LA


# Variaveis globais

plt.rcParams['xtick.labelsize'] = 22
plt.rcParams['ytick.labelsize'] = 22
plt.rcParams['axes.labelsize'] = 26
plt.rcParams['axes.titlesize'] = 30

a = .0        # Inicio do intervalo da variavel independente 
b = 20.       # Final do intervalo da variavel independente
h = 2e-2      # Tamanho inicial de um passo de integracao


# Constantes

M, K, G, L0 = .02, .9, 9.8, .2
VERSOR_x = array([1.,.0,.0])
VERSOR_Y = array([.0,1.,.0])
VERSOR_Z = array([.0,.0,1.])


# # Equacoes de diferenca

def f(r,t):
    l_i = LA.norm(r)
    versor_r = r/l_i
    spring_a = -K*(l_i-L0)/M
    
    fx = spring_a*dot(versor_r,VERSOR_x)
    fy = spring_a*dot(versor_r,VERSOR_Y) -G
    fz = spring_a*dot(versor_r,VERSOR_Z)
    
    return array([fx, fy, fz], float)


def passo_vrl(f,r,v,vpm,t,h):   # Calcula um passo no método de Verlet
    r += h*vpm                  # Calculando a posição no próximo t "inteiro"
    k = h*f(r,t+h)
    v = vpm + 0.5*k             # Calculando a velocidade no próximo t "inteiro"
    vpm += k                    # Calculando a velocidade no próximo t "médio"
    return r, v, vpm            # Retorna os NOVOS VALORES de r, v e vpm


# Condicoes iniciais

ra = array([.05,-.25,.0, -.1,-.15,.2],float)
r_vrl, v_vrl = ra[:3], ra[3:]
t = a

vpm = v_vrl + 0.5*h*f(r_vrl,t)  # Inicializando valores no ponto "médio"
def main():
    global t,r_vrl,v_vrl,vpm

    t_lista, e_k_lista, e_pl_lista, e_pg_lista, e_lista = [], [], [], [], []

    while t<= b:   # Realizando a integração numérica
        e_k = M*dot(v_vrl,v_vrl)/2
        e_pl = K*((LA.norm(r_vrl)-L0)**2)/2
        e_pg = M*G*r_vrl[1]
        
        t_lista.append(t)
        e_k_lista.append(e_k)
        e_pl_lista.append(e_pl)
        e_pg_lista.append(e_pg) 
        e_lista.append(e_k+e_pl+e_pg) 
        r_vrl, v_vrl, vpm = passo_vrl(f,r_vrl,v_vrl,vpm,t,h)
        t += h

    plt.figure(figsize=(12,9))
    plt.plot(t_lista, e_k_lista, '.', label='Cinetica')
    plt.plot(t_lista, e_pl_lista, '.', label='Pot. elastica')
    plt.plot(t_lista, e_pg_lista, '.', label='Pot. gravitacional')
    plt.plot(t_lista, e_lista, '.', label='Mec. tot.')
    plt.title("Verlet")
    plt.xlabel("t")
    plt.ylabel("energia")
    plt.legend(loc='upper right')
    plt.show()


if __name__ == "__main__":
    # execute only if run as a script
    main()