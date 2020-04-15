#!/usr/bin/env python
# coding: utf-8


# Related third party imports



import matplotlib.pyplot as plt

from numpy import arange, array, concatenate
from numpy import linalg as LA

from vpython import sphere,vector,color,rate,scene,attach_trail


# Local application/library specific imports.



from item2 import *


# Variaveis globais



scene.lights = [] # Removendo todas as luzes da cena
scene.width = 1200 # Ajustando a largura da cena
scene.height =1000 # Ajustando a altura da cena
scene.center = vector(0,0,0) # Ajustando o centro da cena

estrela1 = sphere(color=color.red,radius=150./1e4, make_trail=True) # Esfera para representar o Sol
estrela1.emissive = True

estrela2 = sphere(color=color.yellow,radius=200./1e4, make_trail=True) # Esfera para representar o Sol
estrela2.emissive = True

estrela3 = sphere(color=color.blue,radius=250./1e4, make_trail=True) # Esfera para representar o Sol
estrela3.emissive = True


# Constantes



C = .1





# Condicoes iniciais (e.g.: r(a))



ra = array([ 3., 1.,.0,.0,
            -1.,-2.,.0,.0,
            -1., 1.,.0,.0 ],float)
r = ra
t = a
h_atual = h

x1_lista, y1_lista = [], []
x2_lista, y2_lista = [], []
x3_lista, y3_lista = [], [] 

while t<=b:
    rate(C/h)
    estrela1.pos = vector(r[0],r[1],0)
    estrela2.pos = vector(r[4],r[5],0)
    estrela3.pos = vector(r[8],r[9],0)
    dr, h_atual, h_prox = passo_adapt_extloc(f,r,t,h,prec)
    t, r = t + h_atual, r + dr
    h = min(1e-4,h_prox)
