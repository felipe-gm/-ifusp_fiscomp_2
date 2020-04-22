#!/usr/bin/env python
# coding: utf-8


# Related third party imports

import matplotlib.pyplot as plt

from numpy import arange, array, concatenate
from numpy import linalg as LA

from vpython import box,helix,sphere,vector,color,textures,rate,scene


# Local application/library specific imports.

from Questao2_2 import *


# Variaveis globais

b = 200.      # Final do intervalo da variavel independente

scene.width = 1200 # Ajustando a largura da cena
scene.height =1000 # Ajustando a altura da cena
scene.camera.pos = vector(2,-1,2)
scene.camera.axis = vector(-2,1,-2)

suporte = box(
    pos=vector(.0,.05,.0),
    up=vector(.0,1.,.0),
    length=2,
    height=.1,
    width=2,
    texture=textures.wood
) # Suporte

mola = helix(
    pos=vector(0,0,0),
    axis=vector(*r_vrl),
    color=color.yellow,
    radius=.03
) # Mola do pendulo

bola = sphere(
    axis=vector(*r_vrl),
    radius=.1,
    make_trail=True
) # Bola


# Constantes

C = 1.


while t<=b:
    rate(C/h)
    mola.axis = vector(*r_vrl)
    bola.pos = vector(*r_vrl)
    r_vrl, v_vrl, vpm = passo_vrl(f,r_vrl,v_vrl,vpm,t,h)
    t += h
