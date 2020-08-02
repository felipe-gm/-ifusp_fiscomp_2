#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np


# # Parâmetros de integração

# In[2]:


N = 10**8       # Número de pontos a utilizar para a integração


# Queremos calcular 
# $$
# I = \int_0^1 \frac{dx}{\sqrt[2]{x}\sqrt[3]{1-x}}
# $$
# a abordagem para utilizar a amostragem por importância para lidar com 
# ambas as divergências é dividir a integral em dois intervalos, cada um
# contendo uma das divergencias
# $$
# \begin{align}
# I &= \int_0^{\frac{1}{2}} \frac{dx}{\sqrt[2]{x}\sqrt[3]{1-x}}
#    + \int_{\frac{1}{2}}^1 \frac{dx}{\sqrt[2]{x}\sqrt[3]{1-x}} \\
#   &= \int_0^{\frac{1}{2}} \frac{dx}{\sqrt[2]{x}\sqrt[3]{1-x}}
#    + \int_0^{\frac{1}{2}} \frac{dx}{\sqrt[2]{1-x}\sqrt[3]{x}}
# \end{align}
# $$
# então, para cada uma das integrais, usamos o método de Monte Carlo do
# valor médio e a amostragem por importância.

# # Definição das funções cujo valor médio será calculado

# In[3]:


# Note que aqui se trata da função de integração dividida por
# w(x)=x**(-1/2)
def f0(x):
    return 1/((1-x)**(1/3))
# Note que aqui se trata da função de integração dividida por
# w(x)=x**(-1/3)
def f1(x):
    return 1/((1-x)**(1/2))


# In[4]:


print(f"O valor exato da integral entre x=0 e x=1 é {2.58711}...")
print("")


# # integração pelo método do valor médio com amostragem por importância

# In[5]:


rng = np.random.default_rng()
z = rng.random((N,))


# In[6]:


intw0 = 2*(1/2)**(1/2)    # Integral de w(x)=x**(-1/2) entre x=0 e x=0.5
intw1 = (3*(1/2)**(2/3))/2# Integral de w(x)=x**(-1/3) entre x=0 e x=0.5

fsw_medio0 = np.mean(f0((z**2)/2))
fsw_medio1 = np.mean(f1((z**(3/2))/2))


# # Estimativa para a integral I

# In[7]:


integral = fsw_medio0*intw0 + fsw_medio1*intw1


# In[8]:


print("A estimativa para a integral pelo método do valor médio é",integral)
print("")

