#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt


# # Parâmetros do problema

# In[2]:


T0 = 10      # Temperatura inicial do recozimento
Tmin = 1e-2  # Temperatura final do recozimento
tau = 1e4     # Tempo característico da redução da temperatura

passos = int(-tau*np.log(Tmin/T0))


# # Função para calcular a energia

# In[3]:


def f(x): return -(np.cos(x) + np.sin(x*2**(1/2)) + np.cos(x*3**(1/2)))


# # Resfriamento

# In[4]:


T = T0*np.exp(-np.arange(0, passos)/tau) # Inicializando a temperatura


# # Escolha dos deltas em potencial

# In[5]:


mu, sigma = 0, 28 # mean and standard deviation
delta = np.random.default_rng().normal(mu, sigma, passos)


# # Sinal de rejeição

# In[6]:


random = np.random.random(passos)


# # Laço principal

# In[7]:


x, fx = np.empty(passos+1), np.empty(passos+1)
x[0] = 2 # Condição inicial
f0 = fx[0] = f(2)

# In[8]:


for t in range(passos):
    f1 = f((x[t]+delta[t])%50)   # Energia em x(t)+delta(x)
    delta_f = f1 - f0  # Variação na energia
    if random[t]>np.exp(-delta_f/T[t]): 
        x[t+1] = x[t]
        fx[t+1] = f0
    else :  # Se o movimento é aceito, aplicamos delta(t)
        x[t+1] = (x[t] + delta[t])%50
        fx[t+1] = f0 = f1


# In[9]:


# plt.scatter(range(passos+1), x, s=.1, label='x(t)')
# plt.title("gráficos do valor de x em função do tempo")
# plt.ylabel("x")
# plt.xlabel("tempo")
# plt.legend(loc='upper right')
# plt.show()


# In[10]:


print("A estimativa para x maximizando f do item 2 é",x[np.argmin(fx)])

