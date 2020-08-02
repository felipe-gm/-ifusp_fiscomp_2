#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt


# # Parâmetros do problema

# In[2]:


T0 = 10      # Temperatura inicial do recozimento
Tmin = 1e-3  # Temperatura final do recozimento
tau = 1e4     # Tempo característico da redução da temperatura

passos = int(-tau*np.log(Tmin/T0))


# # Função para calcular a energia

# In[3]:


def f(x): return x**2 - np.cos(4*np.pi*x)


# # Resfriamento

# In[4]:


T = T0*np.exp(-np.arange(0, passos)/tau) # Inicializando a temperatura


# # Escolha dos deltas em potencial

# In[5]:


mu, sigma = 0, 1 # mean and standard deviation
delta = np.random.default_rng().normal(mu, sigma, passos)


# # Sinal de rejeição

# In[6]:


random = np.random.random(passos)


# # Laço principal

# In[7]:


x = np.empty(passos+1)
x[0] = 2 # Condição inicial


# In[8]:


for t in range(passos):
    f1 = f(x[t]+delta[t])   # Energia em x(t)+delta(x)
    delta_f = f1 - f(x[t])  # Variação na energia
    if random[t]>np.exp(-delta_f/T[t]): 
        x[t+1] = x[t]
    else :  # Se o movimento é aceito, aplicamos delta(t)
        x[t+1] = x[t] + delta[t]


# In[9]:


plt.scatter(range(passos+1), x, s=.1, label='x(t)')
plt.title("gráficos do valor de x em função do tempo")
plt.ylabel("x")
plt.xlabel("tempo")
plt.legend(loc="upper right")
plt.show()


# In[10]:


print("A estimativa para x minimizando f do item 1 é",x[-1])

