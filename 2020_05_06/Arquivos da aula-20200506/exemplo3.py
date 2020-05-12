import matplotlib.pyplot as plt

# Parâmetros da exibição dos gráficos
plt.rcParams['xtick.labelsize'] = 28
plt.rcParams['ytick.labelsize'] = 28
plt.rcParams['axes.labelsize'] = 32
plt.rcParams['text.usetex'] = True

N = 2**8
a = 1664525
c = 1013904223
m = 4294967296
x0 = 1

x = x0
x1, x2 = [], []
for i in range(N):
    x1.append(x)
    x = (a*x+c)%m
    x2.append(x)
    if x==x0:
        print(i)

plt.figure(figsize=(12,9))
plt.xlabel("$j$")
plt.ylabel("$x_j$")
plt.plot(x2,"bo")
plt.plot(x2,"b-")
plt.show()

plt.figure(figsize=(12,9))
plt.xlabel("$x_j$")
plt.ylabel("$x_{j+1}$")
plt.plot(x1,x2,"o")
plt.show()
