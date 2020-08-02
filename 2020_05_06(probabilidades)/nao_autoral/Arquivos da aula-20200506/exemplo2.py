import matplotlib.pyplot as plt

# Parâmetros da exibição dos gráficos
plt.rcParams['xtick.labelsize'] = 28
plt.rcParams['ytick.labelsize'] = 28
plt.rcParams['axes.labelsize'] = 32
plt.rcParams['text.usetex'] = True

N = 256
a = 57
c = 1
m = 256
x0 = 10

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
