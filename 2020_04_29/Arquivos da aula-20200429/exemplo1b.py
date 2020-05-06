# Problema de valores de contorno: que velocidade inicial faz com que uma
# bola lançada para cima no instante 'ta' retorne à altura inicial no
# instante 'tb'?
from math import sin,cos,exp,sqrt,log
from numpy import arange,array
import time

tempo_inicial = time.time()

# Constantes
g = 9.81        # Aceleração da gravidade

ta = 0.0          # Início do intervalo da variável independente
tb = 10.0         # Final do intervalo da variável independente
ya = 0.0          # Valor da variável independente no início do intervalo
yb = 0.0          # Valor da variável independente no final do intervalo
dt = 1e-2         # Tamanho do passo de integração
prec = 1e-10      # Precisão do resultado para a velocidade inicial

def f(r,t):
    y, v = r[0], r[1]
    f0, f1 = v, -g
    return array([f0,f1],float)

def passo_rk4(f,r,t,dt):            # Calcula um passo no método de RK4
    k1 = dt*f(r,t)
    k2 = dt*f(r+0.5*k1,t+0.5*dt)
    k3 = dt*f(r+0.5*k2,t+0.5*dt)
    k4 = dt*f(r+k3,t+dt)
    return (k1+2.0*(k2+k3)+k4)/6.0

# Função que integra a equação diferencial definida por f entre 'ta' e 'tb'
# com condição inicial 'ra' e retorna o valor final do vetor 'r'
def r_final(f,ta,tb,ra,dt):
    r = ra
    for t in arange(ta,tb,dt):
        r += passo_rk4(f,r,t,dt)
    return r

# Solução via método da secante
v1 = 1   # Primeiro palpite de uma velocidade 
v2 = 2   # Segundo palpite de uma velocidade
r1 = array([ya,v1],float)   # Vetor 'r' inicial correspondente a 'v1'
altura2 = r_final(f,ta,tb,r1,dt)[0] - yb # Altura final relativa com 'v1'
while abs(v2-v1) > prec:
    altura1 = altura2
    r2 = array([ya,v2],float)
    altura2 = r_final(f,ta,tb,r2,dt)[0] - yb # Altura final relativa com 'v2' 
    v1, v2 = v2, v2 - altura2*(v2-v1)/(altura2-altura1)
v = v2       # Resultado final do cálculo

print("Tempo de execução (em segundos):",(time.time() - tempo_inicial))
print("A velocidade inicial necessária é",v,"m/s")

