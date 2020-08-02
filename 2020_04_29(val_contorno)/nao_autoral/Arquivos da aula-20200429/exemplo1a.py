# Problema de valores de contorno: que velocidade inicial faz com que uma
# bola lançada para cima no instante 'ta' retorne à altura inicial no
# instante 'tb'?
from math import sin,cos,exp,sqrt,log
from numpy import arange,array
import time
import sys

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

# Solução via método da bissecção (busca binária)
v1 = 1e-3   # Palpite de uma velocidade com altura final abaixo da desejada
v2 = 1e3    # Palpite de uma velocidade com altura final acima da inicial
r1 = array([ya,v1],float)   # Vetor 'r' inicial correspondente a 'v1'
r2 = array([ya,v2],float)   # Vetor 'r' inicial correspondente a 'v2'
altura1 = r_final(f,ta,tb,r1,dt)[0] - yb # Alt. relativa final partindo com 'v1'
altura2 = r_final(f,ta,tb,r2,dt)[0] - yb # Alt. relativa final partindo com 'v2'
if (altura1*altura2 > 0):
    sys.exit("Alturas relativas com mesmo sinal. Modifique velocidades.")
while abs(v2-v1) > prec:
    vp = (v1+v2)/2                           # Média entre 'v1' e 'v2'
    rp = array([ya,vp],float)                # Vetor 'r' inicial correspondente
    alturap = r_final(f,ta,tb,rp,dt)[0] - yb # Altura relativa correspondente
    if (altura1 * alturap) > 0:              # Altura final menor que desejada? 
        v1, altura1 = vp, alturap            # Sim; aumentamos o palpite 'v1'
    else:
        v2, altura2 = vp, alturap            # Não; diminuímos o palpite 'v2'
v = (v1+v2)/2       # Resultado final do cálculo

print("Tempo de execução (em segundos):",(time.time() - tempo_inicial))
print("A velocidade inicial necessária é",v,"m/s")

