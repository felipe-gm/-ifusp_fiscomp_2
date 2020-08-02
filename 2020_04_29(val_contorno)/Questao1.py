"""
Determine de que altura deve ser abandonada a bola.

Funcao do script: usar o metodo do chute para determinar de que altura
deve ser abandonada uma bola (a partir do repouso) para que atinja o 
solo apos 3 segundos. Supondo que durante a queda a bola esta sujeita 
nao apenas a forca gravitacional, mas tambem a uma forca de resistencia 
do ar proporcional a velocidade instantanea, ou seja, que a altura y da 
bola dependa do tempo de acordo com a equacao de movimento
    (d^2 y)/(dt^2) = -G -(RHO/M)(dy/dt),
Nas unidades do SI, supondo G=9.8, M=.4 e RHO=.1.

Sintaxe da linha de comando:
$ python Questao1.py
"""

# Standard library imports.
import sys

# Related third party imports.
from numpy import arange,array

# Constantes
G, M, RHO = 9.8, .4, .1

ta = .0           # Inicio do intervalo da variavel independente
tb = 3.           # Final do intervalo da variavel independente
va = 0.0          # Valor da derivada da variavel dependente no incio do intervalo
yb = 0.0          # Valor da variavel dependente no final do intervalo
dt = 1e-2         # Tamanho do passo de integracao
prec = 1e-10      # Precisao do resultado para a velocidade inicial

def f(r,t):
    v = r[1]
    f0, f1 = v, -G -(RHO/M)*v
    return array([f0,f1],float)

def passo_rk4(f,r,t,dt):            # Calcula um passo no metodo de RK4
    k1 = dt*f(r,t)
    k2 = dt*f(r+0.5*k1,t+0.5*dt)
    k3 = dt*f(r+0.5*k2,t+0.5*dt)
    k4 = dt*f(r+k3,t+dt)
    return (k1+2.0*(k2+k3)+k4)/6.0

# Funcao que integra a equacao diferencial definida por f entre 'ta' e 'tb'
# com condicao inicial 'ra' e retorna o valor final do vetor 'r'
def r_final(f,ta,tb,ra,dt):
    r = ra
    for t in arange(ta,tb,dt):
        r += passo_rk4(f,r,t,dt)
    return r

# Solucao via metodo da bisseccao (busca binaria)
y1 = 1e-3   # Palpite de uma altura final abaixo da desejada
y2 = 1e3    # Palpite de uma altura final acima da desejada
r1 = array([y1,va],float)   # Vetor 'r' inicial correspondente a 'y1'
r2 = array([y2,va],float)   # Vetor 'r' inicial correspondente a 'y2'
altura1 = r_final(f,ta,tb,r1,dt)[0] - yb # Alt. relativa final partindo com 'y1'
altura2 = r_final(f,ta,tb,r2,dt)[0] - yb # Alt. relativa final partindo com 'y2'
if (altura1*altura2 > 0):
    sys.exit("Alturas relativas com mesmo sinal. Modifique palpites.")
while abs(y2-y1) > prec:
    yp = (y1+y2)/2                           # Media entre 'y1' e 'y2'
    rp = array([yp,va],float)                # Vetor 'r' inicial correspondente
    alturap = r_final(f,ta,tb,rp,dt)[0] - yb # Altura relativa correspondente
    if (altura1 * alturap) > 0:              # Altura final menor que desejada? 
        y1, altura1 = yp, alturap            # Sim; aumentamos o palpite 'y1'
    else:
        y2, altura2 = yp, alturap            # Nao; diminuímos o palpite 'y2'
y = (y1+y2)/2       # Resultado final do calculo

print('Altura inicial necessaria eh: ',y,' m')

