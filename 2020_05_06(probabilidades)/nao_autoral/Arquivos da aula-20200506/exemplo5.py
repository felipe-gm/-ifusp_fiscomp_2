from random import randrange,seed

seed(181) # Diferentes argumentos geram diferentes sequências
print("Vamos gerar um número aleatório.")
z = randrange(10)
print("Um segundo número aleatório é",randrange(10))
print("O primeiro número aleatório foi",z)


