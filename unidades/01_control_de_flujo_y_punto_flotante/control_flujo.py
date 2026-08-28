#! /usr/bin/env python3

#############################
# Condiciales
#############################


###############################################
if False:
    print("Dentro de if")

num = 100

### comparadores; <, > ==, <= ,>=

if  0.0 < num:
    print(f"{num} es positivo")
    a = 10 + 1

elif num < 0.0:
    print(f"{num} es negativo")

else:
    print(f"{num} es cero")


#############################
# Ciclos
#############################

while False:
    print("Estamos en ciclo while")
        
i = 0

##while i <= num:
##    print(i)
##    i += 1

i = 0
suma = 0

while i <= num:
    suma +=  i
    i +=  1

print(suma)










