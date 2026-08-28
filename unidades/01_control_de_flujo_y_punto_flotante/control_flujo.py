#!/usr/bin/env python3

###############################################
# Condicionales: if / elif / else
###############################################

# Clasificamos el estado de movimiento de una partícula según su velocidad
velocidad = -3.5   # m/s (negativa: se mueve en dirección -x)

### comparadores: <, >, ==, <=, >=

if velocidad > 0.0:
    print(f"La partícula se mueve hacia adelante a {velocidad} m/s")
elif velocidad < 0.0:
    print(f"La partícula se mueve hacia atrás a {velocidad} m/s")
else:
    print("La partícula está en reposo")

# Clasificamos su energía cinética respecto a un umbral (p. ej. energía de escape)
masa = 2.0             # kg
energia_cinetica = 0.5 * masa * velocidad**2
energia_umbral = 5.0   # J

if energia_cinetica >= energia_umbral:
    print(f"E_cinetica = {energia_cinetica} J supera el umbral: la partícula escapa")
elif energia_cinetica > 0.0:
    print(f"E_cinetica = {energia_cinetica} J es positiva pero no alcanza el umbral")
else:
    print("La partícula no tiene energía cinética")

###############################################
# Ciclos: while
###############################################

# Un auto frena por fricción con desaceleración constante;
# usamos while para saber cuánto recorre hasta detenerse
velocidad_inicial = 20.0   # m/s
desaceleracion = 4.0       # m/s^2 (debida a la fricción)
dt = 0.5                   # s

velocidad_actual = velocidad_inicial
distancia_total = 0.0
tiempo = 0.0

while velocidad_actual > 0.0:
    distancia_total += velocidad_actual * dt
    velocidad_actual -= desaceleracion * dt
    tiempo += dt

print(f"El auto se detuvo tras recorrer {distancia_total:.2f} m en {tiempo:.2f} s")

# while True con condición de paro: una pelota rebota y en cada bote
# recupera solo una fracción de la altura anterior (coeficiente de restitución),
# hasta que el rebote es demasiado pequeño para importar
altura_bote = 10.0        # m, altura del primer bote
coef_restitucion = 0.6    # fracción de altura que se recupera en cada bote
altura_minima = 0.1       # m, umbral para dejar de contar rebotes
num_bote = 0

while True:
    num_bote += 1
    altura_bote *= coef_restitucion
    print(f"Bote {num_bote}: altura alcanzada = {altura_bote:.3f} m")

    if altura_bote < altura_minima:
        print(f"La pelota dejó de rebotar de forma apreciable tras {num_bote} botes")
        break   # condición de paro: sin este break, el ciclo sería infinito

###############################################
# Ciclos: for, range(), break y continue
###############################################

# Tabla de altura vs. tiempo de una pelota en caída libre,
# muestreada a pasos de tiempo fijos con for y range()
altura_inicial = 50.0   # m
gravedad = 9.8          # m/s^2
dt = 0.5                # s
num_pasos = 20

print("Tabla de altura vs tiempo:")
for paso in range(num_pasos):
    tiempo = paso * dt
    altura = altura_inicial - 0.5 * gravedad * tiempo**2

    if altura <= 0.0:
        print(f"t = {tiempo:.1f} s: la pelota ya tocó el suelo, detenemos la simulación")
        break   # ya no tiene sentido seguir calculando después de tocar el suelo

    if paso % 4 != 0:
        continue   # solo mostramos 1 de cada 4 pasos, para no saturar la salida

    print(f"t = {tiempo:.1f} s -> altura = {altura:.2f} m")
