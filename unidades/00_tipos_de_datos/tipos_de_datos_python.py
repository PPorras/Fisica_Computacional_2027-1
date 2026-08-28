#!/usr/bin/env python3

###############################################
# La función print() y el tipo cadena (strings)
###############################################

# Datos de una pelota que se deja caer desde cierta altura
masa_pelota = 2.0            # kg
velocidad_lanzamiento = 3.0  # m/s (rapidez horizontal de lanzamiento)
altura_inicial = 20.0        # m
gravedad = 9.8                # m/s^2

energia_cinetica = 0.5 * masa_pelota * velocidad_lanzamiento**2
energia_potencial = masa_pelota * gravedad * altura_inicial

print(f"La energía cinética de la pelota es {energia_cinetica} J")
print(f"Su energía potencial a {altura_inicial} m de altura es {energia_potencial} J")
print(f"La energía mecánica total es {energia_cinetica + energia_potencial} J")

print(f"masa_pelota elevada al cuadrado: {masa_pelota**2}")

###############################################
# Los tipos entero (int) y flotante (float)
###############################################

numero_pasos = 8       # int: número de instantes que vamos a simular
paso_temporal = 0.3    # float: dt en segundos entre cada instante

print(f"numero_pasos es del tipo {type(numero_pasos)}")
print(f"paso_temporal es del tipo {type(paso_temporal)}")

tiempo_total = numero_pasos * paso_temporal   # producto
mitad_del_tiempo = tiempo_total / 2           # división real
pasos_en_1_segundo = 1 // paso_temporal       # división entera
resto = numero_pasos % 3                      # residuo (módulo)

print(f"tiempo_total={tiempo_total} s, mitad_del_tiempo={mitad_del_tiempo} s")
print(f"pasos_en_1_segundo={pasos_en_1_segundo}, resto={resto}")

# Conversión entre tipos
tiempo_total_entero = int(tiempo_total)     # trunca la parte decimal
numero_pasos_flotante = float(numero_pasos)

print(f"int(tiempo_total) = {tiempo_total_entero}")
print(f"float(numero_pasos) = {numero_pasos_flotante}")

###############################################
# Operadores de comparación
###############################################

# Comparar dos flotantes regresa un valor booleano (True o False)
cayo_mas_de_15_metros = altura_inicial > 15.0
paso_muy_pequeno = paso_temporal <= 0.5
mitad_menor_que_total = mitad_del_tiempo < tiempo_total
energia_potencial_mayor_igual = energia_potencial >= energia_cinetica
tiempo_total_es_3_segundos = tiempo_total == 3.0
energias_distintas = energia_cinetica != energia_potencial

print(f"¿altura_inicial > 15.0? {cayo_mas_de_15_metros}")
print(f"¿paso_temporal <= 0.5? {paso_muy_pequeno}")
print(f"¿mitad_del_tiempo < tiempo_total? {mitad_menor_que_total}")
print(f"¿energia_potencial >= energia_cinetica? {energia_potencial_mayor_igual}")
print(f"¿tiempo_total == 3.0? {tiempo_total_es_3_segundos}")
print(f"¿energia_cinetica != energia_potencial? {energias_distintas}")

###############################################
# El tipo cadena (str) y sus operaciones
###############################################

nombre_experimento = "Caída libre"
responsable = "Prof. Porras"

titulo_reporte = nombre_experimento + " - " + responsable
separador = "=" * 30

print(separador)
print(titulo_reporte.upper())
print(separador)

print(f"El nombre del experimento tiene {len(nombre_experimento)} caracteres")
print(f"Primeras 5 letras: {nombre_experimento[:5]}")
print(f"Últimas 5 letras: {nombre_experimento[-5:]}")

palabras = nombre_experimento.split()   # separa por espacios -> lista
print(f"palabras = {palabras}")
print("_".join(palabras).lower())

###############################################
# La función type(), las tuplas y listas
###############################################

constantes = (3e8, 6.626e-34, 5.67e-8, 2e-23)  # c, h, sigma, hbar (tupla)

print(f"'constantes' es del tipo {type(constantes)} y contiene {constantes}")
print(f"La velocidad de la luz es: {constantes[0]}")
# Las tuplas son inmutables: la siguiente línea daría error si se descomenta
# constantes[0] = 40

print(f"Número de constantes guardadas: {len(constantes)}")
velocidad_luz, planck, sigma, hbar = constantes   # desempaquetado
print(f"velocidad_luz={velocidad_luz}, planck={planck}")

posicion_inicial = (0.0, 0.0, altura_inicial)   # vector (x, y, z), no cambia
print(f"posicion_inicial (x, y, z) = {posicion_inicial}")

trayectoria = [altura_inicial]   # lista: aquí sí queremos poder modificar/crecer

tiempo = 0.0
for paso in range(numero_pasos):
    tiempo += paso_temporal
    posicion = max(0.0, altura_inicial - 0.5 * gravedad * tiempo**2)
    trayectoria.append(round(posicion, 2))

print(f"'trayectoria' es del tipo {type(trayectoria)}")
print(f"trayectoria = {trayectoria}")

print(f"Altura inicial: {trayectoria[0]} m, última altura registrada: {trayectoria[-1]} m")
trayectoria.append(0.0)          # la pelota se queda en el suelo
trayectoria.sort(reverse=True)   # de mayor a menor altura
print(f"trayectoria ordenada de mayor a menor: {trayectoria}")

###############################################
# El tipo conjunto (set)
###############################################

# En el laboratorio se repite el experimento y se mide el tiempo de caída;
# por redondeo, varias mediciones terminan siendo iguales
tiempos_medidos = [2.0, 2.0, 2.1, 1.9, 2.0, 2.1, 2.0]
tiempos_unicos = set(tiempos_medidos)

print(f"tiempos_unicos es del tipo {type(tiempos_unicos)}")
print(f"tiempos_unicos = {tiempos_unicos}")

# Dos detectores registran ciertas frecuencias de resonancia (en Hz)
frecuencias_detector_1 = {50, 60, 120, 180}
frecuencias_detector_2 = {60, 120, 240}

print(f"Detectadas por ambos: {frecuencias_detector_1 & frecuencias_detector_2}")
print(f"Detectadas por al menos uno: {frecuencias_detector_1 | frecuencias_detector_2}")
print(f"Solo por el detector 1: {frecuencias_detector_1 - frecuencias_detector_2}")
print(f"¿60 Hz fue detectada por el detector 1? {60 in frecuencias_detector_1}")

frecuencias_detector_1.add(300)
print(f"frecuencias_detector_1 tras add(300): {frecuencias_detector_1}")

###############################################
# El tipo diccionario (dict)
###############################################

constantes_fisicas = {
    "velocidad_luz": 3e8,        # m/s
    "planck": 6.626e-34,         # J*s
    "gravitacional": 6.674e-11,  # N*m^2/kg^2
}

print(f"constantes_fisicas es del tipo {type(constantes_fisicas)}")
print(f"La constante de Planck es: {constantes_fisicas['planck']}")

constantes_fisicas["stefan_boltzmann"] = 5.67e-8    # agregamos una nueva llave
constantes_fisicas["planck"] = 6.62607015e-34       # valor actualizado (CODATA)

for nombre_constante, valor in constantes_fisicas.items():
    print(f"{nombre_constante} = {valor}")

propiedades_pelota = {
    "nombre": "pelota de goma",
    "masa": masa_pelota,
    "altura_inicial": altura_inicial,
}

print(f"propiedades_pelota es del tipo {type(propiedades_pelota)}")
print(f"Llaves: {list(propiedades_pelota.keys())}")
print(f"Valores: {list(propiedades_pelota.values())}")

###############################################
# Buleanos
###############################################

si = True
no = False

print(f"La variable si es un:{si} y del tipo {type(si)} ademas en entero es{int(si)} ")
print(f"La variable no es un:{no} y es del tipo {type(no)}")

# Operadores lógicos: and, or, not
pelota_toco_el_suelo = trayectoria[-1] == 0.0
hubo_energia_potencial = energia_potencial > 0.0

experimento_valido = pelota_toco_el_suelo and hubo_energia_potencial
print(f"¿pelota_toco_el_suelo and hubo_energia_potencial? {experimento_valido}")

frecuencia_60hz_detectada = (60 in frecuencias_detector_1) or (60 in frecuencias_detector_2)
print(f"¿60 Hz detectada por al menos un detector (or)? {frecuencia_60hz_detectada}")

pelota_sigue_en_el_aire = not pelota_toco_el_suelo
print(f"¿La pelota sigue en el aire (not)? {pelota_sigue_en_el_aire}")

# Combinando comparadores y operadores lógicos
medicion_confiable = (paso_temporal <= 0.5) and (numero_pasos > 5)
print(f"¿La medición es confiable (paso pequeño y suficientes pasos)? {medicion_confiable}")
