#!/usr/bin/env python3
""" Archivos """

from pathlib import Path

# Todos los archivos de ejemplo se guardan en datos/, junto a este
# script, sin importar desde qué carpeta se ejecute.
CARPETA_DATOS = Path(__file__).resolve().parent / "datos"

###############################################
# Abrir y cerrar archivos: open() / close()
###############################################

# open(ruta, modo) regresa un objeto archivo. Hay que cerrarlo con
# close() cuando ya no se necesita, para que el sistema operativo
# libere el archivo y (en escritura) se guarde todo en disco.
ruta_saludo = CARPETA_DATOS / "saludo.txt"

archivo = open(ruta_saludo, "w")
archivo.write("Hola, física computacional\n")
archivo.close()

print(f"Escribimos {ruta_saludo.name} y lo cerramos a mano.")

# Problema: si algo falla entre el open() y el close() (una excepción,
# un return de por medio...), close() nunca se ejecuta y el archivo se
# queda abierto. Por eso casi nunca se usa open()/close() a mano.

###############################################
# El bloque `with`: manejo seguro de archivos
###############################################

# `with` abre el archivo, lo deja disponible dentro del bloque, y
# garantiza que se cierre al salir del bloque -- incluso si hay una
# excepción en medio. Es la forma recomendada de trabajar con archivos.
with open(ruta_saludo, "r") as archivo:
    contenido = archivo.read()

print(f"Contenido de {ruta_saludo.name}: {contenido!r}")
print(f"¿Está cerrado fuera del bloque `with`? {archivo.closed}")

###############################################
# Modos de apertura
###############################################

# 'r'  -- leer (default). Falla si el archivo no existe.
# 'w'  -- escribir. Crea el archivo si no existe; si ya existe, BORRA
#         su contenido antes de escribir (¡cuidado!).
# 'a'  -- append. Crea el archivo si no existe; si ya existe, escribe
#         al final, sin borrar lo que había.
# 'x'  -- crear. Falla si el archivo ya existe (evita sobreescribir
#         por accidente).
# 'r+' -- leer y escribir, sin borrar el contenido existente.
#
# Además se puede combinar con 't' (texto, default) o 'b' (binario,
# por ejemplo para leer imágenes o datos empaquetados con `struct`).

ruta_bitacora = CARPETA_DATOS / "bitacora.txt"

with open(ruta_bitacora, "w") as archivo:
    archivo.write("corrida 1: convergió en 12 iteraciones\n")

with open(ruta_bitacora, "a") as archivo:
    archivo.write("corrida 2: convergió en 9 iteraciones\n")

with open(ruta_bitacora, "r") as archivo:
    print(f"{ruta_bitacora.name} tiene {len(archivo.readlines())} líneas (append no borra)")

###############################################
# Leer archivos: distintas formas
###############################################

with open(ruta_bitacora, "r") as archivo:
    todo = archivo.read()  # todo el archivo como un solo string
print(f"read() completo:\n{todo}")

with open(ruta_bitacora, "r") as archivo:
    primera_linea = archivo.readline()  # una línea a la vez, con '\n' incluido
    segunda_linea = archivo.readline()
print(f"readline() dos veces: {primera_linea!r}, {segunda_linea!r}")

with open(ruta_bitacora, "r") as archivo:
    lineas = archivo.readlines()  # lista de strings, una por línea
print(f"readlines(): {lineas}")

# La forma más común y eficiente para archivos grandes: iterar
# directamente sobre el objeto archivo, línea por línea, sin cargar
# todo en memoria de golpe.
with open(ruta_bitacora, "r") as archivo:
    for numero_de_linea, linea in enumerate(archivo, start=1):
        print(f"  línea {numero_de_linea}: {linea.strip()}")

###############################################
# Escribir archivos: write() / writelines()
###############################################

ruta_potencias = CARPETA_DATOS / "potencias_de_dos.txt"

with open(ruta_potencias, "w") as archivo:
    for exponente in range(8):
        # write() NO agrega saltos de línea automáticamente, a
        # diferencia de print(); hay que ponerlos explícitamente.
        archivo.write(f"2^{exponente} = {2 ** exponente}\n")

print(f"Escribimos {ruta_potencias.name}")

###############################################
# Física computacional: guardar y leer datos de una "simulación"
###############################################

# Ejemplo: caída libre con paso de tiempo fijo (misma idea que en la
# unidad 01/02), guardando tiempo y posición en un archivo de texto
# con formato de columnas separadas por comas (como un .csv simple).

GRAVEDAD = 9.81
altura_inicial = 20.0
paso_de_tiempo = 0.1

ruta_caida_libre = CARPETA_DATOS / "caida_libre.csv"

with open(ruta_caida_libre, "w") as archivo:
    archivo.write("tiempo,posicion\n")  # encabezado de las columnas
    tiempo = 0.0
    posicion = altura_inicial
    while posicion >= 0:
        archivo.write(f"{tiempo:.2f},{posicion:.4f}\n")
        tiempo += paso_de_tiempo
        posicion = altura_inicial - 0.5 * GRAVEDAD * tiempo**2

print(f"Simulación guardada en {ruta_caida_libre.name}")

# Y ahora lo volvemos a leer, para mostrar cómo se recupera de disco
# lo que guardó otra corrida (o incluso otro script):
with open(ruta_caida_libre, "r") as archivo:
    encabezado = archivo.readline()  # nos la saltamos, ya la conocemos
    datos = []
    for linea in archivo:
        texto_tiempo, texto_posicion = linea.strip().split(",")
        datos.append((float(texto_tiempo), float(texto_posicion)))

print(f"Leímos {len(datos)} puntos de vuelta; los primeros 3: {datos[:3]}")

# Nota: para archivos .csv "de verdad" (con comas dentro de campos de
# texto, comillas, etc.) conviene usar el módulo estándar `csv` en vez
# de partir las líneas a mano con split(); aquí lo hicimos manual
# porque nuestros datos son solo dos números por línea.
