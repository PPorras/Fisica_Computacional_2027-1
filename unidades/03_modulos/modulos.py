#!/usr/bin/env python3
""" Módulos """

###############################################
# Más sobre módulos
###############################################

# Importamos nuestro propio módulo: constantes_fisicas.py está en la
# misma carpeta que este script, así que Python lo encuentra directo.
import constantes_fisicas

print(f"constantes_fisicas.GRAVEDAD = {constantes_fisicas.GRAVEDAD}")
print(f"E_potencial = {constantes_fisicas.energia_potencial(2.0, 10.0)} J")

# También se puede importar solo lo que se necesita, con su propio nombre
from constantes_fisicas import VELOCIDAD_LUZ, energia_potencial

print(f"VELOCIDAD_LUZ = {VELOCIDAD_LUZ} m/s")
print(f"E_potencial en la Luna = {energia_potencial(2.0, 10.0, gravedad=1.62)} J")

### Ejecutar módulos como scripts
# Cada módulo tiene un atributo __name__:
# - si el módulo se ejecuta directamente (`python3 archivo.py`),
#   dentro de ese archivo __name__ vale "__main__"
# - si el módulo se importa desde otro script, __name__ vale el
#   nombre del módulo (el nombre del archivo, sin ".py")
print(f"__name__ de este script: {__name__}")
print(f"__name__ del módulo importado: {constantes_fisicas.__name__}")

# Por eso constantes_fisicas.py tiene un bloque `if __name__ == "__main__":`
# al final: ese código solo corre si ejecutas ese archivo directamente.
# No lo viste ejecutarse arriba, porque aquí solo lo importamos.

### La ruta de búsqueda de módulos (Module Search Path)
import sys

print("Python busca los módulos en estas carpetas (sys.path):")
for carpeta in sys.path:
    print(f"  {carpeta}")

### Archivos "compilados" de Python
# Al importar un módulo, Python guarda su bytecode compilado en una
# carpeta __pycache__/ junto al módulo, para no recompilarlo la próxima
# vez si el archivo no cambió. Por eso el .gitignore del curso ya
# ignora __pycache__/ (ver notas_git_basico.md).
print("Tras correr este script debería aparecer unidades/03_modulos/__pycache__/")

###############################################
# Módulos estándar
###############################################

# sys y math son módulos estándar: vienen con Python, no hay que instalarlos
print(f"sys.argv (argumentos con los que se llamó este script): {sys.argv}")
print(f"Versión de Python: {sys.version.split()[0]}")

import math

print(f"math.pi = {math.pi}")
print(f"math.sqrt(2) = {math.sqrt(2)}")

###############################################
# La función dir()
###############################################

print("Nombres definidos en constantes_fisicas (sin los 'dunder'):")
print([nombre for nombre in dir(constantes_fisicas) if not nombre.startswith("__")])

###############################################
# Paquetes
###############################################

# fiscomp es un paquete: una carpeta (fuera de unidades/) con un
# __init__.py que agrupa varios módulos relacionados (ver fiscomp/__init__.py).
# Como fiscomp vive en la raíz del repositorio y no en esta carpeta,
# hay que agregar esa ruta a sys.path antes de poder importarlo:
from pathlib import Path

raiz_del_repo = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(raiz_del_repo))

import fiscomp
from fiscomp.precision_numerica import EPS

print(f"fiscomp.__doc__: {fiscomp.__doc__.splitlines()[0]}")
print(f"EPS importado desde fiscomp.precision_numerica: {EPS}")

# Cuando funciones_especiales.py necesite EPS más adelante, podrá
# importarlo con una referencia dentro del propio paquete:
#     from .precision_numerica import EPS
# (el "." antes del nombre del módulo indica "en este mismo paquete")

# `from fiscomp import *` importaría todos los nombres públicos de
# fiscomp/__init__.py; como aquí no se define `__all__`, conviene evitarlo
# y seguir importando explícitamente (como arriba) para que quede claro
# de dónde viene cada nombre.
