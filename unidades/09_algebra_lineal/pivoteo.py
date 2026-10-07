#!/usr/bin/env python3
"""Pivoteo: inestabilidad sin mal condicionamiento.

Ver notas.md, sección "Pivoteo". Las funciones eliminacion_gaussiana y
eliminacion_gaussiana_pivoteo viven en fiscomp/algebra_lineal.py.
Este script:

1. Resuelve los tres ejemplos del libro (todos bien condicionados) sin
   pivoteo y con pivoteo parcial: un cero en el primer pivote, un cero
   "escondido" que aparece a media eliminación, y un pivote diminuto.
2. Compara los dos métodos con la matriz de prueba del libro.
3. Muestra un caso donde el pivoteo parcial no basta (el tercer
   ejemplo con un renglón multiplicado por 10^-20).
4. Muestra que con una matriz diagonalmente dominante el pivoteo no
   intercambia nada.
"""

from fiscomp.precision_numerica import error_relativo
from fiscomp.algebra_lineal import (
    mat_vec,
    eliminacion_gaussiana,
    eliminacion_gaussiana_pivoteo,
    crear_prueba,
)


def formatear(x):
    return "[" + ", ".join(f"{x_i:.10g}" for x_i in x) + "]"


def resolver_con_ambos(A, b):
    """Imprime la solución sin pivoteo (o el error que levanta) y con
    pivoteo parcial."""
    try:
        print(f"  sin pivoteo: {formatear(eliminacion_gaussiana(A, b))}")
    except ZeroDivisionError:
        print("  sin pivoteo: ZeroDivisionError (un pivote vale cero)")
    print(f"  con pivoteo: {formatear(eliminacion_gaussiana_pivoteo(A, b))}")


###############################################
# 1. Los tres ejemplos del libro
###############################################

print("--- 1a. Un cero en el primer pivote (solución exacta: [3, -1]) ---")
A = [[0.0, -1.0], [1.0, 1.0]]
b = [1.0, 2.0]
# La inversa con LU tampoco funciona sin pivoteo: kappa se calcula a mano.
# A^-1 = [[1, 1], [-1, 0]], así que kappa = ||A|| ||A^-1|| = 2 * 2 = 4.
print("  kappa_inf(A) = 4: un problema perfectamente bien condicionado.")
resolver_con_ambos(A, b)

print("\n--- 1b. Un cero escondido (solución exacta: [4, -2, 2]) ---")
A = [[2.0, 1.0, 1.0], [2.0, 1.0, -4.0], [5.0, 10.0, 5.0]]
b = [8.0, -2.0, 10.0]
print("  El primer pivote es 2, pero tras eliminar la columna 0 el")
print("  elemento A_11 queda en 1 - 1*1 = 0.")
resolver_con_ambos(A, b)

print("\n--- 1c. Un pivote diminuto (solución exacta: ~[3, -1]) ---")
A = [[1e-20, -1.0], [1.0, 1.0]]
b = [1.0, 2.0]
# A^-1 = [[1, 1], [-1, 1e-20]] / (1 + 1e-20): kappa = 2 * 2 / (1 + 1e-20) ~ 4.
# (inversa() da otra cosa: usa LU sin pivoteo y sufre el mismo problema.)
print("  kappa_inf(A) ~ 4: también está perfectamente bien condicionado.")
print("  Sin pivoteo, el coeficiente es 1/1e-20 = 1e20, y en punto")
print(f"  flotante 1 + 1e20 = {1.0 + 1e20:.0e} y 2 - 1e20 = {2.0 - 1e20:.0e}:")
print("  se pierde toda la información del segundo renglón.")
resolver_con_ambos(A, b)

###############################################
# 2. La matriz de prueba del libro
###############################################

print("\n--- 2. Matriz de prueba del libro, con solución conocida [1, 2, 3, 4] ---")
A, _ = crear_prueba(4, 21)
x_exacta = [1.0, 2.0, 3.0, 4.0]
b = mat_vec(A, x_exacta)
for nombre, metodo in [("sin pivoteo", eliminacion_gaussiana), ("con pivoteo", eliminacion_gaussiana_pivoteo)]:
    x = metodo(A, b)
    error = max(error_relativo(xc, xe) for xc, xe in zip(x, x_exacta))
    print(f"  {nombre}: error relativo máximo = {error:.2e}")
print("  Los dos pierden ~8 dígitos: con kappa ~ 4e8, eso es culpa de la")
print("  matriz, no del método, y el pivoteo no lo puede arreglar. (La")
print("  diferencia entre los dos es solo el redondeo, que cae distinto.)")

###############################################
# 3. Cuando el pivoteo parcial no basta
###############################################

print("\n--- 3. El ejemplo 1c con el segundo renglón multiplicado por 1e-20 ---")
A = [[1e-20, -1.0], [1e-20, 1e-20]]
b = [1.0, 2e-20]
print("  Son las mismas ecuaciones, así que la solución sigue siendo ~[3, -1].")
print("  Pero ahora los dos candidatos a pivote miden 1e-20: no hay")
print("  intercambio, y el pivoteo parcial falla igual que sin pivoteo.")
resolver_con_ambos(A, b)
print("  El pivoteo parcial escalado (notas.md) compara cada candidato")
print("  contra el elemento más grande de su renglón, y sí intercambia.")

###############################################
# 4. Diagonal dominante: no hace falta pivotear
###############################################

print("\n--- 4. Matriz diagonalmente dominante (solución exacta: [4, -2, 2]) ---")
A = [[2.0, 1.0, 1.0], [5.0, 10.0, 5.0], [1.0, 1.0, -2.0]]
b = [8.0, 10.0, -2.0]
x_sin = eliminacion_gaussiana(A, b)
x_con = eliminacion_gaussiana_pivoteo(A, b)
print(f"  sin pivoteo: {formatear(x_sin)}")
print(f"  con pivoteo: {formatear(x_con)}")
print(f"  ¿Idénticas bit a bit? {x_sin == x_con}: el pivoteo no intercambió nada.")
