#!/usr/bin/env python3
"""Método iterativo de Jacobi.

Ver notas.md, sección "Método iterativo de Jacobi". Las funciones
paso_jacobi, cambio_relativo y jacobi viven en
fiscomp/algebra_lineal.py. Este script:

1. Resuelve el sistema de prueba del libro, con 21 sumado a la
   diagonal para que sea diagonalmente dominante, imprimiendo cada
   iteración, y compara contra eliminación gaussiana con pivoteo.
2. Muestra que con la matriz original (no diagonalmente dominante)
   Jacobi no converge.
3. Compara el costo de Jacobi y de la eliminación gaussiana al crecer
   n, con una matriz diagonalmente dominante y rala (tridiagonal).
"""

import time

from fiscomp.precision_numerica import error_relativo
from fiscomp.algebra_lineal import (
    paso_jacobi,
    cambio_relativo,
    jacobi,
    eliminacion_gaussiana_pivoteo,
    crear_prueba,
)


def formatear(x):
    return "[" + ", ".join(f"{x_i:.10f}" for x_i in x) + "]"


def es_diagonal_dominante(A):
    """True si |A_ii| >= sum_{j != i} |A_ij| en todos los renglones."""
    return all(
        abs(renglon[i]) >= sum(abs(a) for j, a in enumerate(renglon) if j != i)
        for i, renglon in enumerate(A)
    )


###############################################
# 1. Jacobi, iteración por iteración
###############################################

n, val = 4, 21
A, b = crear_prueba(n, val)
# Como en el libro: sumamos val a la diagonal para hacerla dominante.
A_dominante = [[A[i][j] + (val if i == j else 0.0) for j in range(n)] for i in range(n)]
print(f"--- 1. Matriz de prueba + {val} I (diagonal dominante: {es_diagonal_dominante(A_dominante)}) ---")

x = [0.0] * n
for k in range(1, 50):
    x_nuevo = paso_jacobi(A_dominante, b, x)
    error = cambio_relativo(x, x_nuevo)
    x = x_nuevo
    if k <= 5 or k % 5 == 0 or error < 1e-6:
        print(f"{k:3d} {formatear(x)}  cambio = {error:.1e}")
    if error < 1e-6:
        break

x_jacobi, iteraciones = jacobi(A_dominante, b)
x_gauss = eliminacion_gaussiana_pivoteo(A_dominante, b)
print(f"\njacobi() convergió en {iteraciones} iteraciones.")
print(f"Jacobi:              {formatear(x_jacobi)}")
print(f"Gauss con pivoteo:   {formatear(x_gauss)}")
diferencia = max(error_relativo(xj, xg) for xj, xg in zip(x_jacobi, x_gauss))
print(f"Diferencia relativa máxima: {diferencia:.1e} (del orden de la tolerancia, 1e-6)")

###############################################
# 2. Sin dominancia diagonal
###############################################

print(f"\n--- 2. Matriz de prueba original (diagonal dominante: {es_diagonal_dominante(A)}) ---")
x_jacobi, iteraciones = jacobi(A, b)
print(f"jacobi() regresa {x_jacobi} tras {iteraciones} iteraciones: no convergió.")
x = [0.0] * n
for k in range(1, 6):
    x = paso_jacobi(A, b, x)
    print(f"{k:3d} {formatear(x)}")
print("Las componentes crecen en cada iteración en vez de estabilizarse.")

###############################################
# 3. Costo: Jacobi contra eliminación gaussiana
###############################################


def tridiagonal(n):
    """Matriz n x n con 4 en la diagonal y -1 en las dos vecinas
    (diagonalmente dominante y rala: casi todos sus elementos son 0).
    Aparece al discretizar ecuaciones diferenciales (capítulo 8)."""
    return [[4.0 if i == j else (-1.0 if abs(i - j) == 1 else 0.0) for j in range(n)] for i in range(n)]


print("\n--- 3. Tiempo: Jacobi contra eliminación gaussiana (matriz tridiagonal) ---")
print(f"{'n':>5} {'iteraciones':>12} {'t_jacobi (s)':>13} {'t_gauss (s)':>12}")
for n in [50, 100, 200, 400]:
    A = tridiagonal(n)
    b = [1.0] * n
    inicio = time.perf_counter()
    _, iteraciones = jacobi(A, b, kmax=500, tol=1e-10)
    t_jacobi = time.perf_counter() - inicio
    inicio = time.perf_counter()
    eliminacion_gaussiana_pivoteo(A, b)
    t_gauss = time.perf_counter() - inicio
    print(f"{n:>5} {iteraciones:>12} {t_jacobi:>13.4f} {t_gauss:>12.4f}")
print("Cada iteración de Jacobi cuesta ~2n^2 y aquí el número de iteraciones")
print("casi no cambia con n; la eliminación cuesta ~2n^3/3. Para n grande")
print("gana Jacobi. (Y si no se guardaran los ceros, cada iteración costaría")
print("solo ~5n.)")
