#!/usr/bin/env python3
"""Eliminación gaussiana y descomposición LU.

Ver notas.md, secciones "Eliminación gaussiana" y "Descomposición LU".
Las funciones eliminacion_gaussiana, descomposicion_lu y resolver_lu
viven en fiscomp/algebra_lineal.py. Este script:

1. Resuelve paso a paso el ejemplo 3x3 de las notas, imprimiendo la
   matriz aumentada después de cada operación de renglón.
2. Saca L y U del mismo ejemplo y comprueba que L U = A.
3. Resuelve el sistema de prueba del libro (crear_prueba) con los dos
   métodos y compara: residuo y sistema con solución conocida.
4. Calcula la inversa de A con LU (n sistemas con la misma A) y con
   ella el número de condición kappa(A): explica por qué en el punto 3
   se pierden dígitos.
5. Cuenta las operaciones de la eliminación y de la descomposición LU
   y las compara contra las fórmulas de las notas.
6. Mide el tiempo de resolver muchos sistemas con la misma A: con
   eliminación gaussiana cada vez, contra una sola LU.
7. Muestra que, sin pivoteo, un cero en la diagonal rompe el método.
"""

import time

from fiscomp.matrices import Matrix
from fiscomp.precision_numerica import error_relativo
from fiscomp.algebra_lineal import (
    mat_vec,
    residuo,
    norma_infinito,
    norma_infinito_vec,
    sustitucion_atras,
    eliminacion_gaussiana,
    descomposicion_lu,
    resolver_lu,
    crear_prueba,
)

###############################################
# Utilidades
###############################################


def imprimir_vector(nombre, x):
    print(f"{nombre} = [" + ", ".join(f"{x_i:.10f}" for x_i in x) + "]")


def imprimir_aumentada(A, b):
    """Imprime la matriz aumentada (A|b), un renglón por línea."""
    for renglon, b_i in zip(A, b):
        print("  " + " ".join(f"{a:7.2f}" for a in renglon) + f" | {b_i:7.2f}")


###############################################
# 1. El ejemplo 3x3, paso a paso
###############################################

A_ejemplo = [[2.0, 1.0, 1.0], [1.0, 1.0, -2.0], [5.0, 10.0, 5.0]]
b_ejemplo = [8.0, -2.0, 10.0]

print("--- 1. Eliminación gaussiana, paso a paso ---")
print("Matriz aumentada inicial:")
imprimir_aumentada(A_ejemplo, b_ejemplo)

# Es el mismo ciclo que eliminacion_gaussiana, con un print en medio.
A = [renglon[:] for renglon in A_ejemplo]
b = b_ejemplo[:]
n = len(b)
for j in range(n - 1):
    for i in range(j + 1, n):
        coeficiente = A[i][j] / A[j][j]
        for k in range(j, n):
            A[i][k] -= coeficiente * A[j][k]
        b[i] -= coeficiente * b[j]
        print(f"\nPivote j = {j}: renglón {i} <- renglón {i} - {coeficiente} * renglón {j}")
        imprimir_aumentada(A, b)

x = sustitucion_atras(A, b)
imprimir_vector("\nSustitución hacia atrás: x", x)
imprimir_vector("Con eliminacion_gaussiana:  x", eliminacion_gaussiana(A_ejemplo, b_ejemplo))

###############################################
# 2. L y U del mismo ejemplo
###############################################

print("\n--- 2. Descomposición LU del ejemplo ---")
L, U = descomposicion_lu(A_ejemplo)
print("L (los coeficientes de la eliminación):")
print(Matrix(L))
print("U (lo que queda de A al final de la eliminación):")
print(Matrix(U))
print("L U (debe ser igual a A):")
print(Matrix(L) * Matrix(U))
imprimir_vector("Con resolver_lu: x", resolver_lu(L, U, b_ejemplo))

###############################################
# 3. El sistema de prueba del libro
###############################################

print("\n--- 3. Sistema de prueba del libro (n = 4) ---")
A, b = crear_prueba(4, 21)
x_gauss = eliminacion_gaussiana(A, b)
L, U = descomposicion_lu(A)
x_lu = resolver_lu(L, U, b)
imprimir_vector("Eliminación gaussiana: x", x_gauss)
imprimir_vector("LU:                    x", x_lu)
diferencia = max(error_relativo(xl, xg) for xl, xg in zip(x_lu, x_gauss))
print(f"Diferencia relativa máxima entre las dos: {diferencia:.2e}")
print(f"||b - A x||_inf = {norma_infinito_vec(residuo(A, x_gauss, b)):.2e}")

x_exacta = [1.0, 2.0, 3.0, 4.0]
x_calc = resolver_lu(L, U, mat_vec(A, x_exacta))
errores = [error_relativo(xc, xe) for xc, xe in zip(x_calc, x_exacta)]
imprimir_vector("Con b = A [1, 2, 3, 4]: x", x_calc)
print(f"Error relativo máximo = {max(errores):.2e}")

###############################################
# 4. La inversa con LU y el número de condición
###############################################

print("\n--- 4. Inversa y número de condición ---")
# La columna k de A^-1 resuelve A x = e_k (e_k: columna k de la
# identidad). Son n sistemas con la misma A: una sola LU y n pares de
# sustituciones.
n = len(A)
columnas = []
for k in range(n):
    e_k = [1.0 if i == k else 0.0 for i in range(n)]
    columnas.append(resolver_lu(L, U, e_k))
A_inv = Matrix(columnas).transpose().data

print("A^-1 A (debe ser la identidad):")
for renglon in (Matrix(A_inv) * Matrix(A)).data:
    print("  " + " ".join(f"{a:10.2e}" for a in renglon))
kappa = norma_infinito(A) * norma_infinito(A_inv)
print(f"kappa(A) = ||A||_inf ||A^-1||_inf = {kappa:.2e}")
print("Se pueden perder hasta log10(kappa) dígitos de los ~16 de un float:")
print("es una propiedad de A, no un defecto del método.")

###############################################
# 5. Conteo de operaciones
###############################################


def eliminacion_contando(A, b):
    """Igual que la fase de eliminación de eliminacion_gaussiana, pero
    cuenta operaciones. Regresa el número de operaciones."""
    A = [renglon[:] for renglon in A]
    b = b[:]
    n = len(b)
    operaciones = 0
    for j in range(n - 1):
        for i in range(j + 1, n):
            coeficiente = A[i][j] / A[j][j]
            operaciones += 1  # la división
            for k in range(j, n):
                A[i][k] -= coeficiente * A[j][k]
                operaciones += 2  # multiplicación y resta
            b[i] -= coeficiente * b[j]
            operaciones += 2
    return operaciones


def lu_contando(A):
    """Igual que descomposicion_lu, pero cuenta operaciones."""
    U = [renglon[:] for renglon in A]
    n = len(A)
    operaciones = 0
    for j in range(n - 1):
        for i in range(j + 1, n):
            coeficiente = U[i][j] / U[j][j]
            operaciones += 1
            for k in range(j, n):
                U[i][k] -= coeficiente * U[j][k]
                operaciones += 2
    return operaciones


print("\n--- 5. Conteo de operaciones (sin las sustituciones) ---")
print(f"{'n':>4} {'eliminación':>12} {'fórmula':>10} {'LU':>10} {'fórmula':>10} {'2n^3/3':>10}")
for n in [2, 4, 8, 16, 32, 64]:
    A, b = crear_prueba(n, 21)
    formula_gauss = (4 * n**3 + 9 * n**2 - 13 * n) // 6  # 2n^3/3 + 3n^2/2 - 13n/6
    formula_lu = (4 * n**3 + 3 * n**2 - 7 * n) // 6  # 2n^3/3 + n^2/2 - 7n/6
    print(
        f"{n:>4} {eliminacion_contando(A, b):>12} {formula_gauss:>10} "
        f"{lu_contando(A):>10} {formula_lu:>10} {2 * n**3 / 3:>10.0f}"
    )

###############################################
# 6. Muchos sistemas con la misma A
###############################################

print("\n--- 6. Resolver m sistemas con la misma A (n = 100) ---")
n, m = 100, 20
# Matriz con diagonal dominante: no necesita pivoteo y es rápida de
# construir (crear_prueba calcula n^2 raíces con raiz_cuadrada).
A = [[float(n) if i == j else 1.0 / (i + j + 1) for j in range(n)] for i in range(n)]
lados_derechos = [[float(k + i) for i in range(n)] for k in range(m)]

inicio = time.perf_counter()
for b in lados_derechos:
    eliminacion_gaussiana(A, b)
t_gauss = time.perf_counter() - inicio

inicio = time.perf_counter()
L, U = descomposicion_lu(A)
for b in lados_derechos:
    resolver_lu(L, U, b)
t_lu = time.perf_counter() - inicio

print(f"{m} eliminaciones gaussianas:    {t_gauss:.3f} s")
print(f"1 LU + {m} pares de sustituciones: {t_lu:.3f} s  ({t_gauss / t_lu:.1f} veces más rápido)")
print(f"Predicción con flops: m (2n^3/3) / (2n^3/3 + 2 m n^2) = "
      f"{m * (2 * n**3 / 3) / (2 * n**3 / 3 + 2 * m * n**2):.1f}")

###############################################
# 7. Un pivote cero
###############################################

print("\n--- 7. Un cero en la diagonal ---")
A = [[0.0, 1.0], [1.0, 1.0]]
b = [1.0, 2.0]
print("A = [[0, 1], [1, 1]], b = [1, 2]; la solución es x = [1, 1] y det(A) = -1.")
try:
    eliminacion_gaussiana(A, b)
except ZeroDivisionError as error:
    print(f"eliminacion_gaussiana falla: ZeroDivisionError ({error}).")
    print("La matriz no es singular: solo hay que intercambiar los renglones (pivoteo).")
