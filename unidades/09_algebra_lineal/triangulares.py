#!/usr/bin/env python3
"""Sistemas triangulares: sustitución hacia adelante y hacia atrás.

Ver notas.md, sección "Matrices triangulares". Las funciones
sustitucion_adelante y sustitucion_atras viven en
fiscomp/algebra_lineal.py, porque las reutilizan la eliminación
gaussiana y la descomposición LU (eliminacion_gaussiana_lu.py).

Este script:

1. Arma la matriz de prueba del libro (código 4.1, testcreate;
   crear_prueba en fiscomp/algebra_lineal.py), toma su parte
   triangular inferior y superior, y resuelve L x = b y U x = b.
2. Comprueba las soluciones sin cajas negras: con el residuo, y
   resolviendo un sistema cuya solución exacta conocemos de antemano.
3. Cuenta las operaciones de punto flotante de la sustitución hacia
   adelante y del producto matriz-vector, y las compara contra n^2 y
   2n^2 - n (sección "Conteo de operaciones").
4. Mide el tiempo de la sustitución hacia adelante al duplicar n.
"""

import time

from fiscomp.matrices import Matrix
from fiscomp.precision_numerica import error_relativo
from fiscomp.algebra_lineal import (
    mat_vec,
    residuo,
    norma_infinito_vec,
    triangular_inferior,
    triangular_superior,
    sustitucion_adelante,
    sustitucion_atras,
    crear_prueba,
)

###############################################
# Utilidades
###############################################


def imprimir_vector(nombre, x):
    print(f"{nombre} = [" + ", ".join(f"{x_i:.10f}" for x_i in x) + "]")


###############################################
# 1 y 2. Resolver y comprobar
###############################################

A, b = crear_prueba(4, 21)
print("Matriz de prueba A:")
print(Matrix(A))
imprimir_vector("b", b)

for nombre, T, resolver in [
    ("inferior (sustitución hacia adelante)", triangular_inferior(A), sustitucion_adelante),
    ("superior (sustitución hacia atrás)", triangular_superior(A), sustitucion_atras),
]:
    print(f"\n--- Triangular {nombre} ---")
    x = resolver(T, b)
    imprimir_vector("x", x)
    print(f"||b - T x||_inf = {norma_infinito_vec(residuo(T, x, b)):.2e}")

    # Sistema con solución conocida: elegimos x_exacta y construimos b = T x_exacta
    x_exacta = [1.0, 2.0, 3.0, 4.0]
    x_calc = resolver(T, mat_vec(T, x_exacta))
    errores = [error_relativo(xc, xe) for xc, xe in zip(x_calc, x_exacta)]
    print(f"Con b = T [1, 2, 3, 4], error relativo máximo = {max(errores):.2e}")

###############################################
# 3. Conteo de operaciones
###############################################


def sustitucion_adelante_contando(L, b):
    """Igual que sustitucion_adelante, pero además cuenta las operaciones
    de punto flotante. Regresa (x, operaciones).
    """
    n = len(b)
    x = [0.0] * n
    operaciones = 0
    for i in range(n):
        suma = 0.0
        for j in range(i):
            suma += L[i][j] * x[j]
            # Una multiplicación y una suma/resta: sumar los i productos
            # toma i - 1 sumas, y restar el resultado de b[i] una más,
            # así que son i sumas/restas en total (ver notas.md).
            operaciones += 2
        x[i] = (b[i] - suma) / L[i][i]
        operaciones += 1  # la división
    return x, operaciones


def mat_vec_contando(A, x):
    """Igual que mat_vec, pero cuenta operaciones. Regresa (y, operaciones)."""
    y = []
    operaciones = 0
    for renglon in A:
        suma = renglon[0] * x[0]
        operaciones += 1
        for j in range(1, len(x)):
            suma += renglon[j] * x[j]
            operaciones += 2
        y.append(suma)
    return y, operaciones


print("\n--- Conteo de operaciones ---")
print(f"{'n':>6} {'sust. adelante':>15} {'n^2':>10} {'mat_vec':>10} {'2n^2 - n':>10}")
for n in [1, 2, 4, 8, 16, 32]:
    A, b = crear_prueba(n, 21)
    _, ops_sust = sustitucion_adelante_contando(triangular_inferior(A), b)
    _, ops_mv = mat_vec_contando(A, b)
    print(f"{n:>6} {ops_sust:>15} {n**2:>10} {ops_mv:>10} {2 * n**2 - n:>10}")

###############################################
# 4. Tiempo de ejecución
###############################################

print("\n--- Tiempo de la sustitución hacia adelante ---")
print("Si el costo es ~n^2, duplicar n debe cuadruplicar el tiempo.")
tiempo_anterior = None
for n in [250, 500, 1000, 2000]:
    # Aquí no usamos crear_prueba: calcular millones de raíces con
    # raiz_cuadrada (desde cero) tardaría más que la sustitución misma.
    # Cualquier L triangular con diagonal distinta de cero sirve.
    L = [[2.0 if i == j else 1.0 / n for j in range(i + 1)] + [0.0] * (n - i - 1)
         for i in range(n)]
    b = [1.0] * n
    tiempo = None
    for repeticion in range(5):  # nos quedamos con la más rápida de 5 corridas
        inicio = time.perf_counter()
        sustitucion_adelante(L, b)
        duracion = time.perf_counter() - inicio
        if tiempo is None or duracion < tiempo:
            tiempo = duracion
    if tiempo_anterior is None:
        print(f"n = {n:4d}: {tiempo:.4f} s")
    else:
        print(f"n = {n:4d}: {tiempo:.4f} s  (x{tiempo / tiempo_anterior:.1f} respecto a n/2)")
    tiempo_anterior = tiempo
