#!/usr/bin/env python3
"""Análisis de error para sistemas lineales: el ejemplo de Kahan.

Ver notas.md, sección "Análisis de error". Este script reproduce, con
números, lo que ahí se discute:

1. Un residuo diminuto no garantiza una buena solución (a posteriori).
2. Cambiar un solo elemento de A en menos de 0.1% cambia la solución
   por completo (a priori: el problema es muy sensible).
3. El criterio "determinante chico" no sirve para detectar esto.
4. El número de condición kappa(A) = ||A|| ||A^-1|| sí lo detecta, y
   la cota ||dx||/||x|| <= kappa(A) ||dA||/||A|| se cumple.

Este script va antes de la eliminación gaussiana y LU en las notas,
así que todo es 2x2 y se resuelve con fórmulas cerradas (regla de
Cramer). Con LU, la inversa y kappa(A) de una matriz de cualquier
tamaño están en eliminacion_gaussiana_lu.py.
"""

from fiscomp.matrices import Matrix
from fiscomp.algebra_lineal import (
    residuo,
    norma_infinito,
    norma_infinito_vec,
    norma_frobenius,
)

###############################################
# Herramientas para matrices 2x2
###############################################


def determinante_2x2(A):
    """det(A) = A00 A11 - A01 A10."""
    return A[0][0] * A[1][1] - A[0][1] * A[1][0]


def inversa_2x2(A):
    """Inversa de una matriz 2x2:

        A^-1 = (1/det(A)) [[ A11, -A01],
                           [-A10,  A00]]
    """
    det = determinante_2x2(A)
    return [
        [A[1][1] / det, -A[0][1] / det],
        [-A[1][0] / det, A[0][0] / det],
    ]


def resolver_2x2(A, b):
    """Resuelve A x = b (2x2) con la regla de Cramer."""
    det = determinante_2x2(A)
    x0 = (b[0] * A[1][1] - A[0][1] * b[1]) / det
    x1 = (A[0][0] * b[1] - b[0] * A[1][0]) / det
    return [x0, x1]


def numero_condicion(A, norma):
    """kappa(A) = ||A|| ||A^-1|| para una matriz 2x2.

    norma: función que recibe una matriz y regresa su norma (por
    ejemplo norma_infinito o norma_frobenius).
    """
    return norma(A) * norma(inversa_2x2(A))


###############################################
# El sistema de Kahan (ec. 4.15 del libro)
###############################################

A = [[1.2969, 0.8648], [0.2161, 0.1441]]
b = [0.8642, 0.1440]

print("Sistema de Kahan, (A|b):")
print(Matrix([A[0] + [b[0]], A[1] + [b[1]]]))

# 1. A posteriori: el residuo engaña
x_aprox = [0.9911, -0.4870]
x_exacta = [2.0, -2.0]

print("\n--- 1. A posteriori: residuo contra error ---")
print(f"solución aproximada x~ = {x_aprox}")
print(f"residuo r = b - A x~   = {residuo(A, x_aprox, b)}")
print(f"solución exacta x      = {x_exacta}")
print(f"residuo de la exacta   = {residuo(A, x_exacta, b)}")
print("El residuo de x~ es ~1e-8 y sin embargo x~ no tiene ni una cifra correcta.")

# 2. A priori: perturbar un solo elemento
print("\n--- 2. A priori: perturbar un solo elemento de A ---")
A_pert = [[1.2970, 0.8648], [0.2161, 0.1441]]  # A00 cambia en ~0.008 %
x_pert = resolver_2x2(A_pert, b)
print(f"A00: {A[0][0]} -> {A_pert[0][0]}  "
      f"(cambio relativo {abs(A_pert[0][0] - A[0][0]) / A[0][0]:.1e})")
print(f"solución con A original:    {resolver_2x2(A, b)}")
print(f"solución con A perturbada:  {x_pert}")
print("Ojo: ni siquiera con A original sale exactamente [2, -2]: el error de")
print("redondeo en los datos (~1e-16) también se amplifica, por un factor ~kappa (ver 4).")

# 3. El determinante no sirve como criterio
print("\n--- 3. ¿Determinante chico? ---")
print(f"det(A) = {determinante_2x2(A):.3e},  ||A||_inf = {norma_infinito(A):.4f}")
print("Aquí |det(A)| << ||A|| y el problema sí es mal condicionado, pero...")

n = 20
D = [[0.1 if i == j else 0.0 for j in range(n)] for i in range(n)]
det_D = 0.1**n  # D es diagonal: su determinante es el producto de la diagonal
print(f"D = 0.1 * I ({n}x{n}):  det(D) = {det_D:.1e},  ||D||_inf = {norma_infinito(D)}")
print("También |det(D)| << ||D||, y sin embargo D x = b se resuelve sin problema")
print("(x = 10 b): D es perfectamente bien condicionada, kappa(D) = 1.")

# 4. Número de condición y la cota (4.30)
print("\n--- 4. Número de condición ---")
kappa_inf = numero_condicion(A, norma_infinito)
kappa_F = numero_condicion(A, norma_frobenius)
print(f"kappa_inf(A) = {kappa_inf:.3e}")
print(f"kappa_F(A)   = {kappa_F:.3e}")

dA = [[A_pert[i][j] - A[i][j] for j in range(2)] for i in range(2)]
dx = [x_pert[i] - x_exacta[i] for i in range(2)]
cambio_relativo_A = norma_infinito(dA) / norma_infinito(A)
cambio_relativo_x = norma_infinito_vec(dx) / norma_infinito_vec(x_exacta)
print(f"||dA||/||A||                = {cambio_relativo_A:.3e}")
print(f"||dx||/||x||                = {cambio_relativo_x:.3e}")
print(f"cota kappa * ||dA||/||A||   = {kappa_inf * cambio_relativo_A:.3e}")
print("El cambio en x es enorme, pero queda por debajo de la cota, como debe.")
