#!/usr/bin/env python3
"""Sensibilidad de eigenvalores y eigenvectores, con matrices 2x2.

Ver notas.md, secciones "Número de condición para eigenvalores
simples" y "Sensibilidad de los eigenvectores". Aquí se comprueba,
con números:

1. kappa_ev = 1 / |u^T v| (u: eigenvector izquierdo, v: derecho, ambos
   normalizados) predice cuánto se mueve un eigenvalor al perturbar A:
   |d lambda| <= kappa_ev ||dA||.
2. Una matriz simétrica tiene kappa_ev = 1 (eigenvalores bien
   condicionados), pero si dos eigenvalores están muy cerca, sus
   eigenvectores pueden girar mucho con una perturbación diminuta.

Como todavía no tenemos un método general para eigenvalores (eso viene
más adelante en el tema), todo es 2x2 y se resuelve con fórmulas
cerradas: los eigenvalores son las raíces del polinomio característico
lambda^2 - tr(A) lambda + det(A) = 0.
"""

from fiscomp.matrices import Matrix
from fiscomp.funciones_especiales import arcotangente, PI
from fiscomp.algebra_lineal import norma_euclidea

###############################################
# Eigenvalores y eigenvectores de una 2x2
###############################################


def eigenvalores_2x2(A):
    """Raíces de lambda^2 - tr(A) lambda + det(A) = 0, de menor a mayor.

    Supone que son reales (discriminante >= 0).
    """
    traza = A[0][0] + A[1][1]
    det = A[0][0] * A[1][1] - A[0][1] * A[1][0]
    raiz = (traza**2 - 4.0 * det) ** 0.5
    return [(traza - raiz) / 2.0, (traza + raiz) / 2.0]


def eigenvector_2x2(A, lam):
    """Eigenvector derecho normalizado de A para el eigenvalor lam.

    Resuelve (A - lam I) v = 0. Del primer renglón,
    (A00 - lam) v0 + A01 v1 = 0, así que v = (A01, lam - A00) sirve si
    A01 != 0; si no, se usa el segundo renglón. Supone eigenvalores
    distintos.
    """
    if A[0][1] != 0.0:
        v = [A[0][1], lam - A[0][0]]
    elif A[1][0] != 0.0:
        v = [lam - A[1][1], A[1][0]]
    elif lam == A[0][0]:  # A diagonal
        v = [1.0, 0.0]
    else:
        v = [0.0, 1.0]
    norma = norma_euclidea(v)
    return [v[0] / norma, v[1] / norma]


def transpuesta(A):
    return Matrix(A).transpose().data


def condicion_eigenvalor(A, lam):
    """kappa_ev = 1 / |u^T v| para el eigenvalor lam (ec. 4.42).

    Los eigenvectores izquierdos de A son los derechos de A^T.
    """
    v = eigenvector_2x2(A, lam)
    u = eigenvector_2x2(transpuesta(A), lam)
    return 1.0 / abs(u[0] * v[0] + u[1] * v[1])


###############################################
# 1. Eigenvalores: matriz no simétrica contra simétrica
###############################################

eps = 1e-6  # perturbación: se suma al elemento A10
matrices = {
    "no simétrica": [[1.0, 1000.0], [0.0, 2.0]],
    "simétrica": [[2.0, 1.0], [1.0, 2.0]],
}

print("--- 1. Condición de eigenvalores ---")
print(f"Perturbación: A10 -> A10 + {eps}  (||dA|| = {eps})")
for nombre, A in matrices.items():
    A_pert = [[A[0][0], A[0][1]], [A[1][0] + eps, A[1][1]]]
    print(f"\nMatriz {nombre}:")
    print(Matrix(A))
    for lam, lam_pert in zip(eigenvalores_2x2(A), eigenvalores_2x2(A_pert)):
        kappa = condicion_eigenvalor(A, lam)
        print(f"  lambda = {lam}:  kappa_ev = {kappa:8.2f},  "
              f"|d lambda| = {abs(lam_pert - lam):.3e},  "
              f"cota kappa ||dA|| = {kappa * eps:.3e}")

###############################################
# 2. Eigenvectores: eigenvalores muy cercanos
###############################################

print("\n--- 2. Sensibilidad de eigenvectores ---")
print("A = [[1, e], [e, 1 + delta]]: simétrica, kappa_ev = 1 para ambos eigenvalores.")
print("Para e = 0 el eigenvector de lambda = 1 es (1, 0). ¿Cuánto gira con e != 0?")
e = 1e-7
for delta in [1.0, 1e-3, 1e-6, 1e-9]:
    A = [[1.0, e], [e, 1.0 + delta]]
    lam = eigenvalores_2x2(A)[0]
    v = eigenvector_2x2(A, lam)
    if v[0] < 0:  # el signo de un eigenvector es arbitrario
        v = [-v[0], -v[1]]
    angulo = arcotangente(v[1] / v[0]) * 180.0 / PI
    print(f"  e = {e:.0e}, delta = {delta:.0e}:  |d lambda| = {abs(lam - 1.0):.1e},  "
          f"el eigenvector gira {abs(angulo):.2e} grados")
print("Con eigenvalores cercanos (delta chico) el eigenvalor casi no se mueve, pero")
print("el eigenvector gira decenas de grados: el denominador lambda_i - lambda_j de")
print("la ec. (4.49) es muy chico.")
