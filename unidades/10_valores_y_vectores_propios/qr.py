#!/usr/bin/env python3
"""Descomposición QR, iteración simultánea y método QR.

Ver README.md, sección "Método QR" y "Todos los eigenvalores y
eigenvectores". Las funciones transpuesta, mat_mat, descomposicion_qr,
metodo_qr y eigenpares viven en fiscomp/algebra_lineal.py. Este script:

1. Descompone la matriz de prueba del libro en Q R y mide qué tan bien
   se reconstruye A y qué tan ortogonal sale Q.
2. Muestra que la diagonal de U (de LU) no da los eigenvalores: A = L U
   no es una transformación de similitud.
3. Iteración simultánea: sin ortonormalizar, todas las columnas acaban
   en el eigenvector dominante; ortonormalizando con QR en cada paso
   funciona, y da lo mismo que el método QR.
4. El método QR paso a paso: la diagonal de A^(k) converge a los
   eigenvalores y la parte de abajo de la diagonal se va a cero. Con
   una matriz simétrica, A^(k) se vuelve diagonal.
5. Todos los eigenpares con eigenpares (QR + potencia inversa),
   comprobados con el residuo, la traza y el determinante.
"""

from fiscomp.algebra_lineal import (
    mat_vec,
    transpuesta,
    mat_mat,
    normalizar,
    norma_euclidea,
    norma_frobenius,
    descomposicion_lu,
    descomposicion_qr,
    metodo_qr,
    eigenpares,
    determinante,
    crear_prueba,
)


def formatear(x):
    return "[" + ", ".join(f"{x_i:13.10f}" for x_i in x) + "]"


def identidad(n):
    return [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]


def resta(A, B):
    return [[a - b for a, b in zip(renglon_A, renglon_B)] for renglon_A, renglon_B in zip(A, B)]


def parte_de_abajo(A):
    """Norma de Frobenius de lo que está abajo de la diagonal: cero si A es
    triangular superior."""
    n = len(A)
    return norma_frobenius([[A[i][j] if i > j else 0.0 for j in range(n)] for i in range(n)])


def fuera_de_la_diagonal(A):
    """Norma de Frobenius de todo lo que no está en la diagonal: cero si A
    es diagonal."""
    n = len(A)
    return norma_frobenius([[A[i][j] if i != j else 0.0 for j in range(n)] for i in range(n)])


def residuo_eigen(A, eigenvalor, v):
    """||A v - lambda v|| / ||v||: cero si (lambda, v) es un eigenpar exacto."""
    Av = mat_vec(A, v)
    r = [Av_i - eigenvalor * v_i for Av_i, v_i in zip(Av, v)]
    return norma_euclidea(r) / norma_euclidea(v)


n, val = 4, 21
A, _ = crear_prueba(n, val)

###############################################
# 1. Descomposición QR
###############################################

print("--- 1. Descomposición QR con Gram-Schmidt clásico ---")
print(f"{'n':>3} {'||A - Q R||':>12} {'||Q^T Q - I||':>14}")
for m in [4, 6, 8]:
    B, _ = crear_prueba(m, val)
    Q, R = descomposicion_qr(B)
    error_A = norma_frobenius(resta(B, mat_mat(Q, R)))
    error_Q = norma_frobenius(resta(mat_mat(transpuesta(Q), Q), identidad(m)))
    print(f"{m:>3} {error_A:12.1e} {error_Q:14.1e}")
print("Q R reconstruye A casi al épsilon de la máquina, pero Q pierde la")
print("ortogonalidad muy rápido: el redondeo le hace daño a Gram-Schmidt clásico.")

###############################################
# 2. LU no es una transformación de similitud
###############################################

print("\n--- 2. ¿La diagonal de U da los eigenvalores? ---")
L, U = descomposicion_lu(A)
print(f"diagonal de U:  {formatear([U[i][i] for i in range(n)])}")
print(f"eigenvalores:   {formatear(metodo_qr(A))}")
print("No: A = L U no es una transformación de similitud. Lo único que")
print("comparten es el producto (el determinante).")

###############################################
# 3. Iteración simultánea
###############################################

print("\n--- 3a. Iteración simultánea sin ortonormalizar: Z = A^10, columnas normalizadas ---")
Z = identidad(n)
for k in range(10):
    Z = mat_mat(A, Z)
for j, columna in enumerate(transpuesta(Z)):
    print(f"columna {j}: {formatear(normalizar(columna))}")
print("Las cuatro columnas son (casi) el mismo vector: el eigenvector dominante v_0.")

print("\n--- 3b. Iteración simultánea con QR contra el método QR ---")
print("Iteración simultánea: Z = A Q, Z = Q R; luego A^(k) = Q^T A Q.")
print("Método QR:            A^(k-1) = Q R, A^(k) = R Q.")
print(f"{'k':>3} {'diagonal de A^(k) (método QR)':^60} {'dif. diagonales':>16}")
Q_sim = identidad(n)
A_k = [renglon[:] for renglon in A]
for k in range(1, 7):
    # Iteración simultánea
    Z = mat_mat(A, Q_sim)
    Q_sim, _ = descomposicion_qr(Z)
    A_sim = mat_mat(mat_mat(transpuesta(Q_sim), A), Q_sim)
    # Método QR
    Q, R = descomposicion_qr(A_k)
    A_k = mat_mat(R, Q)
    diagonal = [A_k[i][i] for i in range(n)]
    diferencia = max(abs(A_sim[i][i] - A_k[i][i]) for i in range(n))
    print(f"{k:>3} {formatear(diagonal)} {diferencia:16.1e}")
print("Los dos métodos dan la misma A^(k): son el mismo método. (Los elementos")
print("de arriba de la diagonal difieren más, ~1e-4, porque la Q de Gram-Schmidt")
print("clásico no es del todo ortogonal y Q^T ya no es exactamente su inversa.)")

###############################################
# 4. El método QR paso a paso
###############################################

print("\n--- 4. Método QR: A^(k) se vuelve triangular superior ---")
print(f"{'k':>3} {'abajo de la diagonal':>21}")
A_k = [renglon[:] for renglon in A]
for k in range(1, 9):
    Q, R = descomposicion_qr(A_k)
    A_k = mat_mat(R, Q)
    print(f"{k:>3} {parte_de_abajo(A_k):21.1e}")
print("A^(8):")
for renglon in A_k:
    print("  " + "  ".join(f"{a:11.3e}" for a in renglon))
print("Cada iteración divide lo de abajo de la diagonal por ~300: el peor cociente")
print("|lambda_(i+1) / lambda_i| de la matriz de prueba es ~3e-3, igual que en la")
print("potencia.")

print("\nCon la matriz simétrica S = A + A^T, A^(k) debería volverse diagonal:")
S = [[A[i][j] + A[j][i] for j in range(n)] for i in range(n)]
S_k = [renglon[:] for renglon in S]
for k in range(1, 9):
    Q, R = descomposicion_qr(S_k)
    S_k = mat_mat(R, Q)
abajo = parte_de_abajo(S_k)
arriba = (fuera_de_la_diagonal(S_k) ** 2 - abajo**2) ** 0.5
print(f"  después de 8 iteraciones: abajo de la diagonal {abajo:.1e},"
      f" arriba de la diagonal {arriba:.1e}")
print(f"  diagonal: {formatear([S_k[i][i] for i in range(n)])}")
print("Lo de abajo se va a cero, pero lo de arriba se atora en ~1e-7: en")
print("aritmética exacta A^(k) seguiría siendo simétrica, pero la Q de Gram-Schmidt")
print("clásico no es ortogonal y R Q ya no es exactamente una transformación de")
print("similitud ortogonal. Los eigenvalores (la diagonal) sí salen bien.")

###############################################
# 5. Todos los eigenvalores y eigenvectores
###############################################

print("\n--- 5. eigenpares(A): método QR + potencia inversa ---")
eigenvalores_qr = metodo_qr(A)
eigenvalores, eigenvectores = eigenpares(A)
print(f"{'QR':>22} {'potencia inversa':>22} {'residuo':>9}")
for lam_qr, lam, v in zip(eigenvalores_qr, eigenvalores, eigenvectores):
    print(f"{lam_qr:22.15e} {lam:22.15e} {residuo_eigen(A, lam, v):9.1e}")
print("La potencia inversa cambia un poco el eigenvalor más chico.")
print("\nEigenvectores (el primero es el de la potencia, sección anterior):")
for v in eigenvectores:
    print(f"  {formatear(v)}")

traza = sum(A[i][i] for i in range(n))
producto = 1.0
for lam in eigenvalores:
    producto *= lam
print("\nComprobación:")
print(f"  suma     = {sum(eigenvalores):.10f}   traza(A) = {traza:.10f}")
print(f"  producto = {producto:.6e}   det(A)   = {determinante(A):.6e}")
