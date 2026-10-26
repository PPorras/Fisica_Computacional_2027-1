#!/usr/bin/env python3
"""Descomposición en valores singulares (SVD).

Ver notas.md, sección "Descomposición en valores singulares". Las
funciones matriz_h y svd viven en fiscomp/algebra_lineal.py. Este
script usa la matriz de prueba del libro menos la identidad,
A = crear_prueba(4, 21) - I (como el libro):

1. Los eigenvalores de H vienen en pares +-sigma: el método QR no
   converge con H, pero sí con H - s I.
2. La SVD de A con svd: reconstrucción U S V^T, ortogonalidad de U y
   V, y comprobación con la norma de Frobenius y el determinante.
3. La norma 2 (espectral) de A es sigma_0: ningún vector x estira más
   que v_0. El número de condición en norma 2 es sigma_0 / sigma_(n-1).
4. Resolver A x = b con la SVD: en las bases de V y U el sistema es
   diagonal.
5. Cuando los valores singulares abarcan muchos órdenes de magnitud
   (la matriz de prueba sin restarle I), svd falla: es culpa de
   metodo_qr, no de la idea de usar H.
"""

import random

from fiscomp.algebra_lineal import (
    mat_vec,
    transpuesta,
    mat_mat,
    producto_punto,
    norma_euclidea,
    norma_frobenius,
    descomposicion_lu,
    resolver_lu,
    metodo_qr,
    determinante,
    matriz_h,
    svd,
    crear_prueba,
)


def formatear(x):
    return "[" + ", ".join(f"{x_i:13.10f}" for x_i in x) + "]"


def identidad(n):
    return [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]


def resta(A, B):
    return [[a - b for a, b in zip(renglon_A, renglon_B)] for renglon_A, renglon_B in zip(A, B)]


def producto_svd(U, sigmas, V_T):
    """U S V^T, con S la matriz diagonal de los sigmas."""
    US = [[U_ij * sigma for U_ij, sigma in zip(renglon, sigmas)] for renglon in U]
    return mat_mat(US, V_T)


n, val = 4, 21
A, b = crear_prueba(n, val)
A = resta(A, identidad(n))

###############################################
# 1. Los eigenvalores de H
###############################################

print("--- 1. Eigenvalores de H = [[0, A^T], [A, 0]] con el método QR ---")
H = matriz_h(A)
print(f"H es de {len(H)} x {len(H)} y simétrica: {H == transpuesta(H)}")
print("Sin desplazamiento (diagonal de H^(k)):")
print("  " + formatear(metodo_qr(H)))
desplazamiento = 1.0
H_desplazada = resta(H, [[desplazamiento * I_ij for I_ij in renglon] for renglon in identidad(2 * n)])
eigenvalores_h = [lam + desplazamiento for lam in metodo_qr(H_desplazada)]
print(f"Con desplazamiento s = {desplazamiento} (diagonal de (H - s I)^(k), más s):")
print("  " + formatear(eigenvalores_h))
print("Sin desplazamiento, cada par +sigma, -sigma tiene el mismo valor absoluto y")
print("QR no los separa: la diagonal se queda en ceros. Con el desplazamiento, los")
print("positivos salen bien (los negativos, que se descartan, convergen mucho más")
print("lento porque quedan muy juntos: -sigma - 1 ~ -2 para los tres más chicos).")

###############################################
# 2. La SVD de A
###############################################

print("\n--- 2. svd(A): A = U S V^T ---")
U, sigmas, V_T = svd(A)
print(f"valores singulares: {formatear(sigmas)}")
print(f"||A - U S V^T||  = {norma_frobenius(resta(A, producto_svd(U, sigmas, V_T))):.1e}")
print(f"||U^T U - I||    = {norma_frobenius(resta(mat_mat(transpuesta(U), U), identidad(n))):.1e}")
print(f"||V^T V - I||    = {norma_frobenius(resta(mat_mat(V_T, transpuesta(V_T)), identidad(n))):.1e}")
print("El error ~1e-10 no es el épsilon de la máquina: sigma_2 y sigma_3 están a")
print("~3e-5 de distancia, y un eigenvector de H con un vecino tan cercano es")
print("sensible (sección \"Sensibilidad de los eigenvectores\"): eps ||H|| / 3e-5 ~ 1e-10.")

print("\nA v_i = sigma_i u_i (las columnas de V y U):")
columnas_u = transpuesta(U)
for i in range(n):
    Av = mat_vec(A, V_T[i])  # v_i, la columna i de V, es el renglón i de V^T
    r = [Av_j - sigmas[i] * u_j for Av_j, u_j in zip(Av, columnas_u[i])]
    print(f"  i = {i}: ||A v_i - sigma_i u_i|| = {norma_euclidea(r):.1e}")

suma_cuadrados = sum(sigma**2 for sigma in sigmas)
producto = 1.0
for sigma in sigmas:
    producto *= sigma
print("\nComprobación:")
print(f"  suma de sigma^2   = {suma_cuadrados:.10f}   ||A||_F^2 = {norma_frobenius(A) ** 2:.10f}")
print(f"  producto de sigma = {producto:.10f}   |det(A)|  = {abs(determinante(A)):.10f}")

###############################################
# 3. Norma 2 y número de condición
###############################################

print("\n--- 3. La norma 2 de A es sigma_0 ---")
random.seed(2027)
mayor = 0.0
for _ in range(10000):
    x = [random.uniform(-1.0, 1.0) for _ in range(n)]
    mayor = max(mayor, norma_euclidea(mat_vec(A, x)) / norma_euclidea(x))
v_0 = V_T[0]
print(f"  max ||A x|| / ||x|| con 10000 x al azar: {mayor:.10f}")
print(f"  ||A v_0|| (||v_0|| = 1):                 {norma_euclidea(mat_vec(A, v_0)):.10f}")
print(f"  sigma_0:                                 {sigmas[0]:.10f}")
print(f"Número de condición: kappa_2(A) = sigma_0 / sigma_(n-1) = {sigmas[0] / sigmas[-1]:.4f}")

###############################################
# 4. Resolver A x = b con la SVD
###############################################

print("\n--- 4. A x = b con la SVD: S x' = b', con x' = V^T x, b' = U^T b ---")
b_prima = mat_vec(transpuesta(U), b)
x_prima = [b_i / sigma for b_i, sigma in zip(b_prima, sigmas)]
x_svd = mat_vec(transpuesta(V_T), x_prima)
L, U_lu = descomposicion_lu(A)
x_lu = resolver_lu(L, U_lu, b)
print(f"  x (SVD): {formatear(x_svd)}")
print(f"  x (LU):  {formatear(x_lu)}")
print(f"  x' = V^T x: {formatear(x_prima)}")
print("En las bases de V y U el sistema es diagonal: cada x'_i = b'_i / sigma_i.")

###############################################
# 5. Cuando svd falla
###############################################

print("\n--- 5. La matriz de prueba sin restarle I ---")
B, _ = crear_prueba(n, val)
U_B, sigmas_B, V_T_B = svd(B)
producto = 1.0
for sigma in sigmas_B:
    producto *= sigma
print(f"valores singulares: {[f'{sigma:.6e}' for sigma in sigmas_B]}")
print(f"||B - U S V^T|| = {norma_frobenius(resta(B, producto_svd(U_B, sigmas_B, V_T_B))):.1e}")
print(f"producto de sigma = {producto:.6e}   |det(B)| = {abs(determinante(B)):.6e}")
print(f"v_2 . v_3 = {producto_punto(V_T_B[2], V_T_B[3]):.6f}  (debería ser 0)")
print("Los valores singulares van de ~21 a ~1e-7. Con el desplazamiento s = 1, los")
print("dos más chicos quedan en ~-1 y casi empatan en valor absoluto con los")
print("negativos -sigma - 1 ~ -1: metodo_qr no los separa en 100 iteraciones, y la")
print("potencia inversa encuentra dos veces el mismo vector. La idea de usar H está")
print("bien; lo que falla es el método para sus eigenvalores.")
