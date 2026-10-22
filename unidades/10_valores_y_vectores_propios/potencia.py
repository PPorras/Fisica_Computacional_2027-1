#!/usr/bin/env python3
"""Métodos de la potencia y de la potencia inversa con desplazamiento.

Ver notas.md, secciones "El problema de eigenvalores, numéricamente",
"Método de la potencia" y "Método de la potencia inversa con
desplazamiento". Las funciones producto_punto, normalizar,
cociente_rayleigh, potencia y potencia_inversa viven en
fiscomp/algebra_lineal.py. Este script:

1. Aplica el método de la potencia a la matriz de prueba del libro,
   iteración por iteración, y muestra que el error del eigenvalor
   baja en cada paso por un factor |lambda_1 / lambda_0|.
2. Con la potencia inversa y varios desplazamientos encuentra los
   cuatro eigenvalores, y los comprueba con la traza y el
   determinante.
3. Muestra que entre más cerca está el desplazamiento del eigenvalor,
   menos iteraciones hacen falta.
4. Muestra por qué potencia_inversa corrige el signo del vector
   cuando el eigenvalor es negativo.
"""

from fiscomp.algebra_lineal import (
    mat_vec,
    normalizar,
    cociente_rayleigh,
    cambio_relativo,
    potencia,
    potencia_inversa,
    determinante,
    norma_euclidea,
    crear_prueba,
)


def formatear(x):
    return "[" + ", ".join(f"{x_i:13.10f}" for x_i in x) + "]"


def residuo_eigen(A, eigenvalor, v):
    """||A v - lambda v|| / ||v||: cero si (lambda, v) es un eigenpar exacto."""
    Av = mat_vec(A, v)
    r = [Av_i - eigenvalor * v_i for Av_i, v_i in zip(Av, v)]
    return norma_euclidea(r) / norma_euclidea(v)


n, val = 4, 21
A, _ = crear_prueba(n, val)

###############################################
# 1. Método de la potencia, iteración por iteración
###############################################

print("--- 1. Método de la potencia con la matriz de prueba del libro ---")
# Referencia: muchas iteraciones, hasta que ya no cambia nada.
lam_ref, _ = potencia(A, kmax=30)

q = normalizar([1.0] * n)
error_anterior = None
print(f"{'k':>3} {'q^(k)':^60} {'mu(q^(k))':>18} {'error rel.':>10} {'cociente':>9}")
for k in range(1, 6):
    z = mat_vec(A, q)
    q = normalizar(z)
    mu = cociente_rayleigh(A, q)
    error = abs(mu - lam_ref) / abs(lam_ref)
    cociente = f"{error / error_anterior:9.2e}" if error_anterior else ""
    print(f"{k:3d} {formatear(q)} {mu:18.13f} {error:10.1e} {cociente}")
    error_anterior = error

lam, q = potencia(A)
print(f"\npotencia(A): lambda_0 = {lam:.10f}")
print(f"             v_0      = {formatear(q)}")
print("Libro:       21.3166626635, [0.44439562 0.48218122 0.51720906 0.55000781]")
print(f"Residuo ||A v - lambda v|| / ||v|| = {residuo_eigen(A, lam, q):.1e}")

###############################################
# 2. Todos los eigenvalores con la potencia inversa
###############################################

print("\n--- 2. Potencia inversa con distintos desplazamientos ---")
print(f"{'desplazamiento s':>17} {'eigenvalor más cercano':>24} {'iteraciones':>12} {'residuo':>9}")
eigenvalores = []
for s in [20.0, 1.0, -0.05, -1e-4, 0.0]:
    lam, v, k = potencia_inversa(A, desplazamiento=s)
    if lam is None:
        print(f"{s:>17g} {'(no convergió)':>24} {k:>12}")
        continue
    eigenvalores.append(lam)
    print(f"{s:>17g} {lam:>24.10e} {k:>12} {residuo_eigen(A, lam, v):>9.1e}")

print("\nCon s = 1 no converge en 200 iteraciones: los dos eigenvalores más")
print("cercanos a 1 (-0.00000008 y -0.00007) están casi a la misma distancia")
print("de 1, y el cociente |lambda*_3 / lambda*_2| es 0.99993.")

traza = sum(A[i][i] for i in range(n))
producto = 1.0
for lam in eigenvalores:
    producto *= lam
print("\nComprobación con los cuatro eigenvalores encontrados:")
print(f"  suma     = {sum(eigenvalores):.10f}   traza(A) = {traza:.10f}")
print(f"  producto = {producto:.6e}   det(A)   = {determinante(A):.6e}")
eigenvalores.sort(key=abs, reverse=True)
print(f"  |lambda_1 / lambda_0| = {abs(eigenvalores[1] / eigenvalores[0]):.2e}"
      " (el cociente de errores de la parte 1)")

###############################################
# 3. Iteraciones contra qué tan bueno es el desplazamiento
###############################################

print("\n--- 3. Buscando lambda_1 = -0.0672...: iteraciones contra desplazamiento ---")
print(f"{'s':>9} {'iteraciones':>12} {'|lambda*_1 / lambda*_2|':>24}")
lam_1, lam_2 = eigenvalores[1], eigenvalores[2]
for s in [-0.5, -0.2, -0.1, -0.07, -0.0672]:
    _, _, k = potencia_inversa(A, desplazamiento=s)
    # En A - s I, el eigenvalor buscado es lam_1 - s y el siguiente
    # más cercano a s es lam_2 - s.
    print(f"{s:>9g} {k:>12} {abs((lam_1 - s) / (lam_2 - s)):>24.3f}")
print("Entre más chico el cociente, más rápido converge.")

###############################################
# 4. Eigenvalor negativo: el signo de q
###############################################

print("\n--- 4. Con -A, el eigenvalor dominante es negativo ---")
menos_A = [[-a for a in renglon] for renglon in A]
q = normalizar([1.0] * n)
for k in range(1, 6):
    q_nuevo = normalizar(mat_vec(menos_A, q))
    print(f"{k:3d} {formatear(q_nuevo)}  cambio relativo = {cambio_relativo(q, q_nuevo):.1e}")
    q = q_nuevo
print("q cambia de signo en cada iteración: el cambio relativo se queda en")
print(f"2n = {2 * n} y un criterio de paro nunca se cumpliría. potencia() no")
print("tiene criterio de paro, así que no le importa:")
lam, _ = potencia(menos_A)
print(f"  potencia(-A): lambda_0 = {lam:.10f}")
lam, _, k = potencia_inversa(menos_A, desplazamiento=-20.0)
print(f"  potencia_inversa(-A, s=-20), que corrige el signo: {lam:.10f} en {k} iteraciones")
