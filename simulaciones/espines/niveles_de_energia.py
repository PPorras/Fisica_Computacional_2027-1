#!/usr/bin/env python3
"""Niveles de energía de dos, tres y cuatro espines 1/2 (ver README.md).

1. Dos espines: H de 4x4 y sus eigenvalores con QR, comparados con
   los exactos del bloque 2x2.
2. Tres espines: H de 8x8 y sus eigenvalores con QR, con
   comprobaciones (traza, estados con todos los espines alineados) y
   comparados con el barrido de det(H - lambda I).
3. Casos límite (problema 4.48 de [1]): sin interacción y sin campo.
   metodo_qr sin desplazar falla en los dos; eigenvalores_qr no. Sin
   campo, las multiplicidades y el espín total S^2.
4. Por qué no usamos eigenpares (problema 4.49 de [1]).
5. Guarda en datos/ las energías contra omega_I para dos, tres y
   cuatro espines (Fig. 4.1 de [1] y problema 4.50).
6. Diez espines (problema 4.51 de [1]): cuánto tardaría.

Unidades: hbar = 1 y gamma = 1 en las gráficas, así que las energías
están en unidades de gamma (en [1], E/(hbar^2 gamma) contra
omega_I/(hbar gamma)).

[1] A. Gezerlis, Numerical Methods in Physics with Python, vol. 1.
    Cambridge University Press, Cambridge, UK, 2023.
"""

import time
from pathlib import Path

from fiscomp.algebra_lineal import eigenpares, metodo_qr
from espines import dos_espines, tres_espines, hamiltoniano, espin_total_cuadrado, es_simetrica
from eigenvalores import (
    eigenvalores_qr,
    eigenvalores_dos_espines,
    eigenvalores_bloque_2x2,
    eigenvalores_por_barrido,
    multiplicidades,
)

CARPETA_DATOS = Path(__file__).resolve().parent / "datos"


def imprimir_matriz(M):
    for renglon in M.data:
        print("  " + " ".join(f"{x:7.3f}" for x in renglon))


def imprimir_lista(nombre, valores):
    print(f"{nombre}: [" + ", ".join(f"{v:.10f}" for v in valores) + "]")


def imprimir_multiplicidades(nombre, valores):
    grupos = multiplicidades(valores)
    textos = [f"{valor:.6f} ({veces} {'vez' if veces == 1 else 'veces'})" for valor, veces in grupos]
    print(f"{nombre}: " + ", ".join(textos))


###############################################
# 1. Dos espines
###############################################

print("--- 1. Dos espines: omega_I = 1, omega_II = 2, gamma = 0.5 ---")
H = dos_espines(1.0, 2.0, 0.5)
imprimir_matriz(H)
print(f"¿Simétrica? {es_simetrica(H)}")

imprimir_lista("QR         ", eigenvalores_qr(H))
imprimir_lista("Bloque 2x2 ", eigenvalores_dos_espines(H))

###############################################
# 2. Tres espines
###############################################

print("\n--- 2. Tres espines: omega = 1, 2, 3, gamma = 0.5 ---")
H = tres_espines(1.0, 2.0, 3.0, 0.5)
imprimir_matriz(H)
print(f"¿Simétrica? {es_simetrica(H)}")

eigenvalores = eigenvalores_qr(H)
imprimir_lista("QR      ", eigenvalores)
imprimir_lista("Barrido ", eigenvalores_por_barrido(H))

# Comprobaciones sin conocer la respuesta:
# - La suma de los eigenvalores es la traza de H.
# - Con los tres espines arriba (estado 0) o abajo (estado 7), H no
#   mezcla con nada: E = -+(omega_I + omega_II + omega_III)/2 + 3 gamma/4.
traza = sum(H[i][i] for i in range(H.rows))
print(f"Suma de eigenvalores = {sum(eigenvalores):.2e}, traza de H = {traza:.2e}")
print(f"Todos arriba: exacta = {-6.0 / 2 + 3 * 0.5 / 4}, calculada = {eigenvalores[0]:.12f}")
print(f"Todos abajo:  exacta = {6.0 / 2 + 3 * 0.5 / 4}, calculada = {eigenvalores[-1]:.12f}")

###############################################
# 3. Casos límite
###############################################

print("\n--- 3a. Sin interacción (gamma = 0): H es diagonal ---")
H = tres_espines(1.0, 2.0, 3.0, 0.0)
imprimir_lista("Diagonal ordenada ", sorted(H[i][i] for i in range(H.rows)))
try:
    metodo_qr(H)
except ZeroDivisionError:
    print("metodo_qr(H):       ZeroDivisionError")
imprimir_lista("eigenvalores_qr(H)", eigenvalores_qr(H))
print("Con omega_I + omega_II = omega_III, los estados arriba-arriba-abajo y")
print("abajo-abajo-arriba tienen energía 0: H es singular, y Gram-Schmidt divide")
print("entre R_jj = 0. Con el desplazamiento, H - s I ya no es singular.")

print("\n--- 3b. Sin campo (omega = 0): niveles degenerados ---")
H = tres_espines(0.0, 0.0, 0.0, 0.5)
imprimir_lista("metodo_qr(H), ordenados", sorted(metodo_qr(H)))
imprimir_multiplicidades("eigenvalores_qr(H)     ", eigenvalores_qr(H))
print(f"El barrido encuentra {len(eigenvalores_por_barrido(H))} de 8.")
print("Sin desplazar, QR no puede separar +3 gamma/4 de -3 gamma/4 (mismo valor")
print("absoluto) y da números que no son eigenvalores. El barrido no ve raíces de")
print("multiplicidad par. Desplazado, QR da los 8.")

print("\nEl espín total: eigenvalores de S^2 = S(S + 1), tres espines:")
imprimir_multiplicidades("  S^2", eigenvalores_qr(espin_total_cuadrado(3)))
print("  15/4 = 3.75 (S = 3/2) 4 veces: un cuadruplete, M = -3/2, ..., 3/2.")
print("  3/4 = 0.75 (S = 1/2) 4 veces: DOS dobletes, M = -1/2, 1/2 cada uno.")
print("  Sin campo, H = (gamma/2)(S^2 - 9/4): E = 3 gamma/4 para S = 3/2 y")
print("  E = -3 gamma/4 para S = 1/2. El cuadruplete y los dos dobletes.")

print("\nDos espines sin campo (gamma = 0.5):")
imprimir_multiplicidades("  H  ", eigenvalores_qr(dos_espines(0.0, 0.0, 0.5)))
imprimir_multiplicidades("  S^2", eigenvalores_qr(espin_total_cuadrado(2)))
print("  El triplete (S = 1, S^2 = 2) en gamma/4 y el singulete (S = 0) en -3 gamma/4.")

###############################################
# 4. ¿Y los eigenvectores?
###############################################

print("\n--- 4. eigenpares (QR + potencia inversa) con tres espines ---")
H = tres_espines(1.0, 2.0, 3.0, 0.5)
try:
    eigenpares(H)
except ZeroDivisionError:
    print("eigenpares(H): ZeroDivisionError")
print("H no mezcla estados con distinto número de espines hacia arriba (bloques")
print("de 1, 3, 3 y 1), así que cada eigenvector tiene ceros exactos fuera de su")
print("bloque. El criterio de paro de potencia_inversa, cambio_relativo, divide")
print("entre cada componente del vector: con un cero, ZeroDivisionError. Y sin")
print("campo, con eigenvalores repetidos, la potencia inversa no podría dar")
print("varios eigenvectores distintos para un mismo eigenvalor.")

###############################################
# 5. Datos para las gráficas (Fig. 4.1 y cuatro espines)
###############################################


def escribir_renglon(archivo, omega: float, energias):
    archivo.write(f"{omega:.4f} " + " ".join(f"{e:.12f}" for e in energias) + "\n")


CARPETA_DATOS.mkdir(exist_ok=True)
gamma = 1.0

ruta = CARPETA_DATOS / "dos_espines.dat"
with open(ruta, "w") as archivo:
    archivo.write("# dos espines, omega_II = 2 omega_I, gamma = hbar = 1\n")
    # Las columnas 2-5 son de QR, ordenadas de menor a mayor. Las 6-9 son
    # las exactas, sin ordenar: cada columna sigue a un mismo estado (los
    # dos alineados y los dos del bloque 2x2), para dibujar los cruces.
    archivo.write("# omega_I  E_0..E_3 (QR)  E_arriba_arriba  E_bloque_menor  E_bloque_mayor  E_abajo_abajo\n")
    for k in range(101):
        omega = 0.02 * k
        H = dos_espines(omega, 2 * omega, gamma)
        menor, mayor = eigenvalores_bloque_2x2(H[1][1], H[1][2], H[2][2])
        escribir_renglon(archivo, omega, eigenvalores_qr(H) + [H[0][0], menor, mayor, H[3][3]])
print(f"\nEnergías de dos espines guardadas en {ruta.name}")

ruta = CARPETA_DATOS / "tres_espines.dat"
inicio = time.perf_counter()
with open(ruta, "w") as archivo:
    archivo.write("# tres espines, omega_II = 2 omega_I, omega_III = 3 omega_I, gamma = hbar = 1\n")
    archivo.write("# omega_I  E_0 ... E_7 (QR, de menor a mayor)\n")
    for k in range(81):
        omega = 0.025 * k
        H = tres_espines(omega, 2 * omega, 3 * omega, gamma)
        escribir_renglon(archivo, omega, eigenvalores_qr(H))
tiempo_tres = (time.perf_counter() - inicio) / 81
print(f"Energías de tres espines guardadas en {ruta.name} ({tiempo_tres:.3f} s por punto)")

ruta = CARPETA_DATOS / "cuatro_espines.dat"
print("Cuatro espines (16x16), unos segundos...")
inicio = time.perf_counter()
with open(ruta, "w") as archivo:
    archivo.write("# cuatro espines, omega_k = k omega_I (k = 1..4), gamma = hbar = 1\n")
    archivo.write("# omega_I  E_0 ... E_15 (QR, de menor a mayor)\n")
    for k in range(81):
        omega = 0.025 * k
        H = hamiltoniano([omega, 2 * omega, 3 * omega, 4 * omega], gamma)
        escribir_renglon(archivo, omega, eigenvalores_qr(H))
tiempo_cuatro = (time.perf_counter() - inicio) / 81
print(f"Energías de cuatro espines guardadas en {ruta.name} ({tiempo_cuatro:.3f} s por punto)")

###############################################
# 6. Diez espines
###############################################

print("\n--- 6. ¿Y diez espines? ---")
cociente = tiempo_cuatro / tiempo_tres
print(f"De 8x8 a 16x16 el tiempo por punto creció {cociente:.1f} veces (QR es O(n^3) por")
print("iteración: se esperaría 2^3 = 8).")
n = 2**10
estimado = tiempo_cuatro * (n / 16) ** 3
print(f"Con 10 espines H es de {n}x{n}: ({n}/16)^3 = {(n // 16) ** 3} veces más que con cuatro,")
print(f"unos {estimado / 3600:.0f} horas por punto con QR en Python puro (sin contar")
print("construir H, que también tiene productos de matrices de 1024x1024).")
