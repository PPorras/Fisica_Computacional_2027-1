#!/usr/bin/env python3
"""Niveles de energía de dos y tres espines 1/2 (ver README.md).

1. Construye H para dos espines, comprueba que es real y simétrica, y
   calcula sus 4 eigenvalores de dos formas: con el bloque 2x2 y con
   el barrido de det(H - lambda I).
2. Lo mismo para tres espines (8x8), donde ya no hay atajo: barrido,
   refinamiento e iteración inversa, con comprobaciones (residuo,
   traza, estados con todos los espines alineados).
3. Dos casos límite: sin interacción (H diagonal) y sin campo
   (niveles degenerados, que el barrido no puede ver).
4. Guarda en datos/ las energías contra omega_I, para las dos gráficas
   de la Fig. 4.1 del libro.

Unidades: hbar = 1 y gamma = 1, así que las energías están en unidades
de gamma (en el libro, E/(hbar^2 gamma) contra omega_I/(hbar gamma)).
"""

from pathlib import Path

from espines import dos_espines, tres_espines, es_simetrica
from eigenvalores import (
    eigenvalores_dos_espines,
    eigenvalores_bloque_2x2,
    eigenvalores_por_barrido,
    iteracion_inversa,
    residuo_eigen,
)

CARPETA_DATOS = Path(__file__).resolve().parent / "datos"


def imprimir_matriz(M):
    for renglon in M.data:
        print("  " + " ".join(f"{x:7.3f}" for x in renglon))


def imprimir_lista(nombre, valores):
    print(f"{nombre}: [" + ", ".join(f"{v:.10f}" for v in valores) + "]")


###############################################
# 1. Dos espines
###############################################

print("--- 1. Dos espines: omega_I = 1, omega_II = 2, gamma = 0.5 ---")
H = dos_espines(1.0, 2.0, 0.5)
imprimir_matriz(H)
print(f"¿Simétrica? {es_simetrica(H)}")

imprimir_lista("Bloque 2x2", eigenvalores_dos_espines(H.data))
imprimir_lista("Barrido   ", eigenvalores_por_barrido(H.data))

###############################################
# 2. Tres espines
###############################################

print("\n--- 2. Tres espines: omega = 1, 2, 3, gamma = 0.5 ---")
H = tres_espines(1.0, 2.0, 3.0, 0.5)
imprimir_matriz(H)
print(f"¿Simétrica? {es_simetrica(H)}")

aproximados = eigenvalores_por_barrido(H.data)
print(f"\nEl barrido encontró {len(aproximados)} eigenvalores:")
imprimir_lista("Barrido           ", aproximados)

eigenvalores = []
residuos = []
for sigma in aproximados:
    eigenvalor, v = iteracion_inversa(H.data, sigma)
    eigenvalores.append(eigenvalor)
    residuos.append(residuo_eigen(H.data, eigenvalor, v))
imprimir_lista("Iteración inversa ", eigenvalores)
print(f"Residuo máximo ||H v - lambda v|| / ||v|| = {max(residuos):.1e}")

# Comprobaciones sin conocer la respuesta:
# - La suma de los eigenvalores es la traza de H.
# - Con los tres espines arriba (estado 0) o abajo (estado 7), H no
#   mezcla con nada: E = -+(omega_I + omega_II + omega_III)/2 + 3 gamma/4.
traza = sum(H.data[i][i] for i in range(H.rows))
print(f"Suma de eigenvalores = {sum(eigenvalores):.2e}, traza de H = {traza:.2e}")
print(f"Todos arriba: exacta = {-6.0 / 2 + 3 * 0.5 / 4}, calculada = {eigenvalores[0]:.12f}")
print(f"Todos abajo:  exacta = {6.0 / 2 + 3 * 0.5 / 4}, calculada = {eigenvalores[-1]:.12f}")

###############################################
# 3. Casos límite
###############################################

print("\n--- 3a. Sin interacción (gamma = 0): H es diagonal ---")
H = tres_espines(1.0, 2.0, 3.0, 0.0)
diagonal = sorted(H.data[i][i] for i in range(H.rows))
imprimir_lista("Diagonal", diagonal)
imprimir_lista("Barrido ", eigenvalores_por_barrido(H.data))
print("Una matriz diagonal es triangular: sus eigenvalores son su diagonal.")
print("El barrido no ve el 0, que está dos veces (multiplicidad par; ver 3b).")

print("\n--- 3b. Sin campo (omega = 0): niveles degenerados ---")
H = tres_espines(0.0, 0.0, 0.0, 0.5)
encontrados = eigenvalores_por_barrido(H.data)
print(f"El barrido encontró {len(encontrados)} de 8 eigenvalores: {[round(e, 6) for e in encontrados]}")
print("Los 8 son 3 gamma/4 = 0.375 (4 veces) y -3 gamma/4 = -0.375 (4 veces).")
print("Una raíz de multiplicidad par no cambia el signo de det(H - lambda I):")
print("el barrido no la ve. Para esto hace falta otro método (QR).")

###############################################
# 4. Datos para las gráficas (Fig. 4.1)
###############################################


def escribir_renglon(archivo, omega: float, energias, total: int):
    """omega y las energías ordenadas; '?' (dato faltante para gnuplot)
    si se encontraron menos de `total`."""
    columnas = [f"{e:.12f}" for e in sorted(energias)] + ["?"] * (total - len(energias))
    archivo.write(f"{omega:.4f} " + " ".join(columnas) + "\n")


CARPETA_DATOS.mkdir(exist_ok=True)
gamma = 1.0

ruta = CARPETA_DATOS / "dos_espines.dat"
with open(ruta, "w") as archivo:
    archivo.write("# dos espines, omega_II = 2 omega_I, gamma = hbar = 1\n")
    # Aquí no ordenamos las energías: cada columna sigue a un mismo
    # estado (los dos alineados y los dos del bloque 2x2), así las
    # líneas de la gráfica se cruzan donde de verdad se cruzan.
    archivo.write("# omega_I  E_arriba_arriba  E_bloque_menor  E_bloque_mayor  E_abajo_abajo\n")
    for k in range(101):
        omega = 0.02 * k
        H = dos_espines(omega, 2 * omega, gamma).data
        menor, mayor = eigenvalores_bloque_2x2(H[1][1], H[1][2], H[2][2])
        archivo.write(f"{omega:.4f} {H[0][0]:.12f} {menor:.12f} {mayor:.12f} {H[3][3]:.12f}\n")
print(f"\nEnergías de dos espines guardadas en {ruta.name}")

ruta = CARPETA_DATOS / "tres_espines.dat"
faltantes = []
with open(ruta, "w") as archivo:
    archivo.write("# tres espines, omega_II = 2 omega_I, omega_III = 3 omega_I, gamma = hbar = 1\n")
    archivo.write("# omega_I  E_0 ... E_7  ('?' = el barrido no lo encontró)\n")
    for k in range(81):
        omega = 0.025 * k
        H = tres_espines(omega, 2 * omega, 3 * omega, gamma)
        aproximados = eigenvalores_por_barrido(H.data, puntos=1000, rondas=3)
        energias = [iteracion_inversa(H.data, sigma)[0] for sigma in aproximados]
        escribir_renglon(archivo, omega, energias, 8)
        if len(energias) < 8:
            faltantes.append((omega, 8 - len(energias)))
print(f"Energías de tres espines guardadas en {ruta.name}")
print("Valores de omega_I donde el barrido no encontró los 8 eigenvalores:")
for omega, cuantos in faltantes:
    print(f"  omega_I = {omega:.2f}: faltan {cuantos}")
