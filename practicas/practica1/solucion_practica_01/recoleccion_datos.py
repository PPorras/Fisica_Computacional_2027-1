#!/usr/bin/env python3
"""Ejercicio 1: recolección de datos del experimento de Millikan
(carga del electrón), usando input() y todos los tipos de datos de la
unidad 00, más una estimación estadística de `e` (Parte B). Genera
reporte_recoleccion.txt (unidad 04) con el resumen.
"""

import math
from pathlib import Path

from fiscomp.precision_numerica import error_relativo

CARGA_ELECTRON_ACEPTADA = 1.602176634e-19  # Coulombs (valor moderno, CODATA)


def estimar_carga_electron(cargas_medidas):
    """Estima `e` a partir de varias cargas medidas (método de Millikan).

    Regresa (estimacion, desviacion_estandar, numeros_de_electrones):
    la carga elemental estimada, qué tan dispersas quedaron las
    estimaciones por gota, y cuántas cargas elementales se le asignó
    a cada gota (mismo orden que `cargas_medidas`).
    """
    e_aproximada = min(cargas_medidas)

    numeros_de_electrones = [round(carga / e_aproximada) for carga in cargas_medidas]
    estimaciones_por_gota = [
        carga / n for carga, n in zip(cargas_medidas, numeros_de_electrones)
    ]

    estimacion = sum(estimaciones_por_gota) / len(estimaciones_por_gota)
    varianza = sum((e - estimacion) ** 2 for e in estimaciones_por_gota) / len(
        estimaciones_por_gota
    )
    desviacion_estandar = math.sqrt(varianza)

    return estimacion, desviacion_estandar, numeros_de_electrones


# --- Parte A: recolección de datos ---

# str
nombre_experimento = input("Nombre del experimento: ")
responsable = input("Responsable: ")

# tupla: condiciones del experimento, fijas durante toda la corrida
voltaje_aplicado = float(input("Voltaje aplicado entre placas (V): "))
distancia_entre_placas = float(input("Distancia entre placas (m): "))
viscosidad_aceite = float(input("Viscosidad del aceite (Pa·s): "))
condiciones_experimento = (voltaje_aplicado, distancia_entre_placas, viscosidad_aceite)

# int: número de gotas a medir
numero_de_gotas = int(input("¿Cuántas gotas vas a medir?: "))

# lista: carga medida de cada gota (en Coulombs)
cargas_medidas = []
for i in range(numero_de_gotas):
    carga = float(input(f"  Carga de la gota {i + 1} (C): "))
    cargas_medidas.append(carga)

# bool: si el experimento se considera válido (cargas positivas y de
# un orden de magnitud razonable)
experimento_valido = all(0 < carga < 1e-17 for carga in cargas_medidas)

# set: valores únicos medidos (por si dos gotas midieron "la misma" carga)
valores_unicos = set(cargas_medidas)

# --- Parte B: estimar la carga del electrón ---

carga_estimada, desviacion_estandar, numeros_de_electrones = estimar_carga_electron(
    cargas_medidas
)
error_vs_aceptado = error_relativo(carga_estimada, CARGA_ELECTRON_ACEPTADA)

# dict: resumen final
resumen = {
    "experimento": nombre_experimento,
    "responsable": responsable,
    "condiciones (voltaje V, distancia m, viscosidad Pa·s)": condiciones_experimento,
    "numero_de_gotas": numero_de_gotas,
    "valores_unicos_medidos": len(valores_unicos),
    "experimento_valido": experimento_valido,
    "carga_electron_estimada (C)": carga_estimada,
    "desviacion_estandar (C)": desviacion_estandar,
    "error_relativo_vs_valor_aceptado": error_vs_aceptado,
}

# --- Reporte ---

ruta_reporte = Path(__file__).resolve().parent / "reporte_recoleccion.txt"
with open(ruta_reporte, "w") as archivo:
    archivo.write(f"Reporte del experimento: {nombre_experimento}\n")
    archivo.write("=" * 40 + "\n")
    for clave, valor in resumen.items():
        archivo.write(f"{clave}: {valor}\n")

    archivo.write("\nDetalle por gota (carga medida, n, estimación individual):\n")
    for i, (carga, n) in enumerate(zip(cargas_medidas, numeros_de_electrones), start=1):
        estimacion_individual = carga / n
        archivo.write(
            f"  {i}. carga={carga:.4e} C  n={n}  "
            f"estimacion_individual={estimacion_individual:.4e} C\n"
        )

print(f"Reporte guardado en {ruta_reporte}")
