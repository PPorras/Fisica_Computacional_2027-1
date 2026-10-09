#!/usr/bin/env python3
"""Ejercicio 4: compara seno/coseno/exponencial/ln contra math, y
guarda el reporte en reporte_ejercicio4.txt (unidad 04).
"""

import math
from pathlib import Path

from fiscomp import funciones_especiales as fe
from fiscomp.precision_numerica import error_relativo

CASOS = {
    "seno": (fe.seno, math.sin, [0.0, 0.5, 1.0, math.pi / 2, math.pi, 3 * math.pi]),
    "coseno": (fe.coseno, math.cos, [0.0, 0.5, 1.0, math.pi / 2, math.pi, 3 * math.pi]),
    "exponencial": (fe.exponencial, math.exp, [0.0, 1.0, -1.0, 2.5, 5.0, -20.0]),
    "ln": (fe.ln, math.log, [1.0, 0.1, 0.5, 2.0, 5.0, 100.0]),
}

ruta_reporte = Path(__file__).resolve().parent / "reporte_ejercicio4.txt"

with open(ruta_reporte, "w") as archivo:
    for nombre_funcion, (funcion_propia, funcion_math, valores) in CASOS.items():
        archivo.write(f"{nombre_funcion}\n")
        archivo.write("-" * len(nombre_funcion) + "\n")
        for x in valores:
            aproximado = funcion_propia(x)
            exacto = funcion_math(x)
            error = error_relativo(aproximado, exacto)
            archivo.write(
                f"  x={x:>12.6f}  propio={aproximado:.12f}  "
                f"math={exacto:.12f}  error_relativo={error:.3e}\n"
            )
        archivo.write("\n")

    archivo.write("Conclusiones\n")
    archivo.write("------------\n")
    archivo.write(
        "El error relativo es sorprendentemente alto en seno(pi), "
        "seno(3 pi) y coseno(pi/2) (de 1e-2 a 1e1, en vez de ~1e-16). "
        "En esos puntos el valor 'real' (math.sin/math.cos) no es "
        "exactamente 0, sino un número muy chico, porque pi no se puede "
        "representar exactamente en binario: math.pi es solo el float "
        "más cercano. Como error_relativo divide entre ese valor "
        "'exacto' casi cero, un error absoluto minúsculo (de los últimos "
        "bits) se convierte en un error relativo enorme. No es que la "
        "aproximación sea mala ahí: es una debilidad del error relativo "
        "como métrica cerca de los ceros de la función.\n"
    )
    archivo.write(
        "En exponencial(-20) el error es del orden de 1e-16 porque "
        "fiscomp calcula e^x para x < 0 como 1/e^|x|. Si se suma la "
        "serie de Taylor directamente con x = -20, los términos llegan "
        "a ~4e7, con signos alternados, y casi se cancelan entre sí "
        "(cancelación catastrófica, unidad 06): el resultado, "
        "e^-20 ~ 2e-9, sale con un error relativo de ~2, es decir, "
        "sin ningún dígito correcto.\n"
    )

print(f"Reporte guardado en {ruta_reporte}")
