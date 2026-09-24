#!/usr/bin/env python3
"""Grafica, con matplotlib, el error relativo de diff2_central (la
diferencia central para la segunda derivada) contra h.

Versión en Python de graficar_segunda_derivada.gp (gnuplot); hace lo
mismo: error relativo contra h en escala log-log, con una línea
vertical en la h_opt teórica (ver notas.md).

matplotlib no es parte del entorno base del curso. Para instalarlo,
con el entorno virtual activado (desde la raíz del repositorio):

    source .venv/bin/activate
    python3 -m pip install matplotlib

Instrucciones detalladas en recursos/notas_matplotlib.md.

Uso (después de correr segunda_derivada.py, que genera el .dat que se
lee aquí):

    python3 graficar_segunda_derivada.py
"""

from pathlib import Path

try:
    import matplotlib.pyplot as plt
except ImportError as error:
    raise SystemExit(
        "Este script necesita matplotlib, que no está instalado en este "
        "entorno virtual. Actívenlo (source .venv/bin/activate, desde la "
        "raíz del repositorio) e instálenlo con: python3 -m pip install matplotlib\n"
        "Instrucciones detalladas en recursos/notas_matplotlib.md."
    ) from error

CARPETA_DATOS = Path(__file__).resolve().parent / "datos"
RUTA_DATOS = CARPETA_DATOS / "segunda_derivada_sin_x2.dat"

# Epsilon de la máquina en doble precisión IEEE 754 (la misma
# constante que fiscomp.precision_numerica.EPS).
EPS = 2.220446049250313e-16

# h_opt "de juguete" de notas.md (f y f^(4) de orden 1).
H_OPT = (768.0 * EPS) ** (1.0 / 4.0)


def leer_datos(ruta):
    """Lee el .dat generado por segunda_derivada.py.

    Columnas: h  diff2_central  error_relativo, separadas por
    espacios; las líneas que empiezan con '#' se ignoran.
    """
    h_vals, errores = [], []
    with open(ruta) as archivo:
        for linea in archivo:
            if linea.startswith("#") or not linea.strip():
                continue
            columnas = linea.split()
            h_vals.append(float(columnas[0]))
            errores.append(float(columnas[2]))
    return h_vals, errores


def main():
    if not RUTA_DATOS.exists():
        raise SystemExit(f"No encontré {RUTA_DATOS}. Corran segunda_derivada.py primero.")

    h_vals, errores = leer_datos(RUTA_DATOS)

    fig, ax = plt.subplots()
    ax.plot(h_vals, errores, "o-", color="darkgreen", markersize=3, label="diferencia central para f'', O(h^2)")
    ax.axvline(H_OPT, color="darkgreen", linestyle="--", label="h_opt")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("h")
    ax.set_ylabel("error relativo")
    ax.set_title("Segunda derivada: error vs. tamaño de paso h")
    ax.grid(True, which="both", linestyle=":", linewidth=0.5)
    ax.legend(loc="lower left")

    plt.show()


if __name__ == "__main__":
    main()
