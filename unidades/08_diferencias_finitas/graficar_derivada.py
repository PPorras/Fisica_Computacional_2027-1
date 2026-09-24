#!/usr/bin/env python3
"""Grafica, con matplotlib, el error relativo de diff_forward y
diff_central contra h.

Versión en Python de graficar_derivada.gp (gnuplot); hace lo mismo:
error relativo contra h en escala log-log, con una línea vertical en
la h_opt teórica de cada método (ver notas.md).

Si además existe datos/derivada_sin_x2_precisas.dat (la parte de
"diferencias finitas más precisas" de notas.md), abre una segunda
ventana que compara los cinco métodos: adelante O(h), central y
segunda adelante O(h^2), segunda central y Richardson O(h^4).

matplotlib no es parte del entorno base del curso (README.md solo
instala fiscomp). Para instalarlo, con el entorno virtual activado
(desde la raíz del repositorio):

    source .venv/bin/activate
    python3 -m pip install matplotlib

Instrucciones detalladas (Linux, macOS y Windows) y solución a
problemas comunes en recursos/notas_matplotlib.md.

Uso (después de completar la práctica en diferencias_finitas.py --
ver practica.md -- y correrlo, para que genere el .dat que se lee
aquí):

    python3 graficar_derivada.py
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
RUTA_DATOS = CARPETA_DATOS / "derivada_sin_x2.dat"
RUTA_DATOS_PRECISAS = CARPETA_DATOS / "derivada_sin_x2_precisas.dat"

# Epsilon de la máquina en doble precisión IEEE 754, escrito a mano
# (es un valor fijo del estándar, no depende de la corrida): la misma
# constante que fiscomp.precision_numerica.EPS.
EPS = 2.220446049250313e-16

# h_opt "de juguete" de notas.md (asumiendo f y sus derivadas de orden
# 1): el mínimo teórico de cada curva de error, donde el error de
# truncamiento y el de redondeo quedan balanceados.
H_OPT_ADELANTE = (4.0 * EPS) ** 0.5
H_OPT_CENTRAL = (24.0 * EPS) ** (1.0 / 3.0)
H_OPT_ADELANTE2 = (48.0 * EPS) ** (1.0 / 3.0)
H_OPT_CENTRAL2 = (1120.0 * EPS / 9.0) ** (1.0 / 5.0)

# Curvas de la segunda gráfica, en el orden de las columnas de
# derivada_sin_x2_precisas.dat (después de h): etiqueta, color y
# estilo de línea. Richardson va punteada porque casi se encima con la
# segunda central (ambas son O(h^4)).
METODOS_PRECISOS = (
    ("adelante, O(h)", "red", "o-"),
    ("central, O(h^2)", "blue", "o-"),
    ("segunda adelante, O(h^2)", "darkorange", "o-"),
    ("segunda central, O(h^4)", "darkgreen", "o-"),
    ("Richardson, O(h^4)", "purple", "--"),
)


def leer_datos(ruta):
    """Lee el .dat generado por la práctica (Ejercicio 3).

    Columnas esperadas: h  diff_forward  error_forward  diff_central
    error_central, separadas por espacios; las líneas que empiezan
    con '#' se ignoran.
    """
    h_vals, error_adelante, error_central = [], [], []
    with open(ruta) as archivo:
        for linea in archivo:
            if linea.startswith("#") or not linea.strip():
                continue
            columnas = linea.split()
            h_vals.append(float(columnas[0]))
            error_adelante.append(float(columnas[2]))
            error_central.append(float(columnas[4]))
    return h_vals, error_adelante, error_central


def leer_datos_precisos(ruta):
    """Lee derivada_sin_x2_precisas.dat.

    Columnas: h  error_forward  error_central  error_forward2
    error_central2  error_richardson. Regresa la lista de h y una lista
    de errores por cada método (en ese orden).
    """
    h_vals = []
    errores = [[] for _ in METODOS_PRECISOS]
    with open(ruta) as archivo:
        for linea in archivo:
            if linea.startswith("#") or not linea.strip():
                continue
            columnas = linea.split()
            h_vals.append(float(columnas[0]))
            for i, lista in enumerate(errores):
                lista.append(float(columnas[i + 1]))
    return h_vals, errores


def graficar_precisos():
    """Segunda figura: el error de los cinco métodos contra h."""
    h_vals, errores = leer_datos_precisos(RUTA_DATOS_PRECISAS)

    fig, ax = plt.subplots()
    for (etiqueta, color, estilo), error in zip(METODOS_PRECISOS, errores):
        ax.plot(h_vals, error, estilo, color=color, markersize=3, label=etiqueta)

    ax.axvline(H_OPT_ADELANTE2, color="darkorange", linestyle="--", label="h_opt segunda adelante")
    ax.axvline(H_OPT_CENTRAL2, color="darkgreen", linestyle="--", label="h_opt segunda central")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("h")
    ax.set_ylabel("error relativo")
    ax.set_title("Diferencias finitas más precisas: error vs. h")
    ax.grid(True, which="both", linestyle=":", linewidth=0.5)
    ax.legend(loc="upper center", fontsize="small")


def main():
    if not RUTA_DATOS.exists():
        raise SystemExit(
            f"No encontré {RUTA_DATOS}. Completen la práctica en "
            "diferencias_finitas.py primero (debe generar ese archivo; "
            "ver practica.md)."
        )

    h_vals, error_adelante, error_central = leer_datos(RUTA_DATOS)

    fig, ax = plt.subplots()
    ax.plot(h_vals, error_adelante, "o-", color="red", markersize=3, label="adelante, O(h)")
    ax.plot(h_vals, error_central, "o-", color="blue", markersize=3, label="central, O(h^2)")

    ax.axvline(H_OPT_ADELANTE, color="red", linestyle="--", label="h_opt adelante")
    ax.axvline(H_OPT_CENTRAL, color="blue", linestyle="--", label="h_opt central")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("h")
    ax.set_ylabel("error relativo")
    ax.set_title("Error de diferencias finitas vs. tamaño de paso h")
    ax.grid(True, which="both", linestyle=":", linewidth=0.5)
    ax.legend()

    if RUTA_DATOS_PRECISAS.exists():
        graficar_precisos()

    plt.show()


if __name__ == "__main__":
    main()
