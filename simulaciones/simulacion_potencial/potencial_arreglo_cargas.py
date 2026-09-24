#!/usr/bin/env python3
"""Potencial electrostático de un arreglo de 36 cargas puntuales.
Un arreglo de cargas de signo alternado y magnitud
variable. Se arman 36 cargas en una rejilla 6x6 (signo tipo tablero de
ajedrez, magnitud creciente hacia una esquina), se evalúa el potencial
total por superposición (potencial_electrostatico, de
simulacion.py en esta misma carpeta) en una rejilla más fina, y el resultado se guarda en
un .dat con el formato de "grid data" que espera gnuplot: bloques de x
constante, separados por un renglón en blanco.

Para graficarlo en gnuplot, por ejemplo:

    set pm3d map
    set cbrange [-2e2:2e2]   # el potencial diverge junto a cada carga;
                              # sin acotar el rango de color, esos picos
                              # dominan la escala y no se ve el resto
    splot 'datos/potencial_arreglo_cargas.dat' using 1:2:3 with pm3d

o, para curvas de nivel:

    set contour base
    set view map
    unset surface
    splot 'datos/potencial_arreglo_cargas.dat' using 1:2:3
"""

from pathlib import Path

from fiscomp.vectores import VectorND
from simulacion import potencial_electrostatico


def arreglo_36_cargas(espaciado=1.0, carga_base=1e-9):
    """Arma 36 cargas puntuales en una rejilla 6x6 centrada en el
    origen: signo alternado (patrón de tablero de ajedrez) y magnitud
    creciente hacia la esquina superior derecha, como en la Fig. 2.4.

    Regresa (cargas, posiciones): dos listas paralelas, del formato
    que espera potencial_electrostatico.
    """
    coordenadas = [(i - 2.5) * espaciado for i in range(6)]

    cargas = []
    posiciones = []
    for i, x in enumerate(coordenadas):
        for j, y in enumerate(coordenadas):
            signo = 1 if (i + j) % 2 == 0 else -1
            magnitud = carga_base * (1 + 0.5 * (i + j))
            cargas.append(signo * magnitud)
            posiciones.append(VectorND([x, y]))
    return cargas, posiciones


def guardar_potencial_dat(ruta_relativa, cargas, posiciones, limite=5.0, n=101):
    """Evalúa el potencial en una rejilla de n x n puntos en
    [-limite, limite]^2 y lo guarda como "x y potencial" por renglón,
    con un renglón en blanco entre cada bloque de x constante (el
    formato de rejilla que gnuplot espera para splot/pm3d/contour).

    Si un punto de evaluación coincide exactamente con una carga (el
    potencial diverge ahí), se escribe "nan"; gnuplot lo interpreta
    como un hueco en la rejilla en vez de tronar.
    """
    ruta = Path(ruta_relativa)
    ruta.parent.mkdir(parents=True, exist_ok=True)

    paso = 2 * limite / (n - 1)
    with open(ruta, "w") as archivo:
        for i in range(n):
            x = -limite + i * paso
            for j in range(n):
                y = -limite + j * paso
                try:
                    v = potencial_electrostatico(VectorND([x, y]), cargas, posiciones)
                except ZeroDivisionError:
                    v = float("nan")
                archivo.write(f"{x:.6f} {y:.6f} {v:.6e}\n")
            archivo.write("\n")

    print(f"Datos guardados en {ruta.resolve()}")


if __name__ == "__main__":
    cargas, posiciones = arreglo_36_cargas()
    print(f"{len(cargas)} cargas armadas en una rejilla 6x6.")
    print(f"carga total del arreglo = {sum(cargas):.3e} C")

    guardar_potencial_dat("datos/potencial_arreglo_cargas.dat", cargas, posiciones)
