#!/usr/bin/env python3
"""Reimplementación propia de funciones matemáticas elementales.

La idea es construir, sin usar el módulo `math` de la librería
estándar, aproximaciones numéricas de funciones como:

- factorial(n)     -- ya implementada
- seno(x)          -- ya implementada, con serie de Taylor
- coseno(x)        -- pendiente (práctica 1)
- exponencial(x)   -- pendiente (práctica 1)
- ln(x)            -- pendiente (práctica 1)
- raiz_cuadrada(x) -- pendiente (práctica 1)

Las funciones basadas en series (seno, coseno, exponencial, ...) usan
EPS (fiscomp.precision_numerica) como criterio de convergencia: se
suman términos mientras el siguiente término siga siendo mayor o
igual que el épsilon de la máquina, y se corta la suma en cuanto deja
de aportar precisión adicional.
"""

from fiscomp.precision_numerica import EPS


def factorial(n):
    """Calcula n! (n factorial) de forma iterativa.

    n debe ser un entero no negativo.
    """
    resultado = 1
    for i in range(2, n + 1):
        resultado *= i
    return resultado


def seno(x, precision=EPS):
    """Aproxima sin(x) con la serie de Taylor alrededor de 0:

        sin(x) = suma_{k=0}^inf (-1)^k * x^(2k+1) / (2k+1)!

    Se suman términos mientras sigan siendo mayores o iguales que
    `precision` (por default, el épsilon de la máquina); en cuanto un
    término es más chico, ya no cambia el resultado y se detiene la suma.

    Nota: esta serie no hace reducción de rango (llevar x a [-pi, pi]
    antes de sumar), así que para |x| grande la precisión se degrada
    por cancelación entre términos grandes de signos alternados.
    """
    suma = 0.0
    k = 0
    while True:
        termino = (-1) ** k * x ** (2 * k + 1) / factorial(2 * k + 1)
        if abs(termino) < precision:
            break
        suma += termino
        k += 1
    return suma


if __name__ == "__main__":
    import math

    from fiscomp.precision_numerica import error_relativo

    print(f"factorial(5) = {factorial(5)}")
    print(f"math.factorial(5) = {math.factorial(5)}")

    for x in (0.0, 0.5, 1.0, math.pi / 2, math.pi, 3 * math.pi):
        aproximado = seno(x)
        exacto = math.sin(x)
        print(
            f"seno({x:.4f}) = {aproximado:.12f}  "
            f"math.sin = {exacto:.12f}  "
            f"error_relativo = {error_relativo(aproximado, exacto):.2e}"
        )
