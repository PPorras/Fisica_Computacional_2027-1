#!/usr/bin/env python3
"""Funciones de onda de una partícula en una dimensión.

Dos casos (ver README.md de esta carpeta):

- psi_oscilador: eigenfunciones del oscilador armónico cuántico.
- psi_caja: ondas planas de una partícula libre en una caja periódica.

Las dos funciones tienen la **misma interfaz**: reciben la posición x
y un diccionario `parametros` con lo que cada una necesita. Así,
cualquier otra función (como energia_cinetica_local, en
energia_cinetica_local.py) puede usarlas sin saber cuál es cuál ni
cuántos parámetros necesita cada una.

Unidades: hbar = m = 1.
"""

from fiscomp.funciones_especiales import PI, coseno, exponencial, factorial, raiz_cuadrada, seno


def hermite(n: int, x: float) -> float:
    """Polinomio de Hermite H_n(x), con la relación de recurrencia

        H_(j+1)(x) = 2x H_j(x) - 2j H_(j-1)(x),

    empezando por H_0(x) = 1 y H_1(x) = 2x.

    n: grado del polinomio (entero no negativo).
    x: punto donde se evalúa.
    """
    if n == 0:
        return 1.0
    h_anterior = 1.0  # H_0
    h_actual = 2.0 * x  # H_1
    for j in range(1, n):
        h_anterior, h_actual = h_actual, 2.0 * x * h_actual - 2.0 * j * h_anterior
    return h_actual


def psi_oscilador(x: float, parametros: dict) -> float:
    """Eigenfunción normalizada del oscilador armónico cuántico:

        psi_n(x) = (m w / hbar)^(1/4) / sqrt(2^n n! sqrt(pi))
                   * H_n(sqrt(m w / hbar) x) * exp(-alfa (m w / hbar) x^2 / 2)

    x: posición.
    parametros: diccionario con
        "n":         número cuántico (entero no negativo),
        "m_w_hbar":  el cociente m w / hbar,
        "alfa":      factor extra en el exponente de la gaussiana; con
                     alfa = 1 es la eigenfunción exacta (con otro
                     valor sirve como función de prueba variacional).
    """
    n = parametros["n"]
    m_w_hbar = parametros["m_w_hbar"]
    alfa = parametros["alfa"]

    normalizacion = m_w_hbar**0.25 / raiz_cuadrada(2**n * factorial(n) * raiz_cuadrada(PI))
    gaussiana = exponencial(-0.5 * alfa * m_w_hbar * x**2)
    return normalizacion * hermite(n, raiz_cuadrada(m_w_hbar) * x) * gaussiana


def psi_caja(x: float, parametros: dict) -> complex:
    """Onda plana normalizada en una caja periódica de longitud L:

        psi_k(x) = exp(i k x) / sqrt(L),   con k = 2 pi n / L

    x: posición (dentro de la caja, -L/2 <= x <= L/2).
    parametros: diccionario con
        "n": número cuántico (entero: positivo, negativo o cero),
        "L": longitud de la caja.
    """
    n = parametros["n"]
    longitud = parametros["L"]

    # exp(i k x) = cos(k x) + i sin(k x) (fórmula de Euler); 1j es la
    # unidad imaginaria en Python.
    k = 2.0 * PI * n / longitud
    return (coseno(k * x) + 1j * seno(k * x)) / raiz_cuadrada(longitud)


if __name__ == "__main__":
    # Las dos llamadas se ven idénticas: solo cambia el diccionario.
    x = 1.0

    parametros = {"n": 100, "m_w_hbar": 1.0, "alfa": 1.0}
    psi_a = psi_oscilador(x, parametros)

    parametros = {"n": -2, "L": 2 * PI}
    psi_b = psi_caja(x, parametros)

    print(f"psi_oscilador(x={x}) = {psi_a}")
    print(f"psi_caja(x={x})      = {psi_b}")
