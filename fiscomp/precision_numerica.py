#! /usr/bin/env python3
"""Herramientas sobre la representación de punto flotante.

Ilustra dos de las fuentes de error numérico más comunes:

- Error de redondeo: epsilon_maquina() lo estima directamente, y
  son_iguales() ofrece una forma de comparar flotantes que lo tolera.
- Error de truncamiento: error_relativo() sirve para medir qué tan
  buena es una aproximación numérica (como una serie truncada) frente
  al valor exacto.

(La incertidumbre experimental, la tercera fuente de error, no es un
tema de aritmética de punto flotante y no se trata en este módulo.)
"""


def epsilon_maquina():
    """Estima el épsilon de la máquina.

    El épsilon de la máquina es el menor número positivo tal que,
    en la aritmética de punto flotante del sistema, 1.0 + epsilon
    ya es distinguible de 1.0. Se calcula partiendo de epsilon=1.0
    y dividiéndolo entre 2 repetidamente, hasta que sumarlo a 1.0
    deja de tener efecto (por el redondeo de punto flotante).

    Esta estimación sirve como criterio de convergencia para las
    aproximaciones numéricas de `fiscomp.funciones_especiales`
    (exp, sin, cos, ...): al sumar los términos de una serie, se
    detiene la suma cuando el siguiente término es menor que el
    épsilon de la máquina, punto en el que ya no aporta precisión
    adicional al resultado.
    """
    epsilon = 1.0
    while 1.0 + epsilon / 2.0 > 1.0:
        epsilon /= 2.0
    return epsilon


EPS = epsilon_maquina()


def son_iguales(a, b, tolerancia=EPS):
    """Compara dos flotantes con una tolerancia, en vez de usar `==`.

    Por los errores de redondeo de la aritmética de punto flotante,
    dos resultados que matemáticamente son iguales pueden diferir en
    los últimos bits (el ejemplo clásico: 0.1 + 0.2 != 0.3). En vez de
    comparar con `==`, se consideran iguales si su diferencia absoluta
    es menor que una tolerancia.
    """
    return abs(a - b) < tolerancia


def error_relativo(aproximado, exacto):
    """Calcula el error relativo entre un valor aproximado y el exacto.

    Sirve para medir qué tan buena es una aproximación numérica (por
    ejemplo, una serie truncada) sin depender de la escala del valor:
    un error absoluto de 0.001 es enorme si el valor exacto es 1e-6,
    pero insignificante si el valor exacto es 1e6.

    Si el valor exacto es 0, se regresa el error absoluto en su lugar
    (para no dividir entre cero).
    """
    if exacto == 0.0:
        return abs(aproximado - exacto)
    return abs(aproximado - exacto) / abs(exacto)


if __name__ == "__main__":
    import math
    import sys

    print(f"epsilon_maquina() = {EPS}")
    print(f"sys.float_info.epsilon = {sys.float_info.epsilon}")
    print(f"¿Coinciden? {EPS == sys.float_info.epsilon}")

    print(f"1.0 + EPS != 1.0: {1.0 + EPS != 1.0}")
    print(f"1.0 + EPS / 2.0 != 1.0: {1.0 + EPS / 2.0 != 1.0}")

    # Error de redondeo: el ejemplo clásico
    print(f"0.1 + 0.2 = {0.1 + 0.2!r}")
    print(f"0.1 + 0.2 == 0.3: {0.1 + 0.2 == 0.3}")
    print(f"son_iguales(0.1 + 0.2, 0.3): {son_iguales(0.1 + 0.2, 0.3)}")

    # Error de truncamiento: una serie de Taylor truncada para exp(1)
    aproximacion_exp = sum(1 / math.factorial(n) for n in range(10))
    print(f"Serie de Taylor truncada para exp(1): {aproximacion_exp}")
    print(f"math.exp(1): {math.exp(1)}")
    print(f"error_relativo: {error_relativo(aproximacion_exp, math.exp(1)):.2e}")
