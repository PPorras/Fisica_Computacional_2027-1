#!/usr/bin/env python3
"""Ejercicio 2: calcula PI con la serie de Leibniz y lo guarda en
fiscomp/constantes.py.

pi/4 = 1 - 1/3 + 1/5 - 1/7 + 1/9 - ...

La serie converge muy lento (el error baja como 1/N), así que llegar
al épsilon de la máquina no es práctico -- aquí se suman 20 millones
de términos (~4 s) y se compara contra math.pi para reportar qué tan
cerca quedó, en vez de sumar hasta EPS.
"""

import math
import time

N = 20_000_000


def calcular_pi(n_terminos):
    suma = 0.0
    signo = 1.0
    for k in range(n_terminos):
        suma += signo / (2 * k + 1)
        signo = -signo
    return 4 * suma


if __name__ == "__main__":
    from fiscomp.precision_numerica import error_relativo

    inicio = time.time()
    pi_aproximado = calcular_pi(N)
    duracion = time.time() - inicio

    print(f"N = {N:,} términos")
    print(f"pi_aproximado = {pi_aproximado!r}")
    print(f"math.pi       = {math.pi!r}")
    print(f"error_relativo = {error_relativo(pi_aproximado, math.pi):.3e}")
    print(f"tiempo = {duracion:.2f} s")
    # error_relativo ~ 1.6e-08 con N=20,000,000 (bien por debajo del
    # 1e-4 que exige la prueba); duplicar N solo reduce el error a la
    # mitad porque la serie converge como 1/N, así que llegar a EPS
    # (~1e-16) tomaría billones de términos -- por eso no se hace.
