"""Derivadas numéricas con diferencias finitas.

Ver unidades/08_diferencias_finitas/notas.md para la deducción de cada
fórmula (a partir de la serie de Taylor) y su análisis de error.

Por ahora solo contiene la segunda derivada, que se usa en varios
lugares del curso (la unidad 08 y la simulación de energía cinética
local). Las diferencias hacia adelante, hacia atrás y central para la
primera derivada no están aquí porque son parte de la práctica 3.
"""

from fiscomp.precision_numerica import EPS


def h_optimo_segunda(f_x0: float = 1.0, f4_x0: float = 1.0) -> float:
    """Paso h óptimo para diff2_central (ver notas.md de la unidad 08):

        h_opt = (768 * EPS * |f(x0)| / |f^(4)(x0)|) ** (1/4)

    f_x0: valor de f(x0) (distinto de 0).
    f4_x0: valor de la cuarta derivada f^(4)(x0) (distinto de 0).
    Por defecto ambos valen 1.0 (la fórmula "de juguete"), que da
    h_opt ~ 6e-4.
    """
    return (768.0 * EPS * abs(f_x0) / abs(f4_x0)) ** (1.0 / 4.0)


def diff2_central(f, x0: float, h: float = h_optimo_segunda()):
    """Diferencia central para la segunda derivada:

        f''(x0) ~ 4 (f(x0+h/2) + f(x0-h/2) - 2 f(x0)) / h^2

    Usa f en x0-h/2, x0 y x0+h/2 (tres evaluaciones); error de
    truncamiento O(h^2). Funciona igual si f regresa números reales o
    complejos.

    f: función de una variable, que recibe un flotante y regresa un
       número.
    x0: punto donde se aproxima la segunda derivada.
    h: tamaño de paso; si no se da, se usa h_optimo_segunda().

    Se evalúa f directamente en x0 + h/2 y x0 - h/2, sin modificar x0.
    Si en cambio se hiciera x0 += h/2, luego x0 -= h, etc., el error de
    redondeo de cada suma se iría acumulando (unidad 06) y al final x0
    ya no sería exactamente el punto original.
    """
    return 4.0 * (f(x0 + h / 2.0) + f(x0 - h / 2.0) - 2.0 * f(x0)) / h**2
