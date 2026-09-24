#!/usr/bin/env python3
"""Segunda derivada con diferencias finitas.

Ejemplos ejecutables de la sección "Segunda derivada" de notas.md:
la diferencia central para f'', que sale de *sumar* las series de
Taylor de f(x0+h/2) y f(x0-h/2) (en vez de restarlas, como para f'),
con error de truncamiento O(h^2). Al final, un barrido de h sobre
f(x) = sin(x^2) muestra que el error mínimo es mucho peor que el de
la primera derivada: el error de redondeo ahora se divide entre h^2.
"""

from pathlib import Path

from fiscomp.funciones_especiales import coseno, seno
from fiscomp.precision_numerica import EPS, error_relativo

CARPETA_DATOS = Path(__file__).resolve().parent / "datos"

###############################################
# Diferencia central para la segunda derivada
###############################################


def h_optimo_segunda(f_x0: float = 1.0, f4_x0: float = 1.0) -> float:
    """Paso h óptimo para diff2_central (ver notas.md):

        h_opt = (768 * EPS * |f(x0)| / |f^(4)(x0)|) ** (1/4)

    f_x0: valor de f(x0) (distinto de 0).
    f4_x0: valor de la cuarta derivada f^(4)(x0) (distinto de 0).
    Por defecto ambos valen 1.0 (la fórmula "de juguete").
    """
    return (768.0 * EPS * abs(f_x0) / abs(f4_x0)) ** (1.0 / 4.0)


def diff2_central(f, x0: float, h: float = h_optimo_segunda()) -> float:
    """Diferencia central para la segunda derivada:

        4 (f(x0+h/2) + f(x0-h/2) - 2 f(x0)) / h^2

    Usa f en x0-h/2, x0 y x0+h/2 (tres evaluaciones); error de
    truncamiento O(h^2).

    f: función de una variable, que recibe un flotante y regresa un
       flotante.
    x0: punto donde se aproxima la segunda derivada.
    h: tamaño de paso; si no se da, se usa h_optimo_segunda().
    """
    return 4.0 * (f(x0 + h / 2.0) + f(x0 - h / 2.0) - 2.0 * f(x0)) / h**2


# Respuesta del "Para pensar" de notas.md: la cuarta derivada.
# Las dos sumas de series de Taylor son
#     S1 = f(x0+h/2) + f(x0-h/2) = 2f + (h^2/4) f'' + (h^4/192) f^(4) + (h^6/23040) f^(6) + ...
#     S2 = f(x0+h)   + f(x0-h)   = 2f +  h^2    f'' + (h^4/12)  f^(4) + (h^6/360)   f^(6) + ...
# El término de f'' de S2 es 4 veces el de S1, así que en S2 - 4 S1 se
# cancela; sumando además 6f se cancela también f:
#     S2 - 4 S1 + 6f = (h^4/16) f^(4) + (h^6/384) f^(6) + ...
# y despejando:
#     f^(4)(x0) = 16 [f(x0+h) + f(x0-h) - 4 f(x0+h/2) - 4 f(x0-h/2) + 6 f(x0)] / h^4
#                 - (h^2/24) f^(6)(x0) - ...
# Error O(h^2). Los coeficientes (1, 1, -4, -4, 6) suman cero. El
# error de redondeo es (1+1+4+4+6) EPS |f| * 16/h^4 = 256 EPS |f| / h^4;
# minimizando E = (h^2/24)|f^(6)| + 256 EPS |f| / h^4 da
# h_opt = (12288 EPS |f| / |f^(6)|)^(1/6), ~1e-2 con f y f^(6) ~ 1.


def diff4_central(f, x0: float, h: float = (12288.0 * EPS) ** (1.0 / 6.0)) -> float:
    """Diferencia central para la cuarta derivada:

        16 (f(x0+h) + f(x0-h) - 4 f(x0+h/2) - 4 f(x0-h/2) + 6 f(x0)) / h^4

    Usa f en x0, x0 ± h/2 y x0 ± h (cinco evaluaciones); error de
    truncamiento O(h^2).

    f: función de una variable, que recibe un flotante y regresa un
       flotante.
    x0: punto donde se aproxima la cuarta derivada.
    h: tamaño de paso; si no se da, se usa (12288 * EPS) ** (1/6).
    """
    return (
        16.0
        * (
            f(x0 + h)
            + f(x0 - h)
            - 4.0 * f(x0 + h / 2.0)
            - 4.0 * f(x0 - h / 2.0)
            + 6.0 * f(x0)
        )
        / h**4
    )


###############################################
# Funciones de prueba y sus segundas derivadas exactas
###############################################


def const_5(x):
    return 5.0


def const_5_segunda(x):
    return 0.0


def ident(x):
    return x


def ident_segunda(x):
    return 0.0


def sqr(x):
    return x**2


def sqr_segunda(x):
    return 2.0


def sin_x2(x):
    return seno(x**2)


def sin_x2_segunda(x):
    # f'(x) = 2x cos(x^2), y derivando otra vez (regla del producto):
    return 2 * coseno(x**2) - 4 * x**2 * seno(x**2)


def sin_x2_cuarta(x):
    return -12 * seno(x**2) - 48 * x**2 * coseno(x**2) + 16 * x**4 * seno(x**2)


FUNCIONES_DE_PRUEBA = (
    ("f(x)=5", const_5, const_5_segunda),
    ("f(x)=x", ident, ident_segunda),
    ("f(x)=x^2", sqr, sqr_segunda),
    ("f(x)=sin(x^2)", sin_x2, sin_x2_segunda),
)

###############################################
# Comparar contra la segunda derivada exacta
###############################################

# x0 = 1 y no 6 (como en solucion_diferencias_finitas.py): seno y
# coseno de fiscomp pierden precisión para argumentos grandes, y en
# x0 = 6 evaluaríamos seno(36).
x0 = 1.0

print(f"Segunda derivada en x0={x0}, con el h óptimo por defecto (h={h_optimo_segunda():.2e}):\n")
print(f"{'función':<16}{'exacta':>14}{'aproximada':>18}{'err. relativo':>16}")
for nombre, f, f_segunda in FUNCIONES_DE_PRUEBA:
    exacta = f_segunda(x0)
    aproximada = diff2_central(f, x0)
    print(f"{nombre:<16}{exacta:14.6f}{aproximada:18.10f}{error_relativo(aproximada, exacta):16.2e}")

# Para f(x)=5 y f(x)=x la segunda derivada es 0, y la aproximación da
# (casi) 0. Para f(x)=x^2 el error de truncamiento es cero (depende de
# f^(4), que para x^2 vale 0) y solo queda redondeo. Para sin(x^2) el
# error queda en ~1e-7/1e-8, como predice notas.md: mucho peor que el
# ~1e-10 de la diferencia central para f', aunque ambas son O(h^2).
print(30 * "=")

###############################################
# Barrido de h
###############################################

# Igual que para la primera derivada: h = 1, 1/2, 1/4, ... Ahora el
# error de redondeo crece como 1/h^2, así que el mínimo llega antes
# (con h más grande) y es más alto.

exacta = sin_x2_segunda(x0)

filas = []
h = 1.0
for _ in range(50):
    aproximada = diff2_central(sin_x2, x0, h)
    filas.append((h, aproximada, error_relativo(aproximada, exacta)))
    h *= 0.5

print(f"\nBarrido de h para f''(x) de f(x)=sin(x^2) en x0={x0} (exacta = {exacta:.12f}):\n")
print(f"{'h':>12}{'aproximada':>20}{'err. relativo':>16}")
for h, aproximada, error in filas:
    print(f"{h:12.2e}{aproximada:20.10f}{error:16.2e}")

mejor = min(filas, key=lambda fila: fila[2])
f_x0 = sin_x2(x0)
f4_x0 = sin_x2_cuarta(x0)
print(f"\nMejor h (barrido):                 {mejor[0]:.2e}  (error {mejor[2]:.2e})")
print(f"h_opt de juguete ~ (768*EPS)^(1/4) = {h_optimo_segunda():.2e}")
print(f"h_opt con f, f^(4) reales          = {h_optimo_segunda(f_x0, f4_x0):.2e}")

# Como en el Ejercicio 4 de la práctica: la h del barrido y la de la
# fórmula coinciden en orden de magnitud, no exactamente. Con los
# valores reales (f^(4)(1) ~ -22.6, lejos de 1) la fórmula se acerca
# más a lo medido.
#
# Para h < ~1e-8 el error vale exactamente 1: f(x0 ± h/2) ya no se
# distingue de f(x0) en punto flotante, el numerador sale exactamente
# 0 y la "segunda derivada" también. Es la cancelación catastrófica
# (unidad 06) llevada al extremo.
print(30 * "=")

###############################################
# Cuarta derivada (respuesta del "Para pensar")
###############################################

exacta_cuarta = sin_x2_cuarta(x0)
aproximada_cuarta = diff4_central(sin_x2, x0)
print(f"\nCuarta derivada de sin(x^2) en x0={x0}:")
print(f"  exacta     = {exacta_cuarta:.10f}")
print(f"  aproximada = {aproximada_cuarta:.10f}  (error {error_relativo(aproximada_cuarta, exacta_cuarta):.2e})")

# El error es todavía peor: ~2e-4 con el h de juguete, y ~1e-5 en el
# mejor de los casos (h ~ 2e-3 en un barrido de potencias de 1/2),
# porque el redondeo ahora se divide entre h^4 y f^(6)(1) ~ 750 está
# muy lejos de 1. Cada derivada adicional cuesta varios dígitos.
print(30 * "=")

###############################################
# Guardar el barrido en un archivo, para graficarlo
###############################################

CARPETA_DATOS.mkdir(exist_ok=True)
ruta_datos = CARPETA_DATOS / "segunda_derivada_sin_x2.dat"
with open(ruta_datos, "w") as archivo:
    archivo.write("# h  diff2_central  error_relativo\n")
    for h, aproximada, error in filas:
        archivo.write(f"{h} {aproximada} {error}\n")

print(f"Barrido guardado en {ruta_datos.name}")
