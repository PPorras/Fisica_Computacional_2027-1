#!/usr/bin/env python3
"""Diferenciación numérica con diferencias finitas.

Ejemplos ejecutables de lo visto en notas.md: diff_forward,
diff_backward y diff_central se derivan cada una de la serie de
Taylor de f alrededor de x0, con distinto orden de error de
truncamiento (O(h), O(h) y O(h^2), respectivamente). Al final, un
barrido de h sobre f(x) = sin(x^2) muestra en números la competencia
entre el error de truncamiento (que baja con h) y el error de
redondeo (que sube con h) descrita en notas.md.
"""

from pathlib import Path

from fiscomp.funciones_especiales import coseno, seno
from fiscomp.precision_numerica import EPS, error_relativo

CARPETA_DATOS = Path(__file__).resolve().parent / "datos"

###############################################
# El paso h óptimo de cada método
###############################################

# Las fórmulas de notas.md piden f(x0), f''(x0) y f'''(x0), que en la
# práctica casi nunca conocemos (¡si conociéramos f''', seguramente
# también conoceríamos f'!). Por eso valen 1.0 por defecto: así se
# obtienen las fórmulas "de juguete" de notas.md, que suponen f y sus
# derivadas de orden 1.
#
# Limitación: si x0 es muy grande (|x0| > ~1e8 para adelante/atrás),
# ese h queda por debajo de la separación entre floats vecinos de x0,
# x0 + h == x0 y la derivada sale 0. No escalamos h con |x0| para
# evitarlo porque eso supone que f cambia en una escala ~|x0|, y
# falla con funciones que oscilan (sin(x) en x0=1e8 cambia en una
# escala ~1, no ~1e8). En esos casos, den h a mano.


def h_optimo_adelante(f_x0: float = 1.0, f2_x0: float = 1.0) -> float:
    """Paso h óptimo para diff_forward y diff_backward (ver notas.md):

        h_opt = sqrt(4 * EPS * |f(x0)| / |f''(x0)|)

    f_x0: valor de f(x0) (distinto de 0).
    f2_x0: valor de la segunda derivada f''(x0) (distinto de 0).
    Por defecto ambos valen 1.0 (la fórmula "de juguete").
    """
    return (4.0 * EPS * abs(f_x0) / abs(f2_x0)) ** 0.5


def h_optimo_central(f_x0: float = 1.0, f3_x0: float = 1.0) -> float:
    """Paso h óptimo para diff_central (ver notas.md):

        h_opt = (24 * EPS * |f(x0)| / |f'''(x0)|) ** (1/3)

    f_x0: valor de f(x0) (distinto de 0).
    f3_x0: valor de la tercera derivada f'''(x0) (distinto de 0).
    Por defecto ambos valen 1.0 (la fórmula "de juguete").
    """
    return (24.0 * EPS * abs(f_x0) / abs(f3_x0)) ** (1.0 / 3.0)


###############################################
# Las tres aproximaciones de diferencias finitas
###############################################

# `x0: float` y `h: float` indican que se esperan flotantes, y
# `-> float` que la función regresa un flotante. Python no lo revisa
# al correr: es documentación. `h: float = h_optimo_adelante()` hace
# que, si no se da h, se use el h óptimo.


def diff_forward(f, x0: float, h: float = h_optimo_adelante()) -> float:
    """Diferencia hacia adelante: (f(x0+h) - f(x0)) / h.

    Usa f en x0 y x0+h; error de truncamiento O(h).

    f: función de una variable, que recibe un flotante y regresa un
       flotante.
    x0: punto donde se aproxima la derivada.
    h: tamaño de paso; si no se da, se usa h_optimo_adelante().
    """
    return (f(x0 + h) - f(x0)) / h


def diff_backward(f, x0: float, h: float = h_optimo_adelante()) -> float:
    """Diferencia hacia atrás: (f(x0) - f(x0-h)) / h.

    Usa f en x0-h y x0; error de truncamiento O(h), igual que
    diff_forward pero moviéndose en la dirección opuesta (por eso
    comparten el mismo h óptimo).

    f: función de una variable, que recibe un flotante y regresa un
       flotante.
    x0: punto donde se aproxima la derivada.
    h: tamaño de paso; si no se da, se usa h_optimo_adelante().
    """
    return (f(x0) - f(x0 - h)) / h


def diff_central(f, x0: float, h: float = h_optimo_central()) -> float:
    """Diferencia central: (f(x0+h/2) - f(x0-h/2)) / h.

    Usa f en x0-h/2 y x0+h/2 (centrados en x0, separados por h, igual
    que diff_forward/diff_backward); error de truncamiento O(h^2), un
    orden mejor porque los términos pares de la serie de Taylor se
    cancelan al restar.

    f: función de una variable, que recibe un flotante y regresa un
       flotante.
    x0: punto donde se aproxima la derivada.
    h: tamaño de paso; si no se da, se usa h_optimo_central().
    """
    return (f(x0 + h / 2.0) - f(x0 - h / 2.0)) / h


###############################################
# Diferencias finitas más precisas (ver notas.md)
###############################################


def h_optimo_adelante2(f_x0: float = 1.0, f3_x0: float = 1.0) -> float:
    """Paso h óptimo para diff_forward2 (ver notas.md):

        h_opt = (48 * EPS * |f(x0)| / |f'''(x0)|) ** (1/3)

    f_x0: valor de f(x0) (distinto de 0).
    f3_x0: valor de la tercera derivada f'''(x0) (distinto de 0).
    Por defecto ambos valen 1.0 (la fórmula "de juguete").
    """
    return (48.0 * EPS * abs(f_x0) / abs(f3_x0)) ** (1.0 / 3.0)


def h_optimo_central2(f_x0: float = 1.0, f5_x0: float = 1.0) -> float:
    """Paso h óptimo para diff_central2 (ver notas.md):

        h_opt = (1120 * EPS * |f(x0)| / (9 * |f^(5)(x0)|)) ** (1/5)

    f_x0: valor de f(x0) (distinto de 0).
    f5_x0: valor de la quinta derivada f^(5)(x0) (distinto de 0).
    Por defecto ambos valen 1.0 (la fórmula "de juguete").
    """
    return (1120.0 * EPS * abs(f_x0) / (9.0 * abs(f5_x0))) ** (1.0 / 5.0)


def diff_forward2(f, x0: float, h: float = h_optimo_adelante2()) -> float:
    """Segunda diferencia hacia adelante:

        (4 f(x0+h/2) - f(x0+h) - 3 f(x0)) / h

    Usa f en x0, x0+h/2 y x0+h (tres evaluaciones, todas a la derecha
    de x0 o en x0); error de truncamiento O(h^2), como diff_central,
    pero sin necesitar puntos a la izquierda de x0.

    f: función de una variable, que recibe un flotante y regresa un
       flotante.
    x0: punto donde se aproxima la derivada.
    h: tamaño de paso; si no se da, se usa h_optimo_adelante2().
    """
    return (4.0 * f(x0 + h / 2.0) - f(x0 + h) - 3.0 * f(x0)) / h


def diff_central2(f, x0: float, h: float = h_optimo_central2()) -> float:
    """Segunda diferencia central:

        (27 f(x0+h/2) + f(x0-3h/2) - 27 f(x0-h/2) - f(x0+3h/2)) / (24 h)

    Usa f en x0 ± h/2 y x0 ± 3h/2 (cuatro evaluaciones, simétricas
    alrededor de x0); error de truncamiento O(h^4).

    f: función de una variable, que recibe un flotante y regresa un
       flotante.
    x0: punto donde se aproxima la derivada.
    h: tamaño de paso; si no se da, se usa h_optimo_central2().
    """
    return (
        27.0 * f(x0 + h / 2.0)
        + f(x0 - 3.0 * h / 2.0)
        - 27.0 * f(x0 - h / 2.0)
        - f(x0 + 3.0 * h / 2.0)
    ) / (24.0 * h)


# Respuesta del "Para pensar" de notas.md: otra fórmula O(h^4),
# combinando la diferencia central con paso h,
#     D1 = (f(x0+h/2) - f(x0-h/2)) / h  = f' + (h^2/24) f''' + (h^4/1920) f^(5) + ...
# con la diferencia central con paso 2h,
#     D2 = (f(x0+h) - f(x0-h)) / (2h)   = f' + (h^2/6)  f''' + (h^4/120)  f^(5) + ...
# El término h^2 de D2 es 4 veces el de D1, así que en (4 D1 - D2)
# se cancela, y queda 3 f' - (h^4/160) f^(5). Dividiendo entre 3:
#     f'(x0) = (8 [f(x0+h/2) - f(x0-h/2)] - [f(x0+h) - f(x0-h)]) / (6h)
#              + (h^4/480) f^(5)(x0) + ...
# (Esta idea -- combinar dos aproximaciones con distinto h para
# cancelar el término principal del error -- se llama extrapolación
# de Richardson.) El error de redondeo es (8+8+1+1) EPS |f| / (6h)
# = 3 EPS |f| / h; minimizando E = (h^4/480)|f^(5)| + 3 EPS |f| / h
# da h_opt = (360 EPS |f| / |f^(5)|)^(1/5), ~2e-3 con f y f^(5) ~ 1.


def diff_central_richardson(f, x0: float, h: float = (360.0 * EPS) ** (1.0 / 5.0)) -> float:
    """Diferencia central O(h^4) por extrapolación de Richardson:

        (8 [f(x0+h/2) - f(x0-h/2)] - [f(x0+h) - f(x0-h)]) / (6 h)

    Usa f en x0 ± h/2 y x0 ± h (cuatro evaluaciones); error de
    truncamiento O(h^4), como diff_central2.

    f: función de una variable, que recibe un flotante y regresa un
       flotante.
    x0: punto donde se aproxima la derivada.
    h: tamaño de paso; si no se da, se usa (360 * EPS) ** (1/5).
    """
    return (8.0 * (f(x0 + h / 2.0) - f(x0 - h / 2.0)) - (f(x0 + h) - f(x0 - h))) / (6.0 * h)


###############################################
# Funciones de prueba y sus derivadas exactas
###############################################


def const_5(x):
    return 5.0


def const_5_prima(x):
    return 0.0


def ident(x):
    return x


def ident_prima(x):
    return 1.0


def sqr(x):
    return x**2


def sqr_prima(x):
    return 2 * x


def sin_x2(x):
    return seno(x**2)


def sin_x2_prima(x):
    return 2 * x * coseno(x**2)


FUNCIONES_DE_PRUEBA = (
    ("f(x)=5", const_5, const_5_prima),
    ("f(x)=x", ident, ident_prima),
    ("f(x)=x^2", sqr, sqr_prima),
    ("f(x)=sin(x^2)", sin_x2, sin_x2_prima),
)

###############################################
# Comparar los tres métodos con una h fija
###############################################

x0 = 6.0
h = 0.1

print(f"Error relativo de cada método en x0={x0} con h={h}:\n")
print(f"{'función':<16}{'exacta':>14}{'err. adelante':>16}{'err. atrás':>16}{'err. central':>16}")
for nombre, f, f_prima in FUNCIONES_DE_PRUEBA:
    exacta = f_prima(x0)
    err_adelante = error_relativo(diff_forward(f, x0, h), exacta)
    err_atras = error_relativo(diff_backward(f, x0, h), exacta)
    err_central = error_relativo(diff_central(f, x0, h), exacta)
    print(f"{nombre:<16}{exacta:14.6f}{err_adelante:16.2e}{err_atras:16.2e}{err_central:16.2e}")

# Para f(x)=5 y f(x)=x, los tres métodos dan (esencialmente) error
# cero: la serie de Taylor de una función constante o lineal se corta
# exactamente en el término de la derivada que buscamos, sin nada que
# truncar (lo poco que queda es error de redondeo). Para f(x)=x^2 y
# f(x)=sin(x^2), que sí tienen derivadas de orden superior distintas
# de cero, aparece el error de truncamiento -- y es más chico en la
# columna central que en adelante/atrás, como predicen los órdenes
# O(h) vs. O(h^2). Para f(x)=x^2 la central incluso sale exacta: su
# error de truncamiento es proporcional a f''', que para x^2 vale 0.
#
# Ojo con sin(x^2) en x0=6, por dos razones:
# - f'(x) = 2x cos(x^2) oscila muy rápido ahí (el factor 2x vale 12),
#   así que h=0.1 es un paso demasiado grande.
# - seno y coseno de fiscomp suman la serie de Taylor sin reducir el
#   argumento a [-pi, pi] (ver su docstring), y en x^2 = 36 pierden
#   precisión: seno(36) se equivoca en ~0.03. Incluso la columna
#   "exacta" está mal: da -1.710, y el valor real es
#   12*cos(36) = -1.536.
print(30 * "=")

###############################################
# Los mismos métodos, sin dar h (usan su h óptimo)
###############################################

# Aquí usamos x0 = 1: en x0 = 6, seno(36) de fiscomp trae un error de
# ~0.03 (ver arriba), que dividido entre un h ~ 1e-8 se vuelve enorme.
x0 = 1.0

print(f"\nError relativo en x0={x0} con el h óptimo de cada método (sin dar h):\n")
print(f"{'función':<16}{'exacta':>14}{'err. adelante':>16}{'err. atrás':>16}{'err. central':>16}")
for nombre, f, f_prima in FUNCIONES_DE_PRUEBA:
    exacta = f_prima(x0)
    err_adelante = error_relativo(diff_forward(f, x0), exacta)
    err_atras = error_relativo(diff_backward(f, x0), exacta)
    err_central = error_relativo(diff_central(f, x0), exacta)
    print(f"{nombre:<16}{exacta:14.6f}{err_adelante:16.2e}{err_atras:16.2e}{err_central:16.2e}")

# Con el h óptimo, el error de sin(x^2) queda en ~1e-8 (adelante/
# atrás) y ~1e-10 (central): los órdenes de magnitud de los
# errores mínimos de notas.md.
print(30 * "=")

###############################################
# Barrido de h: error de truncamiento vs. error de redondeo
###############################################

# Para f(x) = sin(x^2) en x0 = 1.0, vamos reduciendo h a la mitad
# repetidamente y comparamos diff_forward y diff_central contra la
# derivada exacta. Para h grande domina el error de truncamiento (baja
# al reducir h); para h muy chica domina el error de redondeo (sube al
# reducir h, por la resta de valores casi iguales en el numerador).

x0 = 1.0
exacta = sin_x2_prima(x0)

print(f"\nBarrido de h para f(x)=sin(x^2) en x0={x0} (derivada exacta = {exacta:.12f}):\n")
print(f"{'h':>12}{'adelante':>16}{'err. adelante':>16}{'central':>16}{'err. central':>16}")

filas = []
h = 1.0
for _ in range(50):
    adelante = diff_forward(sin_x2, x0, h)
    central = diff_central(sin_x2, x0, h)
    err_adelante = error_relativo(adelante, exacta)
    err_central = error_relativo(central, exacta)
    filas.append((h, adelante, err_adelante, central, err_central))
    h *= 0.5

for h, adelante, err_adelante, central, err_central in filas:
    print(f"{h:12.2e}{adelante:16.10f}{err_adelante:16.2e}{central:16.10f}{err_central:16.2e}")

h_mejor_adelante = min(filas, key=lambda fila: fila[2])
h_mejor_central = min(filas, key=lambda fila: fila[4])
print(
    f"\nMejor h (diferencia hacia adelante): {h_mejor_adelante[0]:.2e}"
    f"  (error {h_mejor_adelante[2]:.2e})"
)
print(
    f"Mejor h (diferencia central):        {h_mejor_central[0]:.2e}"
    f"  (error {h_mejor_central[4]:.2e})"
)

# Comparación contra las h_opt "de juguete" de notas.md, que asumen f
# y sus derivadas de orden 1 (con lo que los cocientes f/f'' y f/f'''
# de las fórmulas completas se vuelven ~1): es lo que dan
# h_optimo_adelante() y h_optimo_central() sin argumentos.
print(f"\nh_opt de juguete (adelante) ~ sqrt(4*EPS)      = {h_optimo_adelante():.2e}")
print(f"h_opt de juguete (central)  ~ (24*EPS)^(1/3)    = {h_optimo_central():.2e}")

# Y con los valores reales de f, f'' y f''' de sin(x^2) en x0 = 1:
#     f''(x)  = 2 cos(x^2) - 4x^2 sin(x^2)
#     f'''(x) = -12x sin(x^2) - 8x^3 cos(x^2)
f_x0 = sin_x2(x0)
f2_x0 = 2 * coseno(x0**2) - 4 * x0**2 * seno(x0**2)
f3_x0 = -12 * x0 * seno(x0**2) - 8 * x0**3 * coseno(x0**2)
print(f"h_opt con f, f'' reales (adelante)              = {h_optimo_adelante(f_x0, f2_x0):.2e}")
print(f"h_opt con f, f''' reales (central)              = {h_optimo_central(f_x0, f3_x0):.2e}")
print(30 * "=")

# Respuesta del Ejercicio 4: ¿coinciden exactamente? No.
#
# - Adelante: la h medida (~1.5e-8) queda a un factor ~2 de la
#   fórmula (~3.0e-8), el mismo orden de magnitud.
# - Central: la h medida (~1.9e-6) es casi 10 veces más chica que la
#   de la fórmula (~1.75e-5).
#
# No tendrían que coincidir exactamente, por dos razones:
#
# 1. La fórmula supone que f y sus derivadas son de orden 1, y para
#    sin(x^2) en x0=1 no lo son:
#        f''(x)  = 2 cos(x^2) - 4x^2 sin(x^2)          -> f''(1)  ~ -2.3
#        f'''(x) = -12x sin(x^2) - 8x^3 cos(x^2)       -> f'''(1) ~ -14.4
#    Con la f'' real, la h_opt de adelante baja en un factor
#    sqrt(2.3) ~ 1.5 (a ~2e-8); con la f''' real, la central baja en
#    un factor (14.4)^(1/3) ~ 2.4 (a ~7e-6). Ambas se acercan a lo
#    medido. La central se ve más afectada porque su f''' es mucho
#    más grande que 1.
# 2. El barrido solo prueba potencias de 1/2 (h = 1, 1/2, 1/4, ...),
#    no cualquier valor de h; y cerca del mínimo el error de redondeo
#    es irregular (depende de cómo caen los bits de cada resta), así
#    que el punto exacto del mínimo "brinca" entre valores vecinos.
#    Lo que importa es el orden de magnitud y la forma de la curva,
#    no el valor exacto.

###############################################
# Guardar el barrido en un archivo, para graficarlo
###############################################

CARPETA_DATOS.mkdir(exist_ok=True)
ruta_datos = CARPETA_DATOS / "derivada_sin_x2.dat"
with open(ruta_datos, "w") as archivo:
    archivo.write("# h  diff_forward  error_forward  diff_central  error_central\n")
    for h, adelante, err_adelante, central, err_central in filas:
        archivo.write(f"{h} {adelante} {err_adelante} {central} {err_central}\n")

print(f"Barrido guardado en {ruta_datos.name}")
print(
    "Graficando log(error) contra log(h) se ven las dos pendientes de "
    "notas.md: el error de adelante cae como O(h) hasta que el "
    "redondeo lo detiene, y el de central cae más rápido, como O(h^2), "
    "antes de que el redondeo también lo alcance (más tarde, con h "
    "más chica, porque parte de un error mucho menor)."
)

###############################################
# Métodos más precisos: comparación y barrido de h
###############################################

# Mismo f(x) = sin(x^2) y mismo x0 = 1.0 que el barrido de arriba,
# ahora comparando los cinco métodos. Para las h_opt con valores
# reales hace falta también la quinta derivada:
#     f^(5)(x) = -120x cos(x^2) + 160x^3 sin(x^2) + 32x^5 cos(x^2)
f5_x0 = -120 * x0 * coseno(x0**2) + 160 * x0**3 * seno(x0**2) + 32 * x0**5 * coseno(x0**2)

METODOS = (
    # nombre, función, evaluaciones de f, h_opt con valores reales
    ("adelante", diff_forward, 2, h_optimo_adelante(f_x0, f2_x0)),
    ("central", diff_central, 2, h_optimo_central(f_x0, f3_x0)),
    ("segunda adelante", diff_forward2, 3, h_optimo_adelante2(f_x0, f3_x0)),
    ("segunda central", diff_central2, 4, h_optimo_central2(f_x0, f5_x0)),
    ("Richardson", diff_central_richardson, 4, (360.0 * EPS * abs(f_x0) / abs(f5_x0)) ** (1.0 / 5.0)),
)

print(f"\nCinco métodos para f(x)=sin(x^2) en x0={x0}, con su h por defecto (de juguete):\n")
print(f"{'método':<18}{'evals.':>7}{'aproximación':>18}{'err. relativo':>16}")
for nombre, metodo, evaluaciones, _ in METODOS:
    aproximacion = metodo(sin_x2, x0)
    print(f"{nombre:<18}{evaluaciones:>7}{aproximacion:18.12f}{error_relativo(aproximacion, exacta):16.2e}")

# Barrido de h, igual que arriba (h = 1, 1/2, 1/4, ...), guardando
# el error relativo de cada método.
filas_precisas = []
h = 1.0
for _ in range(50):
    errores = [error_relativo(metodo(sin_x2, x0, h), exacta) for _, metodo, _, _ in METODOS]
    filas_precisas.append((h, *errores))
    h *= 0.5

print(f"\n{'método':<18}{'mejor h (barrido)':>19}{'error mínimo':>14}{'h_opt (fórmula)':>17}")
for columna, (nombre, _, _, h_opt) in enumerate(METODOS, start=1):
    mejor = min(filas_precisas, key=lambda fila: fila[columna])
    print(f"{nombre:<18}{mejor[0]:19.2e}{mejor[columna]:14.2e}{h_opt:17.2e}")

# Lo que se ve:
# - Cada método de orden más alto permite una h_opt más grande y llega
#   a un error mínimo más chico: ~1e-8 (adelante, O(h)), ~1e-10/1e-11
#   (central y segunda adelante, O(h^2)) y ~1e-13/1e-14 (segunda
#   central y Richardson, O(h^4)), aunque con el h de juguete se
#   quedan en ~5e-12, porque f^(5)(1) ~ 87 está lejos de 1.
# - La segunda diferencia hacia adelante llega al mismo orden de error
#   que la central, con una evaluación más de f, pero sin usar puntos a
#   la izquierda de x0.
# - Igual que en el Ejercicio 4, la h del barrido y la de la fórmula
#   no coinciden exactamente (redondeo irregular cerca del mínimo y
#   solo potencias de 1/2), pero sí en orden de magnitud.
print(30 * "=")

ruta_datos_precisas = CARPETA_DATOS / "derivada_sin_x2_precisas.dat"
with open(ruta_datos_precisas, "w") as archivo:
    archivo.write(
        "# h  error_forward  error_central  error_forward2  error_central2  error_richardson\n"
    )
    for fila in filas_precisas:
        archivo.write(" ".join(str(valor) for valor in fila) + "\n")

print(f"Barrido de los cinco métodos guardado en {ruta_datos_precisas.name}")
