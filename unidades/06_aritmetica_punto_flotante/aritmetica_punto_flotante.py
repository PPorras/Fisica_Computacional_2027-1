#!/usr/bin/env python3
"""Underflow, overflow y cancelación catastrófica en punto flotante.

Ejemplos ejecutables de lo visto en README.md: qué pasa quando un
cálculo se sale del rango representable (UFL, OFL) y qué pasa cuando
se restan dos números casi iguales (cancelación catastrófica).
"""

import math
import sys

from fiscomp.precision_numerica import EPS, error_relativo

###############################################
# Underflow (UFL)
###############################################

# El menor positivo *normalizado* es sys.float_info.min = 2**-1022;
# coincide con UFL = beta**L de las notas (beta=2, L=-1022 en doble
# precisión IEEE 754).
print(f"UFL (menor normalizado) = sys.float_info.min = {sys.float_info.min:.6e}")

# Pero Python (como IEEE 754) sí representa números más chicos que
# UFL: los subnormales (mantisa sin normalizar, van perdiendo dígitos
# de precisión conforme se acercan a cero). El menor positivo
# representable en total es 2**-1074, no 2**-1022.
print(f"Menor subnormal representable = 2**-1074 = {2.0**-1074:.6e}")

# Dividir repetidamente entre 2 ilustra la transición: primero se
# entra a la zona subnormal (se van perdiendo dígitos de precisión),
# y al final el resultado colapsa a 0.0 sin lanzar ningún error.
x = 1.0
pasos = 0
while x > 0.0:
    x_anterior = x
    x /= 2.0
    pasos += 1
print(f"Underflow tras {pasos} divisiones entre 2: {x_anterior:.6e} -> {x}")

###############################################
# Overflow (OFL)
###############################################

# OFL = beta**(U+1) * (1 - beta**-t), el mayor float representable.
print(f"OFL = sys.float_info.max = {sys.float_info.max:.6e}")

# Multiplicar repetidamente por 10 rebasa OFL y el resultado se
# convierte, silenciosamente, en infinito (float('inf')) — a
# diferencia de un int de Python, que crece sin límite.
y = 1.0
pasos = 0
while y < math.inf:
    y_anterior = y
    y *= 10.0
    pasos += 1
print(f"Overflow tras {pasos} multiplicaciones por 10: {y_anterior:.3e} -> {y}")

# No todas las operaciones se quedan calladas al desbordarse: la
# multiplicación (*) regresa inf sin quejarse (como en el ciclo de
# arriba), pero math.exp() y hasta el propio operador ** sí lanzan
# OverflowError en vez de regresar inf silenciosamente.
try:
    math.exp(1000)
except OverflowError as error:
    print(f"math.exp(1000) lanza OverflowError: {error}")
try:
    10.0**1000
except OverflowError as error:
    print(f"10.0**1000 (con **, no con *) también lanza OverflowError: {error}")

###############################################
# Cancelación catastrófica: raíces de una cuadrática
###############################################

# a*x^2 + b*x + c = 0, con b^2 >> 4*a*c: una raíz es del orden de b/a
# (grande) y la otra del orden de c/b (chica). La fórmula general
# resta dos números casi iguales (-b y sqrt(b^2 - 4ac), casi iguales
# porque 4ac es insignificante frente a b^2) para calcular la raíz
# chica, y ahí se cancelan casi todos los dígitos significativos.


def raices_formula_directa(a, b, c):
    """Fórmula general, tal cual, sin cuidar la cancelación."""
    discriminante = b**2 - 4 * a * c
    raiz_disc = math.sqrt(discriminante)
    x1 = (-b + raiz_disc) / (2 * a)
    x2 = (-b - raiz_disc) / (2 * a)
    return x1, x2


def raices_formula_estable(a, b, c):
    """Evita la resta cancelada: primero calcula la raíz "grande"
    (sumando cantidades del mismo signo, sin cancelación), y obtiene
    la raíz "chica" a partir de la identidad x1 * x2 = c / a."""
    discriminante = b**2 - 4 * a * c
    raiz_disc = math.sqrt(discriminante)
    if b >= 0:
        x_grande = (-b - raiz_disc) / (2 * a)
    else:
        x_grande = (-b + raiz_disc) / (2 * a)
    x_chica = c / (a * x_grande)
    return x_grande, x_chica


a, b, c = 1.0, -100000.0, 1.0   # raíces exactas ~= 99999.99999 y ~= 1.00000e-5

x1_directa, x2_directa = raices_formula_directa(a, b, c)
x_grande, x_chica = raices_formula_estable(a, b, c)

print(f"Fórmula directa:  x1={x1_directa!r}, x2={x2_directa!r}")
print(f"Fórmula estable:  x_grande={x_grande!r}, x_chica={x_chica!r}")
print(
    "Error relativo de la raíz chica "
    f"(usando la estable como referencia): {error_relativo(x2_directa, x_chica):.2e}"
)

###############################################
# Cancelación catastrófica: 1 - cos(x) para x chica
###############################################

# 1 - cos(x) resta dos números que, para x chica, son casi iguales
# (cos(x) ~= 1). La identidad trigonométrica 2*sin(x/2)**2 calcula lo
# mismo evitando esa resta.


def uno_menos_coseno_directo(x):
    return 1.0 - math.cos(x)


def uno_menos_coseno_estable(x):
    return 2.0 * math.sin(x / 2.0) ** 2


print(f"{'x':>10}  {'directo':>14}  {'estable':>14}  error_relativo")
for x in (1e-1, 1e-4, 1e-7, 1e-8, 1e-9):
    directo = uno_menos_coseno_directo(x)
    estable = uno_menos_coseno_estable(x)
    print(f"{x:10.1e}  {directo:14.6e}  {estable:14.6e}  {error_relativo(directo, estable):.2e}")

# A partir de x ~ 1e-8, x**2/2 (el término dominante de 1 - cos(x)) es
# menor que EPS: cos(x) se redondea exactamente a 1.0, y la fórmula
# directa da 0.0 -- cancelación total, se pierden TODOS los dígitos
# significativos. Ese umbral es justo sqrt(2 * EPS):
print(f"\nsqrt(2 * EPS) = {math.sqrt(2 * EPS):.3e}  (umbral donde la fórmula directa colapsa a 0.0)")

###############################################
# Error de redondeo acumulado en sumas largas
###############################################

# Sumar 0.1 muchas veces "a mano" (con un acumulador) no da lo mismo
# que N * 0.1: cada += comete un redondeo, y con millones de sumas
# esos redondeos se acumulan. math.fsum() usa un algoritmo de suma
# exacta (con precisión extendida) y sirve aquí como referencia.

N = 10_000_000
terminos = [0.1] * N

suma_naive = 0.0
for termino in terminos:
    suma_naive += termino

referencia = math.fsum(terminos)

print(f"\nSuma naive de {N} veces 0.1: {suma_naive!r}")
print(f"N * 0.1:                     {N * 0.1!r}")
print(f"math.fsum (referencia):      {referencia!r}")
print(f"error de la suma naive:      {abs(suma_naive - referencia):.3e}")

# El orden de la suma también importa: comparamos sumar la serie
# armónica (1 + 1/2 + 1/3 + ... + 1/N) de mayor a menor término contra
# sumarla de menor a mayor.

terminos_de_mayor_a_menor = [1.0 / n for n in range(1, N // 10 + 1)]
terminos_de_menor_a_mayor = list(reversed(terminos_de_mayor_a_menor))


def suma_acumulada(terminos):
    suma = 0.0
    for termino in terminos:
        suma += termino
    return suma


suma_mayor_a_menor = suma_acumulada(terminos_de_mayor_a_menor)
suma_menor_a_mayor = suma_acumulada(terminos_de_menor_a_mayor)
referencia_armonica = math.fsum(terminos_de_mayor_a_menor)

print(f"\nSerie armónica (1 a {N // 10}), sumada de mayor a menor: {suma_mayor_a_menor!r}")
print(f"                              sumada de menor a mayor: {suma_menor_a_mayor!r}")
print(
    f"error (mayor a menor) = {abs(suma_mayor_a_menor - referencia_armonica):.3e}  "
    f"vs.  error (menor a mayor) = {abs(suma_menor_a_mayor - referencia_armonica):.3e}"
)


# Suma compensada (algoritmo de Kahan): en vez de cambiar el orden,
# se lleva un "recordatorio" (compensacion) de los dígitos que se
# perdieron en la suma anterior, y se restan del siguiente término
# antes de sumarlo -- así ese error no se pierde, se reintroduce.
def suma_kahan(terminos):
    suma = 0.0
    compensacion = 0.0
    for termino in terminos:
        y = termino - compensacion
        suma_tentativa = suma + y
        compensacion = (suma_tentativa - suma) - y
        suma = suma_tentativa
    return suma


suma_compensada = suma_kahan(terminos)
print(f"\nSuma de Kahan de {N} veces 0.1: {suma_compensada!r}")
print(f"error de la suma de Kahan:      {abs(suma_compensada - referencia):.3e}")

###############################################
# Ejercicio: la derivada numérica y el error de redondeo
###############################################

# La definición de derivada sugiere aproximar f'(x) con un cociente
# de diferencias: (f(x+h) - f(x)) / h. En teoría, entre más chica h,
# mejor la aproximación -- pero f(x+h) y f(x) son casi iguales para h
# chica, así que f(x+h) - f(x) sufre cancelación catastrófica, y el
# resultado se divide entre una h también chica (que amplifica el
# error). Este ejercicio muestra el resultado de esa pelea entre dos
# fuentes de error opuestas.


def derivada_numerica(f, x, h):
    return (f(x + h) - f(x)) / h


def f(x):
    return math.sin(x)


def f_prima_exacta(x):
    return math.cos(x)


x0 = 1.0
print(f"\nf'(x) exacta en x={x0}: {f_prima_exacta(x0):.10f}")
print(f"{'h':>10}  {'derivada_numerica':>18}  error_relativo")
for exponente in range(1, 17):
    h = 10.0**-exponente
    aproximada = derivada_numerica(f, x0, h)
    error = error_relativo(aproximada, f_prima_exacta(x0))
    print(f"{h:10.0e}  {aproximada:18.12f}  {error:.2e}")

# Preguntas para resolver viendo la tabla anterior (no hay que
# entregar nada, es para discutir en clase):
#
# 1. ¿Para qué valor de h el error relativo es mínimo? ¿Por qué no
#    sigue bajando si h sigue disminuyendo?
# 2. Para h grande (10**-1, 10**-2, ...) el error también es
#    relativamente alto: ¿ese error es por redondeo o por la propia
#    aproximación de la derivada (truncar la definición de límite a
#    una h finita)? (Pista: no tiene nada que ver con EPS.)
# 3. h más chica que EPS (~2.22e-16): x + h ya ni siquiera cambia el
#    valor de x en la aritmética de punto flotante (fl(x + h) == x).
#    ¿Qué le pasa a derivada_numerica(f, x0, h) en ese caso?
