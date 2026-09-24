#!/usr/bin/env python3
"""Reimplementación propia de funciones matemáticas elementales.

La idea es construir, sin usar el módulo `math` de la librería
estándar, aproximaciones numéricas de funciones como:

- factorial(n)     -- iterativa
- seno(x)          -- serie de Taylor, con reducción de rango
- coseno(x)        -- serie de Taylor, con reducción de rango
- tangente(x)      -- seno(x) / coseno(x)
- secante(x)       -- 1 / coseno(x)
- cosecante(x)     -- 1 / seno(x)
- cotangente(x)    -- coseno(x) / seno(x)
- exponencial(x)   -- serie de Taylor (para x < 0, usa 1 / e^|x|)
- ln(x)            -- serie de ln((1+y)/(1-y)), y = (x-1)/(x+1)
- raiz_cuadrada(x) -- método de Newton-Raphson
- arcotangente(x)  -- serie de Taylor, para |x| <= 1
- calcular_pi()    -- fórmula de Machin, a partir de arcotangente

PI (mayúsculas, como constante que es) guarda calcular_pi() ya
evaluada una sola vez al importar el módulo, igual que EPS guarda
epsilon_maquina() en fiscomp.precision_numerica.

Las funciones basadas en series (seno, coseno, exponencial, ln) usan
EPS (fiscomp.precision_numerica) como criterio de convergencia: se
suman términos mientras el siguiente término siga siendo mayor o
igual que el épsilon de la máquina, y se corta la suma en cuanto deja
de aportar precisión adicional. raiz_cuadrada usa el mismo EPS, pero
como criterio de paro de las iteraciones de Newton-Raphson (se detiene
cuando dos aproximaciones sucesivas ya casi no cambian).
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


# Mayor x tal que e^x todavía cabe en un float de doble precisión
# (e^709.78 ~ 1.8e308, el float más grande).
MAXIMO_EXPONENTE = 709.78


def reducir_angulo(x):
    """Lleva el ángulo x (en radianes) al intervalo [-pi, pi],
    restándole un múltiplo entero de 2 pi (seno y coseno son
    periódicos, así que no cambian):

        x_reducido = x - 2 pi * round(x / (2 pi))

    La reducción no es perfecta: PI (calculada con la fórmula de
    Machin) difiere de pi en ~1e-15, y ese error se multiplica por el
    número de vueltas que se restan. Así que el error absoluto crece
    con |x| (~1e-15 para x ~ 1, ~1e-11 para x ~ 1e5), aunque es muchísimo
    menor que sin reducir (seno(36) se equivocaba en ~0.03).
    """
    return x - 2 * PI * round(x / (2 * PI))


def seno(x, precision=EPS):
    """Aproxima sin(x) con la serie de Taylor alrededor de 0:

        sin(x) = suma_{k=0}^inf (-1)^k * x^(2k+1) / (2k+1)!

    Se suman términos mientras sigan siendo mayores o iguales que
    `precision` (por default, el épsilon de la máquina); en cuanto un
    término es más chico, ya no cambia el resultado y se detiene la suma.

    Antes de sumar, x se lleva a [-pi, pi] con reducir_angulo(): para
    |x| grande, los términos de la serie serían enormes y de signos
    alternados, y al sumarlos se perdería casi toda la precisión
    (cancelación catastrófica, unidad 06).
    """
    x = reducir_angulo(x)
    suma = 0.0
    k = 0
    while True:
        termino = (-1) ** k * x ** (2 * k + 1) / factorial(2 * k + 1)
        if abs(termino) < precision:
            break
        suma += termino
        k += 1
    return suma


def coseno(x, precision=EPS):
    """Aproxima cos(x) con la serie de Taylor alrededor de 0:

        cos(x) = suma_{k=0}^inf (-1)^k * x^(2k) / (2k)!

    Mismo criterio de corte y misma reducción de rango que `seno()`.
    """
    x = reducir_angulo(x)
    suma = 0.0
    k = 0
    while True:
        termino = (-1) ** k * x ** (2 * k) / factorial(2 * k)
        if abs(termino) < precision:
            break
        suma += termino
        k += 1
    return suma


def tangente(x, precision=EPS):
    """Aproxima tan(x) = sin(x) / cos(x), reutilizando seno() y coseno().

    Nota: no hay protección especial cerca de x = pi/2 + n*pi (donde
    cos(x) = 0 y tan(x) diverge); ahí `coseno(x)` da un valor cercano
    a 0 pero no exactamente 0, así que el resultado es un número muy
    grande en vez de un error.
    """
    return seno(x, precision) / coseno(x, precision)


def secante(x, precision=EPS):
    """Aproxima sec(x) = 1 / cos(x), reutilizando coseno()."""
    return 1.0 / coseno(x, precision)


def cosecante(x, precision=EPS):
    """Aproxima csc(x) = 1 / sin(x), reutilizando seno()."""
    return 1.0 / seno(x, precision)


def cotangente(x, precision=EPS):
    """Aproxima cot(x) = cos(x) / sin(x), reutilizando seno() y coseno()."""
    return coseno(x, precision) / seno(x, precision)


def exponencial(x, precision=EPS):
    """Aproxima e^x con la serie de Taylor alrededor de 0:

        e^x = suma_{k=0}^inf x^k / k!

    Para x > 0 todos los términos son positivos y no hay cancelación.
    Para x < 0 los términos alternan de signo y, si |x| es grande, se
    pierde casi toda la precisión (por ejemplo, con x = -50 la suma
    daba ~ -7000 en vez de ~2e-22). Por eso, para x < 0 se usa
    e^x = 1 / e^|x|, que solo suma términos positivos.

    Cada término se obtiene del anterior, termino_k = termino_(k-1) * (x / k),
    en vez de calcular x**k / k! directo: para x grande, x**k se
    desborda (no cabe en un float) mucho antes que el cociente.

    Lanza OverflowError si x > 709.78, porque e^x ya no cabe en un
    float; para x < -709.78, regresa 0.0.
    """
    if x > MAXIMO_EXPONENTE:
        raise OverflowError(f"exponencial({x}) no cabe en un float (x > {MAXIMO_EXPONENTE})")
    if x < 0:
        if -x > MAXIMO_EXPONENTE:
            return 0.0
        return 1.0 / exponencial(-x, precision)

    suma = 0.0
    termino = 1.0  # el término k = 0, x^0 / 0!
    k = 0
    while termino >= precision:
        suma += termino
        k += 1
        # (x / k) primero: termino * x podría desbordarse aunque
        # termino * x / k sí quepa en un float.
        termino = termino * (x / k)
    return suma


def ln(x, precision=EPS):
    """Aproxima ln(x), para x > 0, con la serie:

        ln(x) = 2 * suma_{k=0}^inf y^(2k+1) / (2k+1),   y = (x-1)/(x+1)

    A diferencia de la serie de Taylor de ln(x) alrededor de x = 1
    (que solo converge para 0 < x <= 2), esta serie converge para
    cualquier x > 0: entre más lejos esté x de 1, más cerca de 1 (o
    de -1) está y, y más términos hacen falta para que el término
    caiga por debajo de `precision`.
    """
    if x <= 0:
        raise ValueError("ln(x) solo está definido para x > 0")

    y = (x - 1) / (x + 1)
    suma = 0.0
    k = 0
    while True:
        termino = y ** (2 * k + 1) / (2 * k + 1)
        if abs(termino) < precision:
            break
        suma += termino
        k += 1
    return 2 * suma


def raiz_cuadrada(x, precision=EPS):
    """Aproxima sqrt(x), para x >= 0, con el método de Newton-Raphson:

        y_(n+1) = (y_n + x / y_n) / 2

    Partiendo de y_0 = x, cada iteración aproximadamente duplica el
    número de dígitos correctos; se detiene en cuanto dos
    aproximaciones sucesivas difieren en menos que `precision`.
    """
    if x < 0:
        raise ValueError("raiz_cuadrada(x) solo está definida para x >= 0")
    if x == 0:
        return 0.0

    aproximacion = x
    while True:
        siguiente = (aproximacion + x / aproximacion) / 2
        if abs(siguiente - aproximacion) < precision:
            return siguiente
        aproximacion = siguiente


def arcotangente(x, precision=EPS):
    """Aproxima arctan(x), para |x| <= 1, con la serie de Taylor:

        arctan(x) = suma_{k=0}^inf (-1)^k * x^(2k+1) / (2k+1)

    Mismo criterio de corte que seno()/coseno(). Solo se usa aquí con
    argumentos pequeños (1/5, 1/239 en calcular_pi()), donde converge
    rápido; para x cercano a 1 haría falta cada vez más términos.
    """
    if abs(x) > 1:
        raise ValueError("arcotangente(x) solo converge para |x| <= 1")

    suma = 0.0
    k = 0
    while True:
        termino = (-1) ** k * x ** (2 * k + 1) / (2 * k + 1)
        if abs(termino) < precision:
            break
        suma += termino
        k += 1
    return suma


def calcular_pi(precision=EPS):
    """Aproxima pi con la fórmula de Machin:

        pi/4 = 4*arctan(1/5) - arctan(1/239)

    Machin evalúa arcotangente() en argumentos pequeños (1/5, 1/239),
    donde la serie de Taylor converge en pocas decenas de términos. En
    contraste, la serie de Leibniz pi/4 = arctan(1) = 1 - 1/3 + 1/5 -
    ... converge tan lento (el término k-ésimo es ~1/k) que llegar al
    épsilon de la máquina tomaría miles de millones de términos.
    """
    return 16 * arcotangente(1.0 / 5.0, precision) - 4 * arcotangente(
        1.0 / 239.0, precision
    )


PI = calcular_pi()


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

    print()
    for x in (0.0, 0.5, 1.0, -1.0, 2.0, math.pi):
        aproximado = coseno(x)
        exacto = math.cos(x)
        print(
            f"coseno({x:.4f}) = {aproximado:.12f}  "
            f"math.cos = {exacto:.12f}  "
            f"error_relativo = {error_relativo(aproximado, exacto):.2e}"
        )

    print()
    for x in (0.0, 0.5, 1.0, -1.0, math.pi / 4):
        aproximado = tangente(x)
        exacto = math.tan(x)
        print(
            f"tangente({x:.4f}) = {aproximado:.12f}  "
            f"math.tan = {exacto:.12f}  "
            f"error_relativo = {error_relativo(aproximado, exacto):.2e}"
        )

    print()
    for x in (1.0, -1.0, 5.0, -10.0, 2.5):
        aproximado = exponencial(x)
        exacto = math.exp(x)
        print(
            f"exponencial({x:.4f}) = {aproximado:.12f}  "
            f"math.exp = {exacto:.12f}  "
            f"error_relativo = {error_relativo(aproximado, exacto):.2e}"
        )

    print()
    for x in (1.0, 0.1, 0.5, 2.0, 5.0, 100.0):
        aproximado = ln(x)
        exacto = math.log(x)
        print(
            f"ln({x:.4f}) = {aproximado:.12f}  "
            f"math.log = {exacto:.12f}  "
            f"error_relativo = {error_relativo(aproximado, exacto):.2e}"
        )

    print()
    for x in (0.0, 1.0, 2.0, 10.0, 0.5, 1e10):
        aproximado = raiz_cuadrada(x)
        exacto = math.sqrt(x)
        print(
            f"raiz_cuadrada({x:.4f}) = {aproximado:.12f}  "
            f"math.sqrt = {exacto:.12f}  "
            f"error_relativo = {error_relativo(aproximado, exacto):.2e}"
        )

    print()
    for x in (0.0, 1.0 / 5.0, 1.0 / 239.0, -0.5):
        aproximado = arcotangente(x)
        exacto = math.atan(x)
        print(
            f"arcotangente({x:.4f}) = {aproximado:.12f}  "
            f"math.atan = {exacto:.12f}  "
            f"error_relativo = {error_relativo(aproximado, exacto):.2e}"
        )

    print()
    print(f"PI = {PI:.15f}")
    print(f"math.pi = {math.pi:.15f}")
    print(f"error_relativo = {error_relativo(PI, math.pi):.2e}")
