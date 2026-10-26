#!/usr/bin/env python3
"""Aritmética de intervalos como ejemplo de programación orientada a objetos.

Ejemplos ejecutables de lo visto en README.md: la clase Intervalo
sobrecarga los operadores aritméticos para que [a, b] + [c, d], etc.
se comporten como en las fórmulas de la aritmética de intervalos, y
al final se usa esa clase para acotar el rango de una función.
"""


class Intervalo:
    """Representa el intervalo cerrado [lo, hi] de números reales.

    Los operadores aritméticos (+, -, *, /) están sobrecargados para
    que devuelvan el intervalo que contiene todos los resultados
    posibles de operar cualquier par de reales tomados de los
    intervalos originales (ver las fórmulas en README.md).
    """

    def __init__(self, lo, hi=None):
        if hi is None:
            hi = lo
        if lo > hi:
            raise ValueError(f"lo ({lo}) no puede ser mayor que hi ({hi})")
        self.lo = float(lo)
        self.hi = float(hi)

    def __repr__(self):
        return f"Intervalo({self.lo!r}, {self.hi!r})"

    def __eq__(self, otro):
        if not isinstance(otro, Intervalo):
            return NotImplemented
        return self.lo == otro.lo and self.hi == otro.hi

    def __contains__(self, x):
        return self.lo <= x <= self.hi

    @property
    def ancho(self):
        return self.hi - self.lo

    @property
    def centro(self):
        return (self.lo + self.hi) / 2.0

    @staticmethod
    def _coaccionar(valor):
        """Convierte un int/float suelto en un Intervalo de ancho 0,
        para poder operar Intervalo con un número normal (p. ej.
        Intervalo(1, 2) + 3)."""
        if isinstance(valor, Intervalo):
            return valor
        return Intervalo(valor, valor)

    def __add__(self, otro):
        otro = self._coaccionar(otro)
        return Intervalo(self.lo + otro.lo, self.hi + otro.hi)

    __radd__ = __add__

    def __neg__(self):
        return Intervalo(-self.hi, -self.lo)

    def __sub__(self, otro):
        otro = self._coaccionar(otro)
        return Intervalo(self.lo - otro.hi, self.hi - otro.lo)

    def __rsub__(self, otro):
        return self._coaccionar(otro) - self

    def __mul__(self, otro):
        otro = self._coaccionar(otro)
        productos = (
            self.lo * otro.lo,
            self.lo * otro.hi,
            self.hi * otro.lo,
            self.hi * otro.hi,
        )
        return Intervalo(min(productos), max(productos))

    __rmul__ = __mul__

    def __truediv__(self, otro):
        otro = self._coaccionar(otro)
        if otro.lo <= 0.0 <= otro.hi:
            raise ZeroDivisionError(
                f"no se puede dividir entre un intervalo que contiene 0: {otro}"
            )
        inverso = Intervalo(1.0 / otro.hi, 1.0 / otro.lo)
        return self * inverso

    def __rtruediv__(self, otro):
        return self._coaccionar(otro) / self

    def __pow__(self, n):
        """Potencia entera no negativa. A diferencia de encadenar
        self * self * ..., no sufre el problema de dependencia porque
        se calcula directamente sobre lo y hi (ver README.md)."""
        if not isinstance(n, int) or n < 0:
            raise ValueError("solo se soporta potencia entera no negativa")
        if n == 0:
            return Intervalo(1.0, 1.0)
        candidatos = (self.lo**n, self.hi**n)
        if n % 2 == 0 and 0.0 in self:
            return Intervalo(0.0, max(candidatos))
        return Intervalo(min(candidatos), max(candidatos))


###############################################
# Operaciones básicas y contención garantizada
###############################################

X = Intervalo(1.0, 2.0)
Y = Intervalo(3.0, 5.0)

print(f"X = {X}, Y = {Y}")
print(f"X + Y = {X + Y}")
print(f"X - Y = {X - Y}")
print(f"X * Y = {X * Y}")
print(f"X / Y = {X / Y}")
print(f"X + 10 = {X + 10}")
print(30 * "=")

# Contención garantizada: para cualquier x en X y y en Y, x+y siempre
# cae dentro de X+Y (lo comprobamos muestreando, no es una prueba,
# pero ilustra la propiedad).
suma = X + Y
muestras = [
    (x, y, x + y)
    for x in (1.0, 1.3, 1.7, 2.0)
    for y in (3.0, 4.0, 4.5, 5.0)
]
todas_contenidas = all(resultado in suma for _, _, resultado in muestras)
print(f"¿Todas las sumas puntuales caen dentro de X+Y={suma}? {todas_contenidas}")
print(30 * "=")

###############################################
# El problema de dependencia: X - X != [0, 0]
###############################################

# La aritmética de intervalos trata cada aparición de X como si fuera
# independiente, así que resta el peor caso contra el peor caso en vez
# de "darse cuenta" de que es la misma variable.
print(f"X - X = {X - X}  (matemáticamente X - X = 0 para cualquier x)")

# Consecuencia práctica: la misma función da cotas distintas según
# cómo esté escrita la expresión. Con Z que incluye negativos:
Z = Intervalo(-1.0, 2.0)
print(f"\nZ = {Z}")
print(f"Z * Z   = {Z * Z}   (usa __mul__, trata los dos Z como independientes)")
print(f"Z ** 2  = {Z ** 2}   (usa __pow__, calculado directo -> cota más angosta)")
print(30 * "=")

###############################################
# Acotar el rango de una función evaluándola sobre un intervalo
###############################################

# Si f solo usa +, -, *, /, ** con literales, evaluarla con un
# Intervalo en vez de un float acota su rango en ese intervalo, sin
# necesidad de muestrear ni de calcular derivadas.


def f(x):
    return x**2 - 2.0 * x


dominio = Intervalo(0.0, 3.0)
rango_acotado = f(dominio)
print(f"f(x) = x^2 - 2x, evaluada sobre {dominio}: {rango_acotado}")

# Comparamos contra el rango real, obtenido por muestreo denso (solo
# posible aquí porque conocemos f de antemano; el punto de la
# aritmética de intervalos es no depender de esto).
muestras_f = [f(0.0 + 0.001 * i) for i in range(3001)]
rango_real = Intervalo(min(muestras_f), max(muestras_f))
print(f"Rango real (muestreado):                  {rango_real}")
print(
    "La cota por intervalos siempre contiene al rango real, pero aquí "
    "sale bastante más ancha ([-6, 9] vs. [-1, 3]): aunque x**2 se "
    "calcula sin problema (usa __pow__), la resta x**2 - 2*x vuelve a "
    "sufrir el problema de dependencia, porque __sub__ no sabe que la "
    "x de ambos términos es la misma variable."
)
