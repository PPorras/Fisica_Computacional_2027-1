#!/usr/bin/env python3
"""Herencia: una clase hija de Punto.

Ver poo.md, sección "Herencia". La clase Punto viene de punto.py (en
esta misma carpeta). Este script:

1. Define PuntoConMasa(Punto), una partícula puntual: un punto que
   además tiene masa. Agrega un atributo (con super().__init__), un
   método nuevo, y sobreescribe __repr__.
2. Muestra qué se hereda tal cual, incluido un método heredado cuyo
   resultado ya no es de la clase hija (p + q pierde la masa).
3. Muestra la cadena de búsqueda de métodos (__mro__), isinstance e
   issubclass.
4. Muestra el polimorfismo: el mismo print() usa el __repr__ de la
   clase de cada objeto.
5. Muestra la jerarquía de las excepciones de Python, que ya usaron
   en la unidad 05.
"""

from punto import Punto

###############################################
# 1. La clase hija
###############################################


class PuntoConMasa(Punto):
    """Partícula puntual en el plano: un Punto con masa.

    Hereda de Punto: x, y, distancia_al_origen, ==, +, -, * y /.
    Agrega el atributo masa y el método momento_de_inercia, y
    sobreescribe __repr__ para que muestre la masa.
    """

    def __init__(self, x, y, masa):
        # Punto.__init__ guarda x y y; aquí solo agregamos lo nuevo.
        super().__init__(x, y)
        if masa <= 0:
            raise ValueError("La masa debe ser positiva.")
        self.masa = masa

    def momento_de_inercia(self):
        """Momento de inercia respecto al origen, m r^2. Usa
        distancia_al_origen, que se heredó de Punto."""
        return self.masa * self.distancia_al_origen**2

    def __repr__(self):
        # Sobreescribe el __repr__ de Punto, que no mostraría la masa.
        return f"PuntoConMasa({self.x!r}, {self.y!r}, masa={self.masa!r})"


print("--- 1. Una clase hija ---")
p = PuntoConMasa(3, 4, masa=2.0)
print(f"p = {p}")
print(f"p.x, p.y = {p.x}, {p.y}   (los guardó Punto.__init__, vía super())")
print(f"p.masa = {p.masa}   (atributo nuevo)")
print(f"p.distancia_al_origen = {p.distancia_al_origen}   (heredado de Punto)")
print(f"p.momento_de_inercia() = {p.momento_de_inercia()}   (método nuevo: 2 * 5^2)")

try:
    PuntoConMasa(0, 0, masa=-1.0)
except ValueError as error:
    print(f"PuntoConMasa(0, 0, masa=-1.0) -> ValueError: {error}")


class PuntoSinSuper(Punto):
    """Igual que PuntoConMasa pero sin llamar a super().__init__."""

    def __init__(self, x, y, masa):
        self.masa = masa


print("\nSi se olvida super().__init__, Punto.__init__ nunca se ejecuta:")
roto = PuntoSinSuper(3, 4, masa=2.0)
try:
    print(roto.x)
except AttributeError as error:
    print(f"  roto.x -> AttributeError: {error}")

###############################################
# 2. Lo que se hereda tal cual
###############################################

print("\n--- 2. Métodos heredados ---")
q = PuntoConMasa(1, 1, masa=5.0)
print(f"q = {q}")
print(f"p == PuntoConMasa(3, 4, masa=2.0) -> {p == PuntoConMasa(3, 4, masa=2.0)}   (__eq__ de Punto)")
print(f"p == PuntoConMasa(3, 4, masa=9.0) -> {p == PuntoConMasa(3, 4, masa=9.0)}   (¡no compara la masa!)")
suma = p + q
print(f"p + q = {suma}, de tipo {type(suma).__name__}")
print("  Punto.__add__ construye un Punto(...), así que la masa se pierde.")

###############################################
# 3. Búsqueda de métodos, isinstance e issubclass
###############################################

print("\n--- 3. Búsqueda de métodos ---")
print("Orden en que Python busca un método de PuntoConMasa:")
for clase in PuntoConMasa.__mro__:
    print(f"  {clase.__name__}")
print(f"isinstance(p, PuntoConMasa) -> {isinstance(p, PuntoConMasa)}")
print(f"isinstance(p, Punto)        -> {isinstance(p, Punto)}   (un PuntoConMasa también es un Punto)")
print(f"isinstance(Punto(0, 0), PuntoConMasa) -> {isinstance(Punto(0, 0), PuntoConMasa)}   (al revés no)")
print(f"issubclass(PuntoConMasa, Punto) -> {issubclass(PuntoConMasa, Punto)}")
print(f"issubclass(Punto, object)       -> {issubclass(Punto, object)}   (todas las clases heredan de object)")

###############################################
# 4. Polimorfismo
###############################################

print("\n--- 4. Polimorfismo ---")
objetos = [Punto(1, 2), PuntoConMasa(3, 4, masa=2.0), Punto(0, 5)]
for objeto in objetos:
    # La misma línea de código; cada objeto usa el __repr__ de su clase.
    print(f"  {objeto!r:35} distancia al origen = {objeto.distancia_al_origen}")

###############################################
# 5. Las excepciones también forman una jerarquía
###############################################

print("\n--- 5. Jerarquía de excepciones ---")
print("Cadena de ZeroDivisionError:", " -> ".join(c.__name__ for c in ZeroDivisionError.__mro__))
print(f"issubclass(ZeroDivisionError, ArithmeticError) -> {issubclass(ZeroDivisionError, ArithmeticError)}")
print(f"issubclass(ValueError, ArithmeticError)        -> {issubclass(ValueError, ArithmeticError)}")
# "except ArithmeticError" atrapa a ArithmeticError y a todas sus hijas.
try:
    1 / 0
except ArithmeticError as error:
    print(f"  1 / 0: lo atrapó 'except ArithmeticError' ({type(error).__name__})")

try:
    int("hola")
except ArithmeticError as error:
    print(f"  int('hola'): lo atrapó 'except ArithmeticError' ({type(error).__name__})")
except Exception as error:
    print(f"  int('hola'): no es ArithmeticError; lo atrapó 'except Exception' ({type(error).__name__})")
