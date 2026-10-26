#!/usr/bin/env python3
"""Diferenciación numérica con diferencias finitas -- punto de partida.

Ver README.md para la teoría (cada método se deriva de la serie de
Taylor de f alrededor de x0, con un orden de error distinto) y
practica.md para el enunciado completo de lo que falta completar
aquí.

Ya está implementada diff_forward (diferencia hacia adelante), junto
con las funciones de prueba y un ejemplo de uso. El resto -- las
diferencias hacia atrás y central, y el barrido de h que compara el
error de truncamiento contra el de redondeo -- es la práctica.
"""

from fiscomp.funciones_especiales import seno

###############################################
# Diferencias finitas
###############################################


def diff_forward(f, x0, h):
    """Diferencia hacia adelante: (f(x0+h) - f(x0)) / h.

    Usa f en x0 y x0+h; error de truncamiento O(h) (ver README.md).
    """
    return (f(x0 + h) - f(x0)) / h


# TODO (práctica, Ejercicio 1): diff_backward(f, x0, h)
#   diferencia hacia atrás, (f(x0) - f(x0-h)) / h, error O(h).

# TODO (práctica, Ejercicio 1): diff_central(f, x0, h)
#   diferencia central, (f(x0+h/2) - f(x0-h/2)) / h, error O(h^2).


###############################################
# Funciones de prueba
###############################################


def const_5(x):
    return 5.0


def ident(x):
    return x


def sqr(x):
    return x**2


def sin_x2(x):
    return seno(x**2)


###############################################
# Ejemplo de uso de diff_forward
###############################################

x0 = 6.0
h = 0.1

print(f"La derivada de f(x)=5, en x0={x0} con h={h}, es {diff_forward(const_5, x0, h)}")
print(f"La derivada de f(x)=x, en x0={x0} con h={h}, es {diff_forward(ident, x0, h)}")
print(f"La derivada de f(x)=x^2, en x0={x0} con h={h}, es {diff_forward(sqr, x0, h)}")
print(f"La derivada de f(x)=sin(x^2), en x0={x0} con h={h}, es {diff_forward(sin_x2, x0, h)}")

# A partir de aquí empieza la práctica (ver practica.md):
#
#   Ejercicio 1: diff_backward y diff_central.
#   Ejercicio 2: derivadas exactas de las funciones de prueba, y una
#                tabla que compare los tres métodos contra ellas.
#   Ejercicio 3: barrido de h para f(x)=sin(x^2), guardando
#                h, diff_forward, error_forward, diff_central,
#                error_central en datos/derivada_sin_x2.dat.
#   Ejercicio 4: encontrar el h que minimiza cada error y comparalo
#                contra la formula de h_opt de README.md.
#   Ejercicio 5 (opcional): graficar con graficar_derivada.gp o
#                graficar_derivada.py.
