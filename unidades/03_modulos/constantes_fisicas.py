#!/usr/bin/env python3
"""Módulo de ejemplo: algunas constantes físicas y una función simple.

Se usa desde modulos.py para mostrar cómo se importa un módulo propio
(en vez de uno de la librería estándar o del paquete fiscomp).
"""

GRAVEDAD = 9.81             # m/s^2, en la Tierra
VELOCIDAD_LUZ = 299_792.458 # m/s
PLANCK = 6.626_070_15e-34   # J*s


def energia_potencial(masa, altura, gravedad=GRAVEDAD):
    """Calcula la energía potencial gravitacional (J)."""
    return masa * gravedad * altura


if __name__ == "__main__":
    # Este bloque solo corre si el módulo se ejecuta directamente
    # (`python3 constantes_fisicas.py`), no cuando alguien más lo importa.
    print(f"GRAVEDAD = {GRAVEDAD} m/s^2")
    print(f"energia_potencial(2 kg, 10 m) = {energia_potencial(2.0, 10.0)} J")
