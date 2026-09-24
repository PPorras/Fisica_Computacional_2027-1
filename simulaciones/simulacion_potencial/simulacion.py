"""Cálculo del potencial electrostático, para la simulación del
arreglo de cargas de potencial_arreglo_cargas.py (misma carpeta):

    from simulacion import potencial_electrostatico
"""

from fiscomp.funciones_especiales import PI

EPSILON_0 = 8.8541878128e-12  # F/m, permitividad del vacío (CODATA)
COULOMB_K = 1.0 / (4 * PI * EPSILON_0)  # N*m^2/C^2


def potencial_electrostatico(r, cargas, posiciones, k=COULOMB_K):
    """Potencial electrostático en el punto r debido a n cargas puntuales.

    r: VectorND, punto donde se evalúa el potencial
    cargas: lista de n cargas q_i (Coulombs)
    posiciones: lista de n VectorND, la posición r_i de cada carga
    k: constante de Coulomb; por default 1/(4*pi*epsilon_0) en SI

    Por el principio de superposición, la contribución de cada carga
    es k*q_i / |r - r_i|, y el potencial total es la suma escalar de
    esas contribuciones (a diferencia del campo eléctrico, que es
    vectorial). Lanza ZeroDivisionError si r coincide exactamente con
    la posición de alguna carga, donde el potencial diverge.
    """
    if len(cargas) != len(posiciones):
        raise ValueError(
            "cargas y posiciones deben tener la misma longitud: "
            f"{len(cargas)} != {len(posiciones)}"
        )

    total = 0.0
    for q, r_i in zip(cargas, posiciones):
        distancia = (r - r_i).norm()
        if distancia == 0.0:
            raise ZeroDivisionError(
                f"el punto {r} coincide con la carga en {r_i}: "
                "el potencial diverge ahí"
            )
        total += k * q / distancia
    return total
