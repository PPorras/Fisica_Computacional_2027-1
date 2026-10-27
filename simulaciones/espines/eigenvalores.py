#!/usr/bin/env python3
"""Eigenvalores de los hamiltonianos de espín (ver README.md).

1. El método QR (metodo_qr de fiscomp/algebra_lineal.py), con un
   desplazamiento que deja todos los eigenvalores positivos: así QR
   no se atora con eigenvalores de igual valor absoluto y signo
   contrario, ni divide entre cero si H es singular.
2. Para comprobar: multiplicidades (el residuo de un eigenpar es
   residuo_eigen, de fiscomp/algebra_lineal.py).
3. La primera versión del proyecto, antes de tener QR: las raíces de
   det(H - lambda I), calculado con LU, buscando cambios de signo en
   una malla. Se queda para comparar: no ve los eigenvalores de
   multiplicidad par ni los que caen juntos en una celda.
4. Para dos espines, los eigenvalores exactos con el bloque 2x2.
"""

from fiscomp.algebra_lineal import (
    determinante,
    menos_lambda_identidad,
    metodo_qr,
    norma_infinito,
)
from fiscomp.funciones_especiales import raiz_cuadrada


###############################################
# 1. El método QR, desplazado
###############################################


def eigenvalores_qr(H, kmax: int = 500):
    """Todos los eigenvalores de H (real y simétrica) con el método QR.

    metodo_qr solo funciona si los eigenvalores tienen valores absolutos
    distintos, y divide entre cero si H es singular. Los hamiltonianos
    de espín tienen las dos cosas: sin campo hay eigenvalores +-lambda,
    y sin interacción puede haber un cero. Por eso se aplica a

        H - s I,   con s = -(||H||_inf + 0.1),

    cuyos eigenvalores lambda_i - s están entre 0.1 y 2||H||_inf + 0.1:
    todos positivos, así que dos eigenvalores distintos tienen valores
    absolutos distintos, y ninguno es cero. Al final se le suma s.
    (Los eigenvalores repetidos de verdad no le estorban a QR.)

    H: Matrix de n x n. kmax: iteraciones de metodo_qr;
    con eigenvalores muy juntos hacen falta muchas.
    Regresa la lista de eigenvalores, de menor a mayor.
    """
    s = -(norma_infinito(H) + 0.1)
    return sorted(lam + s for lam in metodo_qr(menos_lambda_identidad(H, s), kmax))


###############################################
# 2. Comprobaciones
###############################################


def multiplicidades(valores, tolerancia: float = 1e-6):
    """Agrupa los valores (ordenados de menor a mayor) que difieren en
    menos de `tolerancia`. Regresa una lista de tuplas (valor, veces),
    con el valor promedio de cada grupo.
    """
    grupos = []
    for v in valores:
        if grupos and abs(v - grupos[-1][-1]) < tolerancia:
            grupos[-1].append(v)
        else:
            grupos.append([v])
    return [(sum(grupo) / len(grupo), len(grupo)) for grupo in grupos]


###############################################
# 3. Sin QR: las raíces de det(H - lambda I)
###############################################


def polinomio_caracteristico(H, lam: float):
    """p(lambda) = det(H - lambda I), o None si la LU sin pivoteo
    encuentra un pivote cero en ese lambda."""
    try:
        return determinante(menos_lambda_identidad(H, lam))
    except ZeroDivisionError:
        return None


def barrido_de_signo(H, a: float, b: float, puntos: int):
    """Busca intervalos [lambda_k, lambda_k+1] de una malla uniforme en
    [a, b] donde det(H - lambda I) cambia de signo: en cada uno hay (al
    menos) un eigenvalor. Regresa la lista de intervalos (tuplas).
    """
    paso = (b - a) / puntos
    intervalos = []
    lam_anterior = a
    p_anterior = polinomio_caracteristico(H, a)
    for k in range(1, puntos + 1):
        lam = a + k * paso
        p = polinomio_caracteristico(H, lam)
        if p is not None and p_anterior is not None and p * p_anterior < 0:
            intervalos.append((lam_anterior, lam))
        if p is not None:
            # Si p es None (pivote cero) seguimos comparando contra el
            # último valor que sí se pudo calcular.
            lam_anterior, p_anterior = lam, p
    return intervalos


def refinar(H, a: float, b: float, rondas: int = 10):
    """Afina un eigenvalor en [a, b] (donde det(H - lambda I) cambia de
    signo): divide el intervalo en 10, se queda con el pedazo donde
    sigue cambiando el signo, y repite. Cada ronda gana un dígito.
    Regresa el punto medio del último intervalo.
    """
    for _ in range(rondas):
        sub = barrido_de_signo(H, a, b, 10)
        if not sub:
            break  # no se encontró el cambio de signo (pivote cero)
        a, b = sub[0]
    return (a + b) / 2


def eigenvalores_por_barrido(H, puntos: int = 2000, rondas: int = 6):
    """Los eigenvalores de H (simétrica) que el barrido encuentra.

    Como H es simétrica, sus eigenvalores son reales y cumplen
    |lambda| <= ||H||_inf, así que basta barrer ese intervalo. (Se
    agranda un poquito para no caer justo en un eigenvalor en el borde.)
    H: Matrix. Regresa una lista ordenada de menor a mayor.
    """
    cota = 1.001 * norma_infinito(H) + 1e-12
    intervalos = barrido_de_signo(H, -cota, cota, puntos)
    return [refinar(H, a, b, rondas) for a, b in intervalos]


###############################################
# 4. Dos espines: el bloque 2x2
###############################################


def eigenvalores_bloque_2x2(a: float, b: float, d: float):
    """Eigenvalores de la matriz simétrica [[a, b], [b, d]], con el
    polinomio característico lambda^2 - (a + d) lambda + (a d - b^2):

        lambda = (a + d)/2 +- sqrt(((a - d)/2)^2 + b^2)
    """
    centro = (a + d) / 2
    radio = raiz_cuadrada(((a - d) / 2) ** 2 + b**2)
    return centro - radio, centro + radio


def eigenvalores_dos_espines(H):
    """Los 4 eigenvalores exactos del hamiltoniano de dos espines.

    En la base (arriba-arriba, arriba-abajo, abajo-arriba, abajo-abajo),
    H[0][0] y H[3][3] están solos en su renglón y columna (son
    eigenvalores), y el resto es un bloque 2x2. H: Matrix.
    Regresa los 4 eigenvalores ordenados.
    """
    menor, mayor = eigenvalores_bloque_2x2(H[1][1], H[1][2], H[2][2])
    return sorted([H[0][0], H[3][3], menor, mayor])
