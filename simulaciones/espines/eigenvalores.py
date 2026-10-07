#!/usr/bin/env python3
"""Eigenvalores de una matriz simétrica sin el método QR (ver README.md).

Todavía no hemos visto un método general para eigenvalores (QR). Aquí
usamos lo que sí tenemos, la descomposición LU de la unidad 09:

1. det(H - lambda I) con LU: el producto de la diagonal de U.
2. Los eigenvalores son las raíces de det(H - lambda I): se barre
   lambda en una malla y se buscan cambios de signo.
3. Cada raíz se afina subdividiendo el intervalo donde cambió el signo.
4. (Adelanto del tema de eigenvalores) Iteración inversa: con un
   valor aproximado sigma, resolver (H - sigma I) x_nuevo = x una y
   otra vez con LU converge al eigenvector, y de ahí sale el
   eigenvalor con precisión de máquina.

Limitaciones (ver README.md, "Lo que todavía no se puede"): un
eigenvalor con multiplicidad par no cambia el signo del determinante,
y dos eigenvalores en la misma celda de la malla se cancelan. Para
eso hace falta QR.
"""

from fiscomp.algebra_lineal import (
    determinante,
    descomposicion_lu,
    resolver_lu,
    mat_vec,
    norma_infinito,
    norma_infinito_vec,
)
from fiscomp.funciones_especiales import raiz_cuadrada

###############################################
# El polinomio característico, con el determinante de LU
###############################################


def menos_lambda_identidad(H, lam: float):
    """H - lambda I, como lista de listas (H también es lista de listas)."""
    n = len(H)
    return [[H[i][j] - (lam if i == j else 0.0) for j in range(n)] for i in range(n)]


def polinomio_caracteristico(H, lam: float):
    """p(lambda) = det(H - lambda I), o None si la LU sin pivoteo
    encuentra un pivote cero en ese lambda."""
    try:
        return determinante(menos_lambda_identidad(H, lam))
    except ZeroDivisionError:
        return None


###############################################
# Barrido y refinamiento
###############################################


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
    """Todos los eigenvalores de H (simétrica) que el barrido encuentra.

    Como H es simétrica, sus eigenvalores son reales y cumplen
    |lambda| <= ||H||_inf, así que basta barrer ese intervalo. (Se
    agranda un poquito para no caer justo en un eigenvalor en el borde.)
    H: lista de listas. Regresa una lista ordenada de menor a mayor.
    """
    cota = 1.001 * norma_infinito(H) + 1e-12
    intervalos = barrido_de_signo(H, -cota, cota, puntos)
    return [refinar(H, a, b, rondas) for a, b in intervalos]


###############################################
# Iteración inversa
###############################################


def iteracion_inversa(H, sigma: float, iteraciones: int = 30):
    """Eigenvalor de H más cercano a sigma, y su eigenvector.

    Se factoriza H - sigma I = L U una sola vez y se repite

        (H - sigma I) x_nuevo = x     (sustitución adelante y atrás)
        x <- x_nuevo / (su componente más grande)

    El vector converge al eigenvector cuyo eigenvalor está más cerca de
    sigma (es la componente que (H - sigma I)^-1 amplifica más). El
    eigenvalor se calcula al final con el cociente de Rayleigh,
    x.Hx / x.x, que para H simétrica es muy preciso.

    Regresa (eigenvalor, eigenvector).
    """
    n = len(H)
    L, U = descomposicion_lu(menos_lambda_identidad(H, sigma))
    # Vector inicial con todas sus componentes distintas, para no
    # empezar justo sin componente en el eigenvector buscado.
    x = [1.0 + 0.1 * i for i in range(n)]
    for _ in range(iteraciones):
        x_nuevo = resolver_lu(L, U, x)
        mas_grande = max(x_nuevo, key=abs)
        x = [x_i / mas_grande for x_i in x_nuevo]
    Hx = mat_vec(H, x)
    eigenvalor = sum(a * b for a, b in zip(x, Hx)) / sum(a * a for a in x)
    return eigenvalor, x


def residuo_eigen(H, eigenvalor: float, v):
    """||H v - lambda v||_inf / ||v||_inf: cero si (lambda, v) es un
    eigenpar exacto."""
    Hv = mat_vec(H, v)
    r = [Hv_i - eigenvalor * v_i for Hv_i, v_i in zip(Hv, v)]
    return norma_infinito_vec(r) / norma_infinito_vec(v)


###############################################
# Dos espines: el bloque 2x2
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
    """Los 4 eigenvalores del hamiltoniano de dos espines, sin barrido.

    En la base (arriba-arriba, arriba-abajo, abajo-arriba, abajo-abajo),
    H[0][0] y H[3][3] están solos en su renglón y columna (son
    eigenvalores), y el resto es un bloque 2x2. H: lista de listas.
    Regresa los 4 eigenvalores ordenados.
    """
    menor, mayor = eigenvalores_bloque_2x2(H[1][1], H[1][2], H[2][2])
    return sorted([H[0][0], H[3][3], menor, mayor])
