#!/usr/bin/env python3
"""Partículas de espín 1/2: matrices de Pauli, producto de Kronecker y
hamiltonianos de dos y tres espines (ver README.md de esta carpeta;
libro de Gezerlis, sección 4.6, códigos 4.12 a 4.14).

Todo se construye con la clase Matrix (fiscomp/matrices.py), sin
NumPy. Las matrices pueden tener números complejos (sigma_y los
tiene): Python los trae integrados, se escriben como 1j.

Unidades: hbar = 1.
"""

from fiscomp.matrices import Matrix

HBAR = 1.0

###############################################
# Matrices de Pauli y producto de Kronecker
###############################################


def matrices_de_pauli():
    """Regresa (sigma_x, sigma_y, sigma_z) como Matrix de 2x2:

        sigma_x = [[0, 1], [1, 0]]
        sigma_y = [[0, -i], [i, 0]]
        sigma_z = [[1, 0], [0, -1]]
    """
    sigma_x = Matrix([[0.0, 1.0], [1.0, 0.0]])
    sigma_y = Matrix([[0.0, -1j], [1j, 0.0]])
    sigma_z = Matrix([[1.0, 0.0], [0.0, -1.0]])
    return sigma_x, sigma_y, sigma_z


def identidad(n: int):
    """Matriz identidad n x n, como Matrix."""
    return Matrix([[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)])


def kron(U, V):
    """Producto de Kronecker W = U (x) V de dos matrices cuadradas:

        W[a][b] = U[i][k] * V[j][l],   con a = p i + j,  b = p k + l,

    donde U es n x n y V es p x p, así que W es (n p) x (n p).

    U, V: Matrix cuadradas. Regresa una Matrix.
    Cuesta n^2 p^2 multiplicaciones (cuatro ciclos anidados).
    """
    n = U.rows
    p = V.rows
    W = [[0.0] * (n * p) for _ in range(n * p)]
    for i in range(n):
        for k in range(n):
            for j in range(p):
                for l in range(p):
                    W[p * i + j][p * k + l] = U.data[i][k] * V.data[j][l]
    return Matrix(W)


def parte_real(M):
    """Copia de M con solo la parte real de cada elemento.

    El libro hace H.real con NumPy; Matrix no tiene ese atributo.
    complex(x).real funciona igual si x es float o complejo.
    """
    return Matrix([[complex(x).real for x in renglon] for renglon in M.data])


def es_simetrica(M, tolerancia: float = 1e-14):
    """True si |M[i][j] - M[j][i]| <= tolerancia para todo i, j."""
    return all(
        abs(M.data[i][j] - M.data[j][i]) <= tolerancia
        for i in range(M.rows)
        for j in range(M.cols)
    )


###############################################
# Operadores de espín
###############################################


def operadores_de_espin(posicion: int, numero_de_particulas: int):
    """Operadores (S_x, S_y, S_z) de una partícula dentro de un sistema
    de varias partículas: S = (hbar/2) sigma en la posición de esa
    partícula, e identidades 2x2 en las demás, unidas con kron.

    Por ejemplo, con 3 partículas, la segunda (posicion = 1) da
    S_IIy = I (x) S_y (x) I.

    posicion: 0 para la partícula I, 1 para la II, etc.
    numero_de_particulas: cuántas partículas hay en total.
    Regresa una lista de tres Matrix de 2^N x 2^N.
    """
    iden = identidad(2)
    operadores = []
    for sigma in matrices_de_pauli():
        factores = [sigma * (HBAR / 2) if k == posicion else iden for k in range(numero_de_particulas)]
        S = factores[0]
        for factor in factores[1:]:
            S = kron(S, factor)
        operadores.append(S)
    return operadores


def producto_punto_espin(S_A, S_B):
    """S_A . S_B = S_Ax S_Bx + S_Ay S_By + S_Az S_Bz.

    El libro usa sum([...]) con NumPy; con Matrix eso no funciona:
    sum() empieza sumando 0 + matriz, y Matrix no sabe sumarse a un
    número (no tiene __radd__). Por eso sumamos los tres términos a
    mano.
    """
    return S_A[0] * S_B[0] + S_A[1] * S_B[1] + S_A[2] * S_B[2]


###############################################
# Hamiltonianos
###############################################


def dos_espines(omega_I: float, omega_II: float, gamma: float):
    """Hamiltoniano de dos espines 1/2 (ec. 4.232 del libro), 4x4:

        H = -omega_I S_Iz - omega_II S_IIz + gamma S_I . S_II

    Cada espín interactúa con el campo magnético externo (omega_I,
    omega_II) y los dos interactúan entre sí (gamma).
    Regresa una Matrix real.
    """
    S_I = operadores_de_espin(0, 2)
    S_II = operadores_de_espin(1, 2)
    H = S_I[2] * (-omega_I) + S_II[2] * (-omega_II) + producto_punto_espin(S_I, S_II) * gamma
    # Los productos S_Iy S_IIy tienen i * i = -1: todo queda real.
    return parte_real(H)


def tres_espines(omega_I: float, omega_II: float, omega_III: float, gamma: float):
    """Hamiltoniano de tres espines 1/2 (ec. 4.234 del libro), 8x8:

        H = -omega_I S_Iz - omega_II S_IIz - omega_III S_IIIz
            + gamma (S_I . S_II + S_I . S_III + S_II . S_III)

    Regresa una Matrix real.
    """
    S_I = operadores_de_espin(0, 3)
    S_II = operadores_de_espin(1, 3)
    S_III = operadores_de_espin(2, 3)
    H = S_I[2] * (-omega_I) + S_II[2] * (-omega_II) + S_III[2] * (-omega_III)
    interaccion = (
        producto_punto_espin(S_I, S_II)
        + producto_punto_espin(S_I, S_III)
        + producto_punto_espin(S_II, S_III)
    )
    return parte_real(H + interaccion * gamma)


if __name__ == "__main__":
    # La prueba del libro (código 4.12): sigma_x (x) (matriz de unos 3x3).
    sigma_x, sigma_y, sigma_z = matrices_de_pauli()
    unos = Matrix([[1.0] * 3 for _ in range(3)])
    print("sigma_x (x) unos:")
    print(parte_real(kron(sigma_x, unos)))
    print("\nunos (x) sigma_x (el producto de Kronecker no conmuta):")
    print(parte_real(kron(unos, sigma_x)))
