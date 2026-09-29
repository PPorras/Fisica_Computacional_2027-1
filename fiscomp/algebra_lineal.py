"""Álgebra lineal numérica, escrita desde cero con listas de listas.

Ver unidades/09_matrices/notas.md para la teoría (normas, número de
condición, sustitución hacia adelante y hacia atrás).

Convenciones:

- Una matriz es una lista de listas rectangular, un renglón por
  sublista: A[i][j] es el elemento del renglón i, columna j (índices
  desde 0, como en el libro y en Python). Si tienen una Matrix
  (fiscomp/matrices.py), pasen su atributo .data.
- Un vector es una lista simple de números: x[i].

Aquí vive lo que se reutiliza en el resto del tema (eliminación
gaussiana, LU, ...), empezando por las sustituciones para matrices
triangulares.
"""

###############################################
# Producto matriz-vector y residuo
###############################################


def mat_vec(A, x):
    """Producto matriz-vector y = A x, con y[i] = sum_j A[i][j] x[j].

    A: matriz n x m (lista de listas).
    x: vector de m componentes (lista).
    Regresa una lista de n componentes. Cuesta 2n^2 - n operaciones
    de punto flotante cuando A es cuadrada (ver notas.md).
    """
    if len(A[0]) != len(x):
        raise ValueError(
            f"dimensiones incompatibles: A tiene {len(A[0])} columnas "
            f"y x tiene {len(x)} componentes"
        )
    y = []
    for renglon in A:
        suma = 0.0
        for a_ij, x_j in zip(renglon, x):
            suma += a_ij * x_j
        y.append(suma)
    return y


def residuo(A, x, b):
    """Vector residuo r = b - A x (ec. 4.17 del libro).

    Si x es la solución exacta de A x = b, r es el vector cero. Un
    residuo chico NO garantiza que x esté cerca de la solución
    exacta (ver el ejemplo de Kahan en notas.md).
    """
    Ax = mat_vec(A, x)
    return [b_i - Ax_i for b_i, Ax_i in zip(b, Ax)]


###############################################
# Normas de vectores y de matrices
###############################################


def norma_euclidea(x):
    """Norma euclídea de un vector: sqrt(sum_i |x_i|^2)."""
    return sum(abs(x_i) ** 2 for x_i in x) ** 0.5


def norma_infinito_vec(x):
    """Norma infinito (de máxima magnitud) de un vector: max_i |x_i|."""
    return max(abs(x_i) for x_i in x)


def norma_frobenius(A):
    """Norma de Frobenius de una matriz: sqrt(sum_i sum_j |A_ij|^2)."""
    return sum(abs(a_ij) ** 2 for renglon in A for a_ij in renglon) ** 0.5


def norma_infinito(A):
    """Norma infinito de una matriz (máxima suma de renglón):

        max_i sum_j |A_ij|
    """
    return max(sum(abs(a_ij) for a_ij in renglon) for renglon in A)


###############################################
# Matrices triangulares
###############################################


def triangular_inferior(A):
    """Copia de A con ceros arriba de la diagonal (como numpy.tril)."""
    n = len(A)
    return [[A[i][j] if j <= i else 0.0 for j in range(n)] for i in range(n)]


def triangular_superior(A):
    """Copia de A con ceros abajo de la diagonal (como numpy.triu)."""
    n = len(A)
    return [[A[i][j] if j >= i else 0.0 for j in range(n)] for i in range(n)]


def sustitucion_adelante(L, b):
    """Resuelve L x = b con L triangular inferior (forsub en el libro):

        x_i = (b_i - sum_{j=0}^{i-1} L_ij x_j) / L_ii,  i = 0, 1, ..., n-1

    L: matriz n x n triangular inferior (solo se leen los elementos
       de la diagonal y de abajo de ella).
    b: vector de n componentes.
    Regresa la solución x como lista. Cuesta n^2 operaciones.
    """
    n = len(b)
    x = [0.0] * n
    for i in range(n):
        suma = 0.0
        for j in range(i):  # para i = 0 este ciclo no hace nada
            suma += L[i][j] * x[j]
        x[i] = (b[i] - suma) / L[i][i]
    return x


def sustitucion_atras(U, b):
    """Resuelve U x = b con U triangular superior (backsub en el libro):

        x_i = (b_i - sum_{j=i+1}^{n-1} U_ij x_j) / U_ii,  i = n-1, ..., 1, 0

    U: matriz n x n triangular superior (solo se leen los elementos
       de la diagonal y de arriba de ella).
    b: vector de n componentes.
    Regresa la solución x como lista. Cuesta n^2 operaciones.
    """
    n = len(b)
    x = [0.0] * n
    for i in reversed(range(n)):
        suma = 0.0
        for j in range(i + 1, n):  # para i = n-1 este ciclo no hace nada
            suma += U[i][j] * x[j]
        x[i] = (b[i] - suma) / U[i][i]
    return x
