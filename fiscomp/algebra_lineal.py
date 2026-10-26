"""Álgebra lineal numérica, escrita desde cero con listas de listas.

Ver unidades/09_algebra_lineal/notas.md para la teoría de sistemas
lineales (normas, número de condición, sustituciones, eliminación
gaussiana, LU, Jacobi) y unidades/10_valores_y_vectores_propios/notas.md
para la de eigenvalores (potencia, potencia inversa, método QR) y la
SVD.

Convenciones:

- Una matriz es una lista de listas rectangular, un renglón por
  sublista: A[i][j] es el elemento del renglón i, columna j (índices
  desde 0, como en el libro y en Python). Si tienen una Matrix
  (fiscomp/matrices.py), pasen su atributo .data.
- Un vector es una lista simple de números: x[i].

Aquí vive lo que se reutiliza en el resto del tema: las sustituciones
para matrices triangulares, la eliminación gaussiana (con y sin
pivoteo), la descomposición LU con la inversa y el determinante, el
método de Jacobi, los métodos de la potencia (directo e inverso con
desplazamiento) y el método QR para eigenvalores, la descomposición
en valores singulares y la matriz de prueba del libro.
"""

from fiscomp.funciones_especiales import raiz_cuadrada

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


###############################################
# Eliminación gaussiana y descomposición LU
###############################################


def eliminacion_gaussiana(A, b):
    """Resuelve A x = b con eliminación gaussiana (gauelim en el libro).

    Fase de eliminación: para cada renglón pivote j = 0, ..., n-2 y
    cada renglón i = j+1, ..., n-1 de abajo,

        coeficiente = A_ij / A_jj
        renglón i  <-  renglón i - coeficiente * renglón j

    (en A y en b), hasta que A queda triangular superior. Luego,
    sustitución hacia atrás.

    A: matriz n x n (lista de listas). b: vector de n componentes.
    No modifica A ni b: trabaja sobre copias. Sin pivoteo: si algún
    pivote A_jj vale cero, levanta ZeroDivisionError (ver notas.md).
    Cuesta 2n^3/3 + 3n^2/2 - 13n/6 operaciones la eliminación, más n^2
    la sustitución hacia atrás.
    """
    A = [renglon[:] for renglon in A]
    b = b[:]
    n = len(b)
    for j in range(n - 1):
        for i in range(j + 1, n):
            coeficiente = A[i][j] / A[j][j]
            # En el libro, con NumPy: A[i,j:] -= coeficiente*A[j,j:]
            for k in range(j, n):
                A[i][k] -= coeficiente * A[j][k]
            b[i] -= coeficiente * b[j]
    return sustitucion_atras(A, b)


def descomposicion_lu(A):
    """Descomposición LU de Doolittle: A = L U (ludec en el libro).

    Es la misma eliminación gaussiana, pero sin b: lo que queda de A
    al final es U (triangular superior), y los coeficientes que se
    usaron para eliminar se guardan en L (triangular inferior, con
    unos en la diagonal): L_ij = coeficiente que anuló a A_ij.

    A: matriz n x n (lista de listas); no se modifica.
    Regresa (L, U), dos listas de listas. Sin pivoteo, como
    eliminacion_gaussiana. Cuesta 2n^3/3 + n^2/2 - 7n/6 operaciones.
    """
    n = len(A)
    U = [renglon[:] for renglon in A]
    L = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    for j in range(n - 1):
        for i in range(j + 1, n):
            coeficiente = U[i][j] / U[j][j]
            for k in range(j, n):
                U[i][k] -= coeficiente * U[j][k]
            L[i][j] = coeficiente
    return L, U


def resolver_lu(L, U, b):
    """Resuelve A x = b a partir de la descomposición A = L U:

        L y = b   (sustitución hacia adelante)
        U x = y   (sustitución hacia atrás)

    L, U: lo que regresa descomposicion_lu(A). b: vector de n
    componentes. Cuesta 2n^2 operaciones: con la misma L y U se pueden
    resolver muchos sistemas con distinto b sin repetir la
    descomposición, que es la parte cara.
    """
    y = sustitucion_adelante(L, b)
    return sustitucion_atras(U, y)



def inversa(A):
    """Matriz inversa con LU: la columna k de A^-1 es la solución de

        A x_k = e_k,

    con e_k la columna k de la identidad (porque A A^-1 = I). Son n
    sistemas con la misma A: una sola descomposición LU y n pares de
    sustituciones. Cuesta ~2n^3/3 + 2n^3 = 8n^3/3 operaciones.

    A: matriz n x n (lista de listas); no se modifica.
    Regresa A^-1 como lista de listas. Sin pivoteo (ver
    descomposicion_lu).
    """
    n = len(A)
    L, U = descomposicion_lu(A)
    columnas = []
    for k in range(n):
        e_k = [1.0 if i == k else 0.0 for i in range(n)]
        columnas.append(resolver_lu(L, U, e_k))
    # columnas[k][i] es el elemento (i, k) de la inversa: transponemos.
    return [[columnas[k][i] for k in range(n)] for i in range(n)]


def determinante(A):
    """Determinante con LU: det(A) = det(L) det(U), y como L tiene unos
    en la diagonal,

        det(A) = U_00 U_11 ... U_(n-1)(n-1).

    Cuesta ~2n^3/3 operaciones (contra ~n! con cofactores). Sin
    pivoteo: si aparece un pivote cero levanta ZeroDivisionError,
    aunque A no sea singular.
    """
    _, U = descomposicion_lu(A)
    producto = 1.0
    for i in range(len(U)):
        producto *= U[i][i]
    return producto


###############################################
# Pivoteo parcial
###############################################


def eliminacion_gaussiana_pivoteo(A, b):
    """Eliminación gaussiana con pivoteo parcial (gauelim_pivot en el
    libro).

    Igual que eliminacion_gaussiana, pero antes de usar el renglón j
    como pivote se busca, en la columna j y del renglón j para abajo,
    el elemento de mayor valor absoluto (renglón k), y se intercambian
    los renglones j y k de A y de b. Así el pivote nunca es cero (si A
    no es singular) y los coeficientes nunca pasan de 1 en valor
    absoluto.

    A: matriz n x n (lista de listas). b: vector de n componentes.
    No modifica A ni b. Regresa la solución x.
    """
    A = [renglon[:] for renglon in A]
    b = b[:]
    n = len(b)
    for j in range(n - 1):
        # k: el renglón (de j para abajo) con el |A[k][j]| más grande.
        # Con ">" estricto, si hay empate se queda el primero.
        k = j
        for m in range(j + 1, n):
            if abs(A[m][j]) > abs(A[k][j]):
                k = m
        if k != j:
            A[j], A[k] = A[k], A[j]
            b[j], b[k] = b[k], b[j]
        for i in range(j + 1, n):
            coeficiente = A[i][j] / A[j][j]
            for k_col in range(j, n):
                A[i][k_col] -= coeficiente * A[j][k_col]
            b[i] -= coeficiente * b[j]
    return sustitucion_atras(A, b)


###############################################
# Método iterativo de Jacobi
###############################################


def paso_jacobi(A, b, x):
    """Una iteración del método de Jacobi:

        x_nuevo_i = (b_i - sum_{j != i} A_ij x_j) / A_ii

    Todas las componentes nuevas se calculan con el x *anterior*.
    Regresa una lista nueva (no modifica x).
    """
    n = len(b)
    x_nuevo = []
    for i in range(n):
        suma = 0.0
        for j in range(n):
            if j != i:
                suma += A[i][j] * x[j]
        x_nuevo.append((b[i] - suma) / A[i][i])
    return x_nuevo


def cambio_relativo(x_viejo, x_nuevo):
    """Criterio de paro del libro (termcrit):

        sum_i |(x_nuevo_i - x_viejo_i) / x_nuevo_i|
    """
    return sum(abs((xn - xv) / xn) for xv, xn in zip(x_viejo, x_nuevo))


def jacobi(A, b, kmax=50, tol=1e-6):
    """Resuelve A x = b con el método iterativo de Jacobi, empezando en
    x = 0, hasta que cambio_relativo sea menor que tol.

    A: matriz n x n (lista de listas), con diagonal sin ceros. Converge
       seguro si A es diagonalmente dominante.
    kmax: número máximo de iteraciones (entero); tol: tolerancia.
    Regresa (x, k): la solución y el número de iteraciones que tomó,
    o (None, kmax) si no convergió.
    """
    x = [0.0] * len(b)
    for k in range(1, kmax):
        x_nuevo = paso_jacobi(A, b, x)
        error = cambio_relativo(x, x_nuevo)
        x = x_nuevo
        if error < tol:
            break
    else:
        # Solo llega aquí si el for terminó sin break: no convergió.
        return None, kmax
    return x, k


###############################################
# Eigenvalores: método de la potencia
###############################################


def producto_punto(x, y):
    """Producto punto de dos vectores: sum_i x_i y_i."""
    suma = 0.0
    for x_i, y_i in zip(x, y):
        suma += x_i * y_i
    return suma


def normalizar(x):
    """Copia de x dividida entre su norma euclídea (norma 1)."""
    norma = norma_euclidea(x)
    return [x_i / norma for x_i in x]


def cociente_rayleigh(A, x):
    """Cociente de Rayleigh (ec. 4.133 del libro):

        mu(x) = x.(A x) / x.x

    Si x es un eigenvector de A, da exactamente su eigenvalor; si x
    está cerca de un eigenvector, da una buena aproximación.
    """
    return producto_punto(x, mat_vec(A, x)) / producto_punto(x, x)


def potencia(A, kmax=6):
    """Método de la potencia (power en el libro): el eigenvalor de mayor
    valor absoluto de A (el dominante) y su eigenvector. Para k = 1, 2, ...

        z = A q
        q = z / ||z||

    empezando con q = (1, 1, ..., 1) normalizado. Al final, el
    eigenvalor es el cociente de Rayleigh de q.

    A: matriz n x n (lista de listas). kmax: se hacen kmax - 1
    iteraciones, sin criterio de paro (como en el libro).
    Regresa (eigenvalor, eigenvector), con el eigenvector de norma 1.
    """
    q = normalizar([1.0] * len(A))
    for k in range(1, kmax):
        z = mat_vec(A, q)
        q = normalizar(z)
    return cociente_rayleigh(A, q), q


def potencia_inversa(A, desplazamiento=20.0, kmax=200, tol=1e-8):
    """Método de la potencia inversa con desplazamiento (invpowershift en
    el libro): el eigenvalor de A más cercano a `desplazamiento` (s) y su
    eigenvector. Se factoriza A - s I = L U una sola vez, y para
    k = 1, 2, ...

        (A - s I) z = q     (con resolver_lu: sustitución adelante y atrás)
        q = z / ||z||

    hasta que cambio_relativo entre dos q seguidos sea menor que tol.
    Si q cambió de signo respecto al anterior se le vuelve a cambiar,
    porque con un eigenvalor de A - s I negativo el signo de q alterna
    en cada iteración y el criterio de paro nunca se cumpliría.

    A: matriz n x n (lista de listas). desplazamiento: número real s;
    no debe ser exactamente un eigenvalor (A - s I sería singular).
    kmax: número máximo de iteraciones (entero); tol: tolerancia.
    Regresa (eigenvalor, eigenvector, k): el eigenvalor (cociente de
    Rayleigh con la A original), el eigenvector de norma 1 y el número
    de iteraciones, o (None, None, kmax) si no convergió.
    """
    n = len(A)
    A_desplazada = [[A[i][j] - (desplazamiento if i == j else 0.0) for j in range(n)] for i in range(n)]
    L, U = descomposicion_lu(A_desplazada)
    q = normalizar([1.0] * n)
    for k in range(1, kmax):
        z = resolver_lu(L, U, q)
        q_nuevo = normalizar(z)
        if producto_punto(q, q_nuevo) < 0:
            q_nuevo = [-q_i for q_i in q_nuevo]
        error = cambio_relativo(q, q_nuevo)
        q = q_nuevo
        if error < tol:
            break
    else:
        # Solo llega aquí si el for terminó sin break: no convergió.
        return None, None, kmax
    return cociente_rayleigh(A, q), q, k


###############################################
# Eigenvalores: método QR
###############################################


def transpuesta(A):
    """Transpuesta de una matriz: (A^T)[i][j] = A[j][i]. También sirve
    para pasar de "lista de renglones" a "lista de columnas".
    """
    return [list(columna) for columna in zip(*A)]


def mat_mat(A, B):
    """Producto de matrices C = A B, con C[i][j] = sum_k A[i][k] B[k][j].
    Para matrices n x n cuesta ~2n^3 operaciones.
    """
    B_T = transpuesta(B)
    return [[producto_punto(renglon, columna) for columna in B_T] for renglon in A]


def descomposicion_qr(A):
    """Descomposición QR (qrdec en el libro) con Gram-Schmidt clásico:
    A = Q R, con Q ortogonal y R triangular superior. Para cada columna
    a_j de A (j = 0, 1, ..., n-1):

        R_ij = q_i . a_j                (i = 0, ..., j-1)
        a'_j = a_j - sum_i R_ij q_i
        R_jj = ||a'_j||
        q_j  = a'_j / R_jj

    A: matriz n x n (lista de listas) con columnas linealmente
    independientes. Regresa (Q, R), las dos como listas de listas.
    """
    n = len(A)
    columnas = transpuesta(A)
    q = []  # columnas de Q, ya ortonormales
    R = [[0.0] * n for _ in range(n)]
    for j in range(n):
        a_prima = columnas[j][:]
        for i in range(j):
            R[i][j] = producto_punto(q[i], columnas[j])
            a_prima = [a - R[i][j] * q_i for a, q_i in zip(a_prima, q[i])]
        R[j][j] = norma_euclidea(a_prima)
        q.append([a / R[j][j] for a in a_prima])
    return transpuesta(q), R


def metodo_qr(A, kmax=100):
    """Método QR (qrmet en el libro): todos los eigenvalores de A. Empezando
    con A^(0) = A, para k = 1, 2, ...

        A^(k-1) = Q R     (descomposicion_qr)
        A^(k)   = R Q     (multiplicar en orden inverso)

    A^(k) tiene los mismos eigenvalores que A y tiende a una matriz
    triangular superior, así que los eigenvalores se leen de su diagonal.

    A: matriz n x n (lista de listas) con eigenvalores reales y
    distintos. kmax: se hacen kmax - 1 iteraciones, sin criterio de
    paro (como en el libro). Regresa la lista de la diagonal de A^(k).
    """
    A_k = [renglon[:] for renglon in A]
    for k in range(1, kmax):
        Q, R = descomposicion_qr(A_k)
        A_k = mat_mat(R, Q)
    return [A_k[i][i] for i in range(len(A))]


def eigenpares(A, eps=1e-12):
    """Todos los eigenvalores y eigenvectores de A (eig en el libro).
    Los eigenvalores de metodo_qr se usan como desplazamientos de
    potencia_inversa, que da el eigenvector (y refina el eigenvalor).
    Al desplazamiento se le suma eps para que A - s I no sea singular.

    A: matriz n x n (lista de listas) con eigenvalores reales y distintos.
    Regresa (eigenvalores, eigenvectores): una lista de n números y una
    lista de n vectores de norma 1; eigenvectores[i] es el de
    eigenvalores[i]. (El libro los regresa como columnas de una matriz.)
    """
    eigenvalores = []
    eigenvectores = []
    for lam_qr in metodo_qr(A):
        lam, v, _ = potencia_inversa(A, desplazamiento=lam_qr + eps)
        eigenvalores.append(lam)
        eigenvectores.append(v)
    return eigenvalores, eigenvectores


###############################################
# Descomposición en valores singulares (SVD)
###############################################


def matriz_h(A):
    """Matriz simétrica H de 2n x 2n (makeH en el libro):

        H = | 0  A^T |
            | A   0  |

    Sus eigenvalores son +sigma_i y -sigma_i (los valores singulares de
    A con los dos signos), y el eigenvector de +sigma_i es (v_i, u_i)
    / sqrt(2): arriba el vector singular derecho, abajo el izquierdo.

    A: matriz n x n (lista de listas). Regresa H como lista de listas.
    """
    n = len(A)
    H = [[0.0] * (2 * n) for _ in range(2 * n)]
    for i in range(n):
        for j in range(n):
            H[n + i][j] = A[i][j]  # bloque de abajo a la izquierda: A
            H[j][n + i] = A[i][j]  # bloque de arriba a la derecha: A^T
    return H


def svd(A, desplazamiento=1.0, eps=1e-12):
    """Descomposición en valores singulares (svd en el libro): A = U S V^T,
    con U y V ortogonales y S diagonal, sin pasar por A^T A.

    1. Arma H = matriz_h(A), cuyos eigenvalores son +-sigma_i.
    2. Los eigenvalores de H salen de metodo_qr aplicado a H - s I: el
       desplazamiento s separa los pares +sigma y -sigma, que tienen el
       mismo valor absoluto y con los que QR no converge. Luego se le
       vuelve a sumar s.
    3. Se quedan los eigenvalores >= 0 (los sigma_i, de mayor a menor),
       y para cada uno potencia_inversa da el eigenvector w de H (y
       refina sigma). Los primeros n elementos de w, por sqrt(2), son
       v_i; los últimos n, por sqrt(2), son u_i.

    A: matriz n x n (lista de listas) con valores singulares distintos
    y no muy separados en órdenes de magnitud (limitación de metodo_qr).
    desplazamiento: número real s, el corrimiento del paso 2.
    eps: como en eigenpares, para que H - sigma I no sea singular.
    Regresa (U, sigmas, V_T): U como lista de listas, la lista de los
    valores singulares y V^T como lista de listas. (El libro regresa
    los valores singulares como una matriz diagonal.)
    """
    n = len(A)
    H = matriz_h(A)
    H_desplazada = [[H[i][j] - (desplazamiento if i == j else 0.0) for j in range(2 * n)] for i in range(2 * n)]
    eigenvalores_h = [lam + desplazamiento for lam in metodo_qr(H_desplazada)]
    sigmas_qr = sorted([lam for lam in eigenvalores_h if lam >= 0], reverse=True)
    sigmas = []
    columnas_u = []
    columnas_v = []
    for sigma_qr in sigmas_qr:
        sigma, w, _ = potencia_inversa(H, desplazamiento=sigma_qr + eps)
        sigmas.append(sigma)
        columnas_v.append([raiz_cuadrada(2.0) * w_i for w_i in w[:n]])
        columnas_u.append([raiz_cuadrada(2.0) * w_i for w_i in w[n:]])
    # Las columnas de V son los renglones de V^T: columnas_v ya es V^T.
    return transpuesta(columnas_u), sigmas, columnas_v


###############################################
# Matriz de prueba
###############################################


def crear_prueba(n, val):
    """Matriz de prueba del libro (testcreate): A_ij = sqrt(val + n i + j),
    y b_j = (A_0j)^2.1. No es simétrica. Regresa (A, b).
    """
    A = [[raiz_cuadrada(val + n * i + j) for j in range(n)] for i in range(n)]
    b = [a_0j**2.1 for a_0j in A[0]]
    return A, b
