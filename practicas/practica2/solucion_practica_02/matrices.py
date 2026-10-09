#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
matrices.py
-----------

Este módulo define la clase Matrix para representar matrices y
realizar operaciones básicas.

Incluye:
- Creación de matrices a partir de listas de listas.
- Suma, resta y multiplicación (escalar y matricial).
- Transpuesta.
- Métodos auxiliares para acceder a filas y columnas.
"""


class Matrix:
    """
    Clase para representar matrices y sus operaciones básicas.

    Attributes
    ----------
    data : list of list of numbers
        Almacena los elementos de la matriz, un renglón por sublista.
    rows : int
        Número de filas.
    cols : int
        Número de columnas.

    Methods
    -------
    shape():
        Devuelve una tupla con la dimensión de la matriz (filas, columnas).
    get_row(i):
        Devuelve la fila i (0-indexada).
    get_col(j):
        Devuelve la columna j (0-indexada).
    transpose():
        Devuelve la matriz transpuesta.
    copy():
        Devuelve una copia independiente de la matriz.
    __add__(other):
        Suma dos matrices del mismo tamaño.
    __sub__(other):
        Resta dos matrices del mismo tamaño.
    __mul__(other):
        Multiplica la matriz por un escalar o por otra matriz.
    """

    def __init__(self, data):
        """
        Inicializa una matriz a partir de una lista de listas.

        Parameters
        ----------
        data : list of list of numbers
            Lista de listas rectangular con los valores de la matriz.

        Raises
        ------
        ValueError
            Si `data` está vacía, o si las filas no tienen todas la
            misma longitud.
        """
        if not data:
            raise ValueError("una Matrix necesita al menos un renglón.")
        if not all(len(data[0]) == len(row) for row in data):
            raise ValueError("Todas las filas deben tener la misma longitud.")

        self.data = [list(row) for row in data]
        self.rows = len(self.data)
        self.cols = len(self.data[0])

    def __str__(self):
        """Representación bonita para impresión."""
        return "\n".join(["\t".join(map(str, row)) for row in self.data])

    def shape(self):
        """Devuelve la dimensión de la matriz como (filas, columnas)."""
        return (self.rows, self.cols)

    def get_row(self, i):
        """Devuelve la fila i (0-indexada; también acepta índices
        negativos, como las listas de Python)."""
        if not -self.rows <= i < self.rows:
            raise ValueError("Índice de fila fuera de rango.")
        return self.data[i]

    def get_col(self, j):
        """Devuelve la columna j (0-indexada; también acepta índices
        negativos, como las listas de Python)."""
        if not -self.cols <= j < self.cols:
            raise ValueError("Índice de columna fuera de rango.")
        return [row[j] for row in self.data]

    def copy(self):
        """
        Devuelve una copia profunda (deep copy) de la matriz.
        Útil para evitar modificar el objeto original en operaciones
        destructivas.
        """
        return Matrix([row[:] for row in self.data])

    def transpose(self):
        """Devuelve la transpuesta de la matriz."""
        result = [[self.data[j][i] for j in range(self.rows)] for i in range(self.cols)]
        return Matrix(result)

    def __add__(self, other):
        """
        Suma de matrices del mismo tamaño.

        Parameters
        ----------
        other : Matrix

        Returns
        -------
        Matrix

        Raises
        ------
        ValueError
            Si `self` y `other` no tienen las mismas dimensiones.
        """
        if self.shape() != other.shape():
            raise ValueError("Las matrices deben tener las mismas dimensiones.")
        return Matrix(
            [
                [self.data[i][j] + other.data[i][j] for j in range(self.cols)]
                for i in range(self.rows)
            ]
        )

    def __sub__(self, other):
        """
        Resta de matrices del mismo tamaño.

        Parameters
        ----------
        other : Matrix

        Returns
        -------
        Matrix

        Raises
        ------
        ValueError
            Si `self` y `other` no tienen las mismas dimensiones.
        """
        if self.shape() != other.shape():
            raise ValueError("Las matrices deben tener las mismas dimensiones.")
        return self + other * -1

    def __mul__(self, other):
        """
        Multiplicación por escalar o por otra matriz.

        Parameters
        ----------
        other : int, float, Matrix
            Escalar o matriz a multiplicar.

        Returns
        -------
        Matrix
            Nueva matriz con el resultado.

        Raises
        ------
        ValueError
            Si las dimensiones no son compatibles en multiplicación
            matricial (columnas de self != filas de other).
        TypeError
            Si el tipo de `other` no es soportado.
        """
        if isinstance(other, (int, float)):
            return Matrix(
                [
                    [self.data[i][j] * other for j in range(self.cols)]
                    for i in range(self.rows)
                ]
            )

        if isinstance(other, Matrix):
            if self.cols != other.rows:
                raise ValueError("Dimensiones incompatibles para multiplicación.")
            result = [
                [
                    sum(self.data[i][k] * other.data[k][j] for k in range(self.cols))
                    for j in range(other.cols)
                ]
                for i in range(self.rows)
            ]
            return Matrix(result)

        raise TypeError("Operación no soportada.")

    __rmul__ = __mul__


if __name__ == "__main__":
    A = Matrix([[1, 2, 3], [4, 5, 6]])
    print(f"Las componentes de la matriz A son {A.data}")
    print(f"La matriz A tiene {A.rows} renglones")
    print(f"La matriz A tiene {A.cols} columnas")

    B = Matrix([[7, 8, 9], [10, 11, 12]])

    print("Matriz A:")
    print(A)
    print(f"La forma de la matriz A es {A.shape()}")
    print(f"El primer renglón de la matriz A es {A.get_row(0)}")
    print(f"El segundo renglón de la matriz A es {A.get_row(1)}")
    print(f"La segunda columna de la matriz A es {A.get_col(1)}")

    print("\nMatriz B:")
    print(B)

    print("\nA + B:")
    print(A + B)

    print("\nA - B:")
    print(A - B)

    print("\nA * 2 (escalar):")
    print(A * 2)
    print("\n2 * A (escalar, del otro lado):")
    print(2 * A)

    C = Matrix([[1, 2], [3, 4], [5, 6]])

    print("\nC:")
    print(C)
    print("\nA * C (matricial):")
    print(A * C)

    print("\nTranspuesta de A:")
    print(A.transpose())
