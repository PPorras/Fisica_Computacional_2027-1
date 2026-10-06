#!/usr/bin/env python3
"""Pruebas para la Práctica 2 (ver practica_02.md), con `unittest`.

Se corre con `python3 practicas/practica2/pruebas_practica_02.py`, con el
entorno virtual activado (`source .venv/bin/activate` desde la raíz
del repositorio) y `fiscomp` instalado en modo editable
(`pip install -e .`; ver README.md).

`__add__`, `__sub__` y `__mul__` empiezan con un `raise
NotImplementedError` (el TODO de la práctica), así que hasta que los
completen, las pruebas que dependen de ellos van a fallar con `ERROR`
en vez de `OK` -- es la señal de que a esa parte le falta código, no
de que esté mal. Cada valor de prueba se revisa con `subTest`, así que
si uno falla los demás se siguen corriendo, y los mensajes no dicen el
valor esperado ni el obtenido -- esa parte les toca investigarla a
ustedes.
"""

import unittest

from fiscomp.matrices import Matrix


class TestConstruccion(unittest.TestCase):
    def test_shape_y_datos(self):
        m = Matrix([[1, 2, 3], [4, 5, 6]])
        with self.subTest(caso="shape"):
            self.assertEqual(m.shape(), (2, 3))
        with self.subTest(caso="get_row"):
            self.assertEqual(m.get_row(0), [1, 2, 3])
            self.assertEqual(m.get_row(1), [4, 5, 6])
        with self.subTest(caso="get_col"):
            self.assertEqual(m.get_col(0), [1, 4])
            self.assertEqual(m.get_col(2), [3, 6])

    def test_sin_renglones(self):
        with self.assertRaises(ValueError, msg="Matrix sin renglones debería fallar"):
            Matrix([])

    def test_renglones_de_distinta_longitud(self):
        with self.assertRaises(
            ValueError, msg="Matrix con renglones de distinta longitud debería fallar"
        ):
            Matrix([[1, 2, 3], [4, 5]])

    def test_str_muestra_las_entradas(self):
        m = Matrix([[1, 2], [3, 4]])
        texto = str(m)
        for valor in (1, 2, 3, 4):
            with self.subTest(valor=valor):
                self.assertIn(
                    str(valor), texto, msg=f"print(matriz) no muestra el valor {valor}"
                )


class TestTransposeYCopy(unittest.TestCase):
    def test_transpose(self):
        m = Matrix([[1, 2, 3], [4, 5, 6]])
        t = m.transpose()
        with self.subTest(caso="shape"):
            self.assertEqual(t.shape(), (3, 2))
        with self.subTest(caso="datos"):
            self.assertEqual(t.data, [[1, 4], [2, 5], [3, 6]])

    def test_copy_es_independiente(self):
        m = Matrix([[1, 2], [3, 4]])
        c = m.copy()
        c.data[0][0] = 999
        self.assertEqual(
            m.data[0][0], 1, msg="modificar la copia no debería afectar el original"
        )


class TestSumaResta(unittest.TestCase):
    def setUp(self):
        self.a = Matrix([[1, 2], [3, 4]])
        self.b = Matrix([[5, 6], [7, 8]])

    def test_suma(self):
        self.assertEqual((self.a + self.b).data, [[6, 8], [10, 12]])

    def test_resta(self):
        self.assertEqual((self.a - self.b).data, [[-4, -4], [-4, -4]])

    def test_dimensiones_distintas(self):
        c = Matrix([[1, 2, 3]])
        with self.subTest(operacion="suma"):
            with self.assertRaises(
                ValueError, msg="sumar dimensiones distintas debería fallar"
            ):
                self.a + c
        with self.subTest(operacion="resta"):
            with self.assertRaises(
                ValueError, msg="restar dimensiones distintas debería fallar"
            ):
                self.a - c


class TestMultiplicacionEscalar(unittest.TestCase):
    def test_escalar_por_derecha_e_izquierda(self):
        m = Matrix([[1, 2], [3, 4]])
        esperado = [[2, 4], [6, 8]]
        with self.subTest(orden="matriz * escalar"):
            self.assertEqual((m * 2).data, esperado)
        with self.subTest(orden="escalar * matriz"):
            self.assertEqual((2 * m).data, esperado)

    def test_tipo_no_soportado(self):
        m = Matrix([[1, 2], [3, 4]])
        with self.assertRaises(
            TypeError, msg="multiplicar por algo que no es número ni Matrix debería fallar"
        ):
            m * "hola"


class TestMultiplicacionDeMatrices(unittest.TestCase):
    def test_matrices_cuadradas(self):
        a = Matrix([[1, 2], [3, 4]])
        b = Matrix([[5, 6], [7, 8]])
        self.assertEqual((a * b).data, [[19, 22], [43, 50]])

    def test_matrices_no_cuadradas(self):
        c = Matrix([[1, 2, 3], [4, 5, 6]])  # 2x3
        d = Matrix([[7, 8], [9, 10], [11, 12]])  # 3x2
        self.assertEqual((c * d).data, [[58, 64], [139, 154]])

    def test_dimensiones_internas_incompatibles(self):
        a = Matrix([[1, 2], [3, 4]])  # 2x2
        c = Matrix([[1, 2, 3], [4, 5, 6]])  # 2x3
        with self.assertRaises(
            ValueError,
            msg="multiplicar con dimensiones internas incompatibles debería fallar",
        ):
            c * a  # 2x3 * 2x2: columnas de c (3) != renglones de a (2)


if __name__ == "__main__":
    unittest.main(verbosity=2)
