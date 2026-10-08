#!/usr/bin/env python3
"""Pruebas para la Práctica 4 (ver practica_04.md), con `unittest`.

Se corre desde la raíz del repositorio con

    python3 practicas/practica4/pruebas_practica_04.py

con el entorno virtual activado (`source .venv/bin/activate`) y
`fiscomp` instalado en modo editable (`pip install -e .`; ver
README.md).

Mientras una clase no exista todavía, las pruebas que la usan aparecen
como `skipped` con un mensaje que dice qué falta -- no es un error, es
la señal de que esa parte no está hecha. Como en las prácticas
anteriores, los mensajes no dicen el valor esperado ni el obtenido:
esa parte les toca investigarla a ustedes.

Qué se revisa:

- Ejercicio 1: LowerTriangularMatrix y UpperTriangularMatrix (en
  fiscomp/matrices.py) heredan de Matrix, validan su forma, y
  diagonal(), determinant(), is_singular(), transpose().
- Ejercicio 2: SistemaTriangular y MatrizSingularError (en
  fiscomp/algebra_lineal.py): soluciones correctas, excepciones, y el
  conteo self.flops.
- Ejercicio 3: corre unidades/09_algebra_lineal/complejidad.py y revisa
  datos/conteo_flops.dat (los conteos deben coincidir exactamente con
  las fórmulas) y el formato de datos/tiempos.dat.
- Ejercicio 4: PivoteCeroError, EliminacionGaussiana, SistemaLineal y
  FactorizacionLU (en fiscomp/algebra_lineal.py): U, b transformado y
  L (con la bandera coeficientes), soluciones, L U = A, determinante,
  inversa y excepciones. No revisa el diseño (que cada clase use las
  anteriores): eso se revisa a mano.

No revisa los docstrings, ni las respuestas escritas del Ejercicio 3:
eso se revisa a mano.
"""

import subprocess
import sys
import unittest
from pathlib import Path

from fiscomp.matrices import Matrix

try:
    from fiscomp.matrices import LowerTriangularMatrix, UpperTriangularMatrix
except ImportError:
    LowerTriangularMatrix = None
    UpperTriangularMatrix = None

try:
    from fiscomp.algebra_lineal import SistemaTriangular, MatrizSingularError
except ImportError:
    SistemaTriangular = None
    MatrizSingularError = None

try:
    from fiscomp.algebra_lineal import PivoteCeroError, EliminacionGaussiana
except ImportError:
    PivoteCeroError = None
    EliminacionGaussiana = None

try:
    from fiscomp.algebra_lineal import SistemaLineal
except ImportError:
    SistemaLineal = None

try:
    from fiscomp.algebra_lineal import FactorizacionLU
except ImportError:
    FactorizacionLU = None

RAIZ = Path(__file__).resolve().parent.parent.parent
CARPETA_UNIDAD = RAIZ / "unidades" / "09_algebra_lineal"
RUTA_SCRIPT = CARPETA_UNIDAD / "complejidad.py"
RUTA_FLOPS = CARPETA_UNIDAD / "datos" / "conteo_flops.dat"
RUTA_TIEMPOS = CARPETA_UNIDAD / "datos" / "tiempos.dat"

# Matrices de prueba
L_DATOS = [[2.0, 0.0, 0.0], [1.0, 4.0, 0.0], [-1.0, 3.0, 5.0]]
U_DATOS = [[2.0, 1.0, -1.0], [0.0, 4.0, 3.0], [0.0, 0.0, 5.0]]
X_EXACTA = [1.0, -2.0, 3.0]

# Ejercicio 4: el ejemplo 3x3 de notas.md (solución [4, -2, 2], det 40,
# con su U, su b transformado y su L, todos exactos en punto flotante),
# una matriz donde la eliminación deja un "cero" de -1.1e-16, una con
# pivote cero sin ser singular, y una singular.
A_DATOS = [[2.0, 1.0, 1.0], [1.0, 1.0, -2.0], [5.0, 10.0, 5.0]]
B_DATOS = [8.0, -2.0, 10.0]
X_DATOS = [4.0, -2.0, 2.0]
U_ELIMINADA = [[2.0, 1.0, 1.0], [0.0, 0.5, -2.5], [0.0, 0.0, 40.0]]
B_ELIMINADO = [8.0, -6.0, 80.0]
L_COEFICIENTES = [[1.0, 0.0, 0.0], [0.5, 1.0, 0.0], [2.5, 15.0, 1.0]]
A_REDONDEO = [[0.3, 1.0], [0.7, 2.0]]
A_PIVOTE_CERO = [[0.0, -1.0], [1.0, 1.0]]
A_SINGULAR = [[1.0, 2.0], [2.0, 4.0]]


def producto(A, x):
    return [sum(a_ij * x_j for a_ij, x_j in zip(renglon, x)) for renglon in A]


def requiere_triangulares(prueba):
    if LowerTriangularMatrix is None:
        prueba.skipTest(
            "todavía no existen LowerTriangularMatrix / UpperTriangularMatrix "
            "en fiscomp/matrices.py (Ejercicio 1)"
        )


def requiere_ejercicio_4(prueba, clase, nombre):
    requiere_sistema(prueba)
    if PivoteCeroError is None or clase is None:
        prueba.skipTest(
            f"todavía no existen PivoteCeroError / {nombre} en "
            "fiscomp/algebra_lineal.py (Ejercicio 4)"
        )


def requiere_sistema(prueba):
    requiere_triangulares(prueba)
    if SistemaTriangular is None:
        prueba.skipTest(
            "todavía no existen SistemaTriangular / MatrizSingularError "
            "en fiscomp/algebra_lineal.py (Ejercicio 2)"
        )


###############################################
# Ejercicio 1
###############################################


class TestMatricesTriangulares(unittest.TestCase):
    def setUp(self):
        requiere_triangulares(self)

    def test_heredan_de_matrix(self):
        for clase in (LowerTriangularMatrix, UpperTriangularMatrix):
            with self.subTest(clase=clase.__name__):
                self.assertTrue(
                    issubclass(clase, Matrix), msg=f"{clase.__name__} debe heredar de Matrix"
                )

    def test_construccion_valida(self):
        L = LowerTriangularMatrix(L_DATOS)
        U = UpperTriangularMatrix(U_DATOS)
        with self.subTest(caso="shape inferior"):
            self.assertEqual(L.shape(), (3, 3))
        with self.subTest(caso="data superior"):
            self.assertEqual(U.data, U_DATOS)
        with self.subTest(caso="métodos heredados"):
            self.assertEqual(L.get_row(1), L_DATOS[1])

    def test_no_cuadrada(self):
        for clase in (LowerTriangularMatrix, UpperTriangularMatrix):
            with self.subTest(clase=clase.__name__):
                with self.assertRaises(ValueError, msg="una matriz no cuadrada debería fallar"):
                    clase([[1.0, 0.0, 0.0], [1.0, 1.0, 0.0]])

    def test_lado_equivocado(self):
        with self.subTest(caso="inferior con elemento arriba de la diagonal"):
            with self.assertRaises(ValueError):
                LowerTriangularMatrix(U_DATOS)
        with self.subTest(caso="superior con elemento abajo de la diagonal"):
            with self.assertRaises(ValueError):
                UpperTriangularMatrix(L_DATOS)

    def test_sigue_validando_como_matrix(self):
        for clase in (LowerTriangularMatrix, UpperTriangularMatrix):
            with self.subTest(clase=clase.__name__, caso="vacía"):
                with self.assertRaises(ValueError):
                    clase([])
            with self.subTest(clase=clase.__name__, caso="renglones de distinta longitud"):
                with self.assertRaises(ValueError):
                    clase([[1.0], [1.0, 1.0]])

    def test_diagonal_y_determinante(self):
        L = LowerTriangularMatrix(L_DATOS)
        U = UpperTriangularMatrix(U_DATOS)
        with self.subTest(caso="diagonal"):
            self.assertEqual(L.diagonal(), [2.0, 4.0, 5.0])
        with self.subTest(caso="determinant inferior"):
            self.assertAlmostEqual(L.determinant(), 40.0)
        with self.subTest(caso="determinant superior"):
            self.assertAlmostEqual(U.determinant(), 40.0)

    def test_is_singular(self):
        with self.subTest(caso="no singular"):
            self.assertFalse(LowerTriangularMatrix(L_DATOS).is_singular())
        with self.subTest(caso="singular"):
            singular = UpperTriangularMatrix([[1.0, 2.0], [0.0, 0.0]])
            self.assertTrue(singular.is_singular())

    def test_transpose_cambia_de_tipo(self):
        L = LowerTriangularMatrix(L_DATOS)
        U = UpperTriangularMatrix(U_DATOS)
        with self.subTest(caso="inferior -> superior"):
            self.assertIsInstance(L.transpose(), UpperTriangularMatrix)
        with self.subTest(caso="superior -> inferior"):
            self.assertIsInstance(U.transpose(), LowerTriangularMatrix)
        with self.subTest(caso="datos"):
            self.assertEqual(L.transpose().data, Matrix(L_DATOS).transpose().data)


###############################################
# Ejercicio 2
###############################################


class TestSistemaTriangular(unittest.TestCase):
    def setUp(self):
        requiere_sistema(self)
        self.L = LowerTriangularMatrix(L_DATOS)
        self.U = UpperTriangularMatrix(U_DATOS)

    def assertVectoresCercanos(self, x, y):
        self.assertEqual(len(x), len(y))
        for x_i, y_i in zip(x, y):
            self.assertAlmostEqual(x_i, y_i, places=12)

    def test_singular_error_es_value_error(self):
        self.assertTrue(
            issubclass(MatrizSingularError, ValueError),
            msg="MatrizSingularError debe heredar de ValueError",
        )

    def test_resolver_inferior(self):
        sistema = SistemaTriangular(self.L, producto(L_DATOS, X_EXACTA))
        self.assertVectoresCercanos(sistema.resolver(), X_EXACTA)

    def test_resolver_superior(self):
        sistema = SistemaTriangular(self.U, producto(U_DATOS, X_EXACTA))
        self.assertVectoresCercanos(sistema.resolver(), X_EXACTA)

    def test_metodos_directos(self):
        with self.subTest(caso="sustitucion_adelante"):
            sistema = SistemaTriangular(self.L, producto(L_DATOS, X_EXACTA))
            self.assertVectoresCercanos(sistema.sustitucion_adelante(), X_EXACTA)
        with self.subTest(caso="sustitucion_atras"):
            sistema = SistemaTriangular(self.U, producto(U_DATOS, X_EXACTA))
            self.assertVectoresCercanos(sistema.sustitucion_atras(), X_EXACTA)

    def test_metodo_equivocado(self):
        with self.subTest(caso="adelante con superior"):
            with self.assertRaises(TypeError):
                SistemaTriangular(self.U, [1.0, 1.0, 1.0]).sustitucion_adelante()
        with self.subTest(caso="atrás con inferior"):
            with self.assertRaises(TypeError):
                SistemaTriangular(self.L, [1.0, 1.0, 1.0]).sustitucion_atras()

    def test_residuo(self):
        b = producto(L_DATOS, X_EXACTA)
        sistema = SistemaTriangular(self.L, b)
        r = sistema.residuo(sistema.resolver())
        self.assertVectoresCercanos(r, [0.0, 0.0, 0.0])

    def test_no_modifica_b(self):
        b = [1.0, 2.0, 3.0]
        SistemaTriangular(self.L, b).resolver()
        self.assertEqual(b, [1.0, 2.0, 3.0], msg="resolver() no debe modificar la lista b")

    def test_matriz_no_triangular(self):
        for T in (Matrix(L_DATOS), L_DATOS):
            with self.subTest(tipo=type(T).__name__):
                with self.assertRaises(TypeError):
                    SistemaTriangular(T, [1.0, 1.0, 1.0])

    def test_dimensiones_incompatibles(self):
        with self.assertRaises(ValueError):
            SistemaTriangular(self.L, [1.0, 1.0])

    def test_singular(self):
        casos = [
            ("inferior", LowerTriangularMatrix([[1.0, 0.0], [2.0, 0.0]])),
            ("superior", UpperTriangularMatrix([[0.0, 1.0], [0.0, 2.0]])),
        ]
        for nombre, T in casos:
            with self.subTest(caso=nombre):
                with self.assertRaises(MatrizSingularError):
                    SistemaTriangular(T, [1.0, 1.0]).resolver()

    def test_flops(self):
        for n in (1, 2, 5, 10):
            datos = [[1.0 if j <= i else 0.0 for j in range(n)] for i in range(n)]
            for T in (LowerTriangularMatrix(datos), LowerTriangularMatrix(datos).transpose()):
                with self.subTest(n=n, tipo=type(T).__name__):
                    sistema = SistemaTriangular(T, [1.0] * n)
                    self.assertEqual(sistema.flops, 0, msg="flops debe empezar en 0")
                    sistema.resolver()
                    self.assertEqual(sistema.flops, n**2)
                    sistema.resolver()
                    self.assertEqual(
                        sistema.flops, n**2, msg="flops debe contar solo la última resolución"
                    )


###############################################
# Ejercicio 3
###############################################


class TestComplejidad(unittest.TestCase):
    """Corre complejidad.py una sola vez (setUpClass) y todas las
    pruebas revisan ese mismo resultado."""

    @classmethod
    def setUpClass(cls):
        if not RUTA_SCRIPT.exists():
            cls.resultado = None
            return
        # Borramos los .dat de una corrida anterior, para comprobar
        # que el script en verdad los vuelve a generar.
        RUTA_FLOPS.unlink(missing_ok=True)
        RUTA_TIEMPOS.unlink(missing_ok=True)
        try:
            cls.resultado = subprocess.run(
                [sys.executable, str(RUTA_SCRIPT)],
                cwd=CARPETA_UNIDAD,
                capture_output=True,
                text=True,
                timeout=180,
            )
        except subprocess.TimeoutExpired as error:
            cls.resultado = error

    def _revisar_que_corrio(self):
        if self.resultado is None:
            self.skipTest("unidades/09_algebra_lineal/complejidad.py no existe todavía")
        if isinstance(self.resultado, subprocess.TimeoutExpired):
            self.skipTest("complejidad.py no terminó en 180 s (¿n demasiado grande?)")
        if self.resultado.returncode != 0:
            self.skipTest("el script no corrió; ver test_corre_sin_errores")

    def _leer(self, ruta):
        self._revisar_que_corrio()
        if not ruta.exists():
            self.skipTest(f"datos/{ruta.name} no existe todavía")
        filas = []
        with open(ruta) as archivo:
            for linea in archivo:
                if linea.startswith("#") or not linea.strip():
                    continue
                filas.append([float(c) for c in linea.split()])
        return filas

    def test_corre_sin_errores(self):
        if self.resultado is None:
            self.skipTest("unidades/09_algebra_lineal/complejidad.py no existe todavía")
        if isinstance(self.resultado, subprocess.TimeoutExpired):
            self.fail("complejidad.py no terminó en 180 s (¿n demasiado grande?)")
        self.assertEqual(
            self.resultado.returncode,
            0,
            msg=f"complejidad.py terminó con un error:\n{self.resultado.stderr}",
        )

    def test_conteo_flops(self):
        filas = self._leer(RUTA_FLOPS)
        self.assertGreaterEqual(len(filas), 5, msg="se esperaban al menos 5 valores de n")
        nombres = ["producto_punto", "mat_vec", "mat_mat", "sustitucion"]
        for fila in filas:
            self.assertEqual(
                len(fila), 5, msg="cada línea debe tener 5 columnas: n " + " ".join(nombres)
            )
            n = int(fila[0])
            esperados = [2 * n - 1, 2 * n**2 - n, 2 * n**3 - n**2, n**2]
            for nombre, medido, esperado in zip(nombres, fila[1:], esperados):
                with self.subTest(n=n, operacion=nombre):
                    self.assertEqual(int(medido), esperado)

    def test_formato_tiempos(self):
        filas = self._leer(RUTA_TIEMPOS)
        self.assertGreaterEqual(len(filas), 4, msg="se esperaban al menos 4 valores de n")
        for fila in filas:
            self.assertEqual(
                len(fila), 3, msg="cada línea debe tener 3 columnas: n t_sustitucion t_mat_mat"
            )
        with self.subTest(caso="n se duplica"):
            for anterior, siguiente in zip(filas, filas[1:]):
                self.assertEqual(siguiente[0], 2 * anterior[0])
        with self.subTest(caso="el tiempo crece con n"):
            self.assertGreater(filas[-1][1], filas[0][1])
            self.assertGreater(filas[-1][2], filas[0][2])


###############################################
# Ejercicio 4
###############################################


class CasoEjercicio4(unittest.TestCase):
    """Comparaciones que usan las tres clases de pruebas de abajo."""

    def assertVectoresCercanos(self, x, y):
        self.assertEqual(len(x), len(y))
        for x_i, y_i in zip(x, y):
            self.assertAlmostEqual(x_i, y_i, places=12)

    def assertMatricesCercanas(self, A, B):
        self.assertEqual(len(A), len(B))
        for renglon_A, renglon_B in zip(A, B):
            self.assertVectoresCercanos(renglon_A, renglon_B)


class TestEliminacionGaussiana(CasoEjercicio4):
    def setUp(self):
        requiere_ejercicio_4(self, EliminacionGaussiana, "EliminacionGaussiana")

    def test_pivote_cero_error_es_zero_division_error(self):
        self.assertTrue(
            issubclass(PivoteCeroError, ZeroDivisionError),
            msg="PivoteCeroError debe heredar de ZeroDivisionError",
        )

    def test_solo_U(self):
        e = EliminacionGaussiana(Matrix(A_DATOS))
        with self.subTest(caso="tipo de U"):
            self.assertIsInstance(e.U, UpperTriangularMatrix)
        with self.subTest(caso="U"):
            self.assertMatricesCercanas(e.U.data, U_ELIMINADA)
        with self.subTest(caso="sin b, b_nuevo es None"):
            self.assertIsNone(e.b_nuevo)
        with self.subTest(caso="coeficientes=False, L es None"):
            self.assertIsNone(e.L)

    def test_con_b(self):
        e = EliminacionGaussiana(Matrix(A_DATOS), B_DATOS)
        with self.subTest(caso="U"):
            self.assertMatricesCercanas(e.U.data, U_ELIMINADA)
        with self.subTest(caso="b_nuevo"):
            self.assertVectoresCercanos(e.b_nuevo, B_ELIMINADO)
        with self.subTest(caso="L es None"):
            self.assertIsNone(e.L)

    def test_con_coeficientes(self):
        e = EliminacionGaussiana(Matrix(A_DATOS), coeficientes=True)
        with self.subTest(caso="tipo de L"):
            self.assertIsInstance(e.L, LowerTriangularMatrix)
        with self.subTest(caso="L"):
            self.assertMatricesCercanas(e.L.data, L_COEFICIENTES)
        with self.subTest(caso="U"):
            self.assertMatricesCercanas(e.U.data, U_ELIMINADA)
        with self.subTest(caso="L U = A"):
            self.assertMatricesCercanas((e.L * e.U).data, A_DATOS)

    def test_ceros_de_redondeo(self):
        # Si U[1][0] se calcula como 0.7 - (0.7/0.3)*0.3, queda en -1.1e-16
        # y UpperTriangularMatrix la rechaza.
        try:
            EliminacionGaussiana(Matrix(A_REDONDEO), [1.0, 1.0], coeficientes=True)
        except ValueError:
            self.fail("la U tiene un 'cero' de redondeo abajo de la diagonal (ver practica_04.md)")

    def test_no_modifica_A_ni_b(self):
        A = Matrix(A_DATOS)
        b = list(B_DATOS)
        EliminacionGaussiana(A, b, coeficientes=True)
        with self.subTest(caso="A"):
            self.assertEqual(A.data, A_DATOS, msg="no se debe modificar A")
        with self.subTest(caso="b"):
            self.assertEqual(b, B_DATOS, msg="no se debe modificar b")

    def test_entradas_invalidas(self):
        with self.subTest(caso="lista de listas"):
            with self.assertRaises(TypeError):
                EliminacionGaussiana(A_DATOS)
        with self.subTest(caso="no cuadrada"):
            with self.assertRaises(ValueError):
                EliminacionGaussiana(Matrix([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]))
        with self.subTest(caso="b de otro tamaño"):
            with self.assertRaises(ValueError):
                EliminacionGaussiana(Matrix(A_DATOS), [1.0, 1.0])

    def test_pivote_cero(self):
        for coeficientes in (False, True):
            with self.subTest(coeficientes=coeficientes):
                with self.assertRaises(PivoteCeroError):
                    EliminacionGaussiana(Matrix(A_PIVOTE_CERO), [1.0, 2.0], coeficientes)

    def test_singular_si_termina(self):
        # El pivote cero queda en el último renglón: no se divide entre él.
        e = EliminacionGaussiana(Matrix(A_SINGULAR), [1.0, 1.0])
        self.assertEqual(e.U.data[1][1], 0.0)


class TestSistemaLineal(CasoEjercicio4):
    def setUp(self):
        requiere_ejercicio_4(self, SistemaLineal, "SistemaLineal")

    def test_resolver_gauss(self):
        sistema = SistemaLineal(Matrix(A_DATOS), B_DATOS)
        self.assertVectoresCercanos(sistema.resolver_gauss(), X_DATOS)

    def test_acepta_triangulares(self):
        # Una LowerTriangularMatrix es una Matrix (principio L de SOLID).
        sistema = SistemaLineal(LowerTriangularMatrix(L_DATOS), producto(L_DATOS, X_EXACTA))
        self.assertVectoresCercanos(sistema.resolver_gauss(), X_EXACTA)

    def test_ceros_de_redondeo(self):
        x_exacta = [1.0, 2.0]
        sistema = SistemaLineal(Matrix(A_REDONDEO), producto(A_REDONDEO, x_exacta))
        try:
            x = sistema.resolver_gauss()
        except ValueError:
            self.fail("la U tiene un 'cero' de redondeo abajo de la diagonal (ver practica_04.md)")
        self.assertVectoresCercanos(x, x_exacta)

    def test_residuo(self):
        sistema = SistemaLineal(Matrix(A_DATOS), B_DATOS)
        r = sistema.residuo(sistema.resolver_gauss())
        self.assertVectoresCercanos(r, [0.0, 0.0, 0.0])

    def test_no_modifica_A_ni_b(self):
        A = Matrix(A_DATOS)
        b = list(B_DATOS)
        SistemaLineal(A, b).resolver_gauss()
        with self.subTest(caso="A"):
            self.assertEqual(A.data, A_DATOS, msg="no se debe modificar A")
        with self.subTest(caso="b"):
            self.assertEqual(b, B_DATOS, msg="no se debe modificar b")

    def test_entradas_invalidas(self):
        with self.subTest(caso="lista de listas"):
            with self.assertRaises(TypeError):
                SistemaLineal(A_DATOS, B_DATOS)
        with self.subTest(caso="no cuadrada"):
            with self.assertRaises(ValueError):
                SistemaLineal(Matrix([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]), [1.0, 1.0])
        with self.subTest(caso="b de otro tamaño"):
            with self.assertRaises(ValueError):
                SistemaLineal(Matrix(A_DATOS), [1.0, 1.0])

    def test_pivote_cero(self):
        with self.assertRaises(PivoteCeroError):
            SistemaLineal(Matrix(A_PIVOTE_CERO), [1.0, 2.0]).resolver_gauss()

    def test_singular(self):
        with self.assertRaises(MatrizSingularError):
            SistemaLineal(Matrix(A_SINGULAR), [1.0, 1.0]).resolver_gauss()


class TestFactorizacionLU(CasoEjercicio4):
    def setUp(self):
        requiere_ejercicio_4(self, FactorizacionLU, "FactorizacionLU")

    def test_L_y_U(self):
        lu = FactorizacionLU(Matrix(A_DATOS))
        with self.subTest(caso="tipo de L"):
            self.assertIsInstance(lu.L, LowerTriangularMatrix)
        with self.subTest(caso="tipo de U"):
            self.assertIsInstance(lu.U, UpperTriangularMatrix)
        with self.subTest(caso="L"):
            self.assertMatricesCercanas(lu.L.data, L_COEFICIENTES)
        with self.subTest(caso="L U = A"):
            self.assertMatricesCercanas((lu.L * lu.U).data, A_DATOS)
        with self.subTest(caso="ceros de redondeo"):
            lu_redondeo = FactorizacionLU(Matrix(A_REDONDEO))
            self.assertMatricesCercanas((lu_redondeo.L * lu_redondeo.U).data, A_REDONDEO)

    def test_resolver_varios_b(self):
        lu = FactorizacionLU(Matrix(A_DATOS))
        with self.subTest(caso="b del ejemplo"):
            self.assertVectoresCercanos(lu.resolver(B_DATOS), X_DATOS)
        with self.subTest(caso="otro b, misma factorización"):
            x_otra = [1.0, 2.0, 3.0]
            self.assertVectoresCercanos(lu.resolver(producto(A_DATOS, x_otra)), x_otra)

    def test_resolver_lu_en_sistema_lineal(self):
        if SistemaLineal is None:
            self.skipTest("todavía no existe SistemaLineal (Ejercicio 4b)")
        sistema = SistemaLineal(Matrix(A_DATOS), B_DATOS)
        self.assertVectoresCercanos(sistema.resolver_lu(), X_DATOS)

    def test_determinante(self):
        with self.subTest(caso="ejemplo de notas.md"):
            self.assertAlmostEqual(FactorizacionLU(Matrix(A_DATOS)).determinante(), 40.0)
        with self.subTest(caso="matriz con ceros de redondeo"):
            self.assertAlmostEqual(FactorizacionLU(Matrix(A_REDONDEO)).determinante(), -0.1)
        with self.subTest(caso="singular"):
            self.assertAlmostEqual(FactorizacionLU(Matrix(A_SINGULAR)).determinante(), 0.0)

    def test_inversa(self):
        A = Matrix(A_DATOS)
        inv = FactorizacionLU(A).inversa()
        with self.subTest(caso="tipo"):
            self.assertIsInstance(inv, Matrix)
        with self.subTest(caso="A A^-1 = I"):
            identidad = [[1.0 if i == j else 0.0 for j in range(3)] for i in range(3)]
            self.assertMatricesCercanas((A * inv).data, identidad)

    def test_entradas_invalidas(self):
        with self.subTest(caso="lista de listas"):
            with self.assertRaises(TypeError):
                FactorizacionLU(A_DATOS)
        with self.subTest(caso="no cuadrada"):
            with self.assertRaises(ValueError):
                FactorizacionLU(Matrix([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]))

    def test_pivote_cero(self):
        with self.assertRaises(PivoteCeroError):
            FactorizacionLU(Matrix(A_PIVOTE_CERO))

    def test_singular(self):
        lu = FactorizacionLU(Matrix(A_SINGULAR))
        with self.subTest(caso="resolver"):
            with self.assertRaises(MatrizSingularError):
                lu.resolver([1.0, 1.0])
        with self.subTest(caso="inversa"):
            with self.assertRaises(MatrizSingularError):
                lu.inversa()


if __name__ == "__main__":
    unittest.main(verbosity=2)
