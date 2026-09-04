#!/usr/bin/env python3
"""Pruebas para la Práctica 1 (ver practica_01.md), con `unittest`.

Se corre con `python3 practicas/pruebas_practica_01.py`, con el
entorno virtual activado (`source .venv/bin/activate` desde la raíz
del repositorio) y `fiscomp` instalado en modo editable
(`pip install -e .`; ver README.md).

Cada valor de prueba se revisa con `subTest`, así que si uno falla los
demás se siguen corriendo. Los mensajes de falla solo dicen qué caso
falló (la función y el valor de `x`), no el valor esperado ni el
obtenido — esa parte les toca investigarla a ustedes.
"""

import math
import subprocess
import sys
import unittest
from pathlib import Path

from fiscomp import funciones_especiales as fe
from fiscomp.precision_numerica import error_relativo

TOLERANCIA = 1e-9


class PruebaBase(unittest.TestCase):
    # Con longMessage=False, el mensaje que le pasemos a assertTrue
    # reemplaza el mensaje por default (que sí incluiría los valores).
    longMessage = False

    def revisar(self, nombre_funcion, funcion_real, valores_de_prueba, tolerancia=TOLERANCIA):
        funcion_bajo_prueba = getattr(fe, nombre_funcion, None)
        if funcion_bajo_prueba is None:
            self.skipTest(f"{nombre_funcion}() no está definida todavía en funciones_especiales.py")

        for x in valores_de_prueba:
            with self.subTest(x=x):
                aproximado = funcion_bajo_prueba(x)
                exacto = funcion_real(x)
                self.assertTrue(
                    error_relativo(aproximado, exacto) < tolerancia,
                    msg=f"{nombre_funcion}({x!r}) no es suficientemente preciso",
                )


class TestFactorialYSeno(PruebaBase):
    """Ya resueltas en clase: si esto falla, revisen su entorno (¿está
    activado el entorno virtual? ¿instalaron fiscomp con pip install -e .?)
    antes que su propio código."""

    def test_factorial(self):
        self.revisar("factorial", math.factorial, [0, 1, 5, 10])

    def test_seno(self):
        self.revisar("seno", math.sin, [0.0, 0.5, 1.0, -1.0, 2.0])


class TestCoseno(PruebaBase):
    """Ejercicio 3.

    Los valores de prueba evitan a propósito quedar muy cerca de una
    raíz de coseno (pi/2, 3*pi/2, ...): ahí el valor "exacto" está muy
    cerca de 0, y error_relativo() (que divide entre el valor exacto)
    se dispara aunque el error absoluto sea chiquito. Ese fenómeno es,
    de hecho, el tema del Ejercicio 4.
    """

    def test_coseno(self):
        self.revisar("coseno", math.cos, [0.0, 0.5, 1.0, -1.0, 2.0, -2.5, 3.0, 4.0])


class TestExponencial(PruebaBase):
    """Ejercicio 3."""

    def test_exponencial(self):
        self.revisar("exponencial", math.exp, [0.0, 0.1, 1.0, -1.0, -3.0, 2.5, 5.0])


class TestLn(PruebaBase):
    """Ejercicio 3.

    Incluye valores fuera de (0, 2], donde no converge la serie de
    Taylor de ln(x) alrededor de x=1 (ver la pista del Ejercicio 3):
    si su `ln()` solo funciona en ese intervalo, estos casos van a
    fallar y es la señal de que les falta la serie alternativa.
    """

    def test_ln(self):
        self.revisar("ln", math.log, [1.0, 0.1, 0.5, 2.0, 5.0, 10.0])


class TestPi(PruebaBase):
    """Ejercicio 2."""

    def test_pi_guardado(self):
        try:
            from fiscomp.constantes import PI
        except ImportError:
            self.skipTest("fiscomp/constantes.py no existe o no define PI")

        self.assertTrue(
            error_relativo(PI, math.pi) < 1e-4,
            msg="PI: el valor guardado está lejos del pi real",
        )


class TestRecoleccionDatos(unittest.TestCase):
    """Ejercicio 1.

    A diferencia de las demás pruebas, esta no revisa el contenido del
    reporte (recoleccion_datos.py pide datos por input(), y qué tan
    detallado sea el reporte queda a criterio de cada quien). Solo
    comprueba que el script corra de principio a fin sin crashear, y
    que genere reporte_recoleccion.txt. Le manda, por la entrada
    estándar, un montón de líneas genéricas (texto y números) para que
    no se quede esperando una respuesta que nunca llega, sin importar
    cuántas veces llame a input().
    """

    CARPETA = Path(__file__).resolve().parent
    RUTA_SCRIPT = CARPETA / "recoleccion_datos.py"
    RUTA_REPORTE = CARPETA / "reporte_recoleccion.txt"

    def test_corre_sin_errores_y_genera_el_reporte(self):
        if not self.RUTA_SCRIPT.exists():
            self.skipTest("recoleccion_datos.py no existe todavía")

        # Borramos cualquier reporte de una corrida anterior, para
        # comprobar que este script en verdad lo vuelve a generar
        # (y no que ya estaba ahí de antes).
        self.RUTA_REPORTE.unlink(missing_ok=True)

        # "1" sirve como respuesta genérica sin importar qué pida el
        # prompt: es un str válido (para nombre/responsable), y se
        # puede convertir con int("1") y con float("1") sin tronar.
        # No puede alternar con otro valor (como un texto no numérico)
        # porque no sabemos en qué orden pide cada dato cada quien.
        entradas_de_prueba = "\n".join(["1"] * 30) + "\n"
        resultado = subprocess.run(
            [sys.executable, str(self.RUTA_SCRIPT)],
            cwd=self.CARPETA,  # para que "reporte_recoleccion.txt" caiga aquí
            input=entradas_de_prueba,
            capture_output=True,
            text=True,
            timeout=10,
        )
        self.assertEqual(
            resultado.returncode,
            0,
            msg=f"recoleccion_datos.py terminó con un error:\n{resultado.stderr}",
        )
        self.assertTrue(
            self.RUTA_REPORTE.exists(),
            msg="recoleccion_datos.py no generó reporte_recoleccion.txt",
        )
        self.assertTrue(
            self.RUTA_REPORTE.read_text().strip(),
            msg="reporte_recoleccion.txt se generó pero está vacío",
        )


class TestReporteEjercicio4(unittest.TestCase):
    """Ejercicio 4.

    Tampoco revisa contenido (es un reporte abierto: valores de x,
    conclusiones...), solo que el archivo exista y no esté vacío. Que
    esta prueba pase no quiere decir que las respuestas estén bien —
    esa parte se las revisamos a mano.
    """

    RUTA_REPORTE = Path(__file__).resolve().parent / "reporte_ejercicio4.txt"

    def test_reporte_existe(self):
        if not self.RUTA_REPORTE.exists():
            self.skipTest("reporte_ejercicio4.txt no existe todavía")

        contenido = self.RUTA_REPORTE.read_text().strip()
        self.assertTrue(contenido, msg="reporte_ejercicio4.txt existe pero está vacío")


if __name__ == "__main__":
    unittest.main(verbosity=2)
