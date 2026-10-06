#!/usr/bin/env python3
"""Pruebas para la Práctica 3 (ver practica_03.md), con `unittest`.

Se corre con `python3 practicas/practica3/pruebas_practica_03.py`, con el
entorno virtual activado (`source .venv/bin/activate` desde la raíz
del repositorio) y `fiscomp` instalado en modo editable
(`pip install -e .`; ver README.md).

Corre diferencias_finitas.py completo, como subproceso y con un
límite de tiempo (por si el barrido de h de alguien se queda pegado
en un ciclo infinito), y revisa el archivo datos/derivada_sin_x2.dat
que debe generar (Ejercicio 3): que exista, que tenga el formato de
columnas esperado, y que el error de diff_central sea
consistentemente más chico que el de diff_forward para h grande (ahí
domina el error de truncamiento, y eso sí tiene una respuesta
objetivamente correcta: O(h^2) < O(h) si diff_central está bien
implementada).

No revisa diff_backward por separado (Ejercicio 1), ni el contenido
de la tabla del Ejercicio 2, ni la respuesta del Ejercicio 4: esas
partes son más abiertas y se revisan a mano.
"""

import subprocess
import sys
import unittest
from pathlib import Path

CARPETA_UNIDAD = (
    Path(__file__).resolve().parent.parent.parent / "unidades" / "08_diferencias_finitas"
)
RUTA_SCRIPT = CARPETA_UNIDAD / "diferencias_finitas.py"
RUTA_DATOS = CARPETA_UNIDAD / "datos" / "derivada_sin_x2.dat"


class TestBarridoDeH(unittest.TestCase):
    """Corre diferencias_finitas.py una sola vez (setUpClass) y todas
    las pruebas revisan ese mismo resultado, en vez de volver a correr
    el barrido de h en cada una."""

    @classmethod
    def setUpClass(cls):
        if not RUTA_SCRIPT.exists():
            cls.resultado = None
            return

        # Borramos cualquier .dat de una corrida anterior, para
        # comprobar que el script en verdad lo vuelve a generar.
        RUTA_DATOS.unlink(missing_ok=True)

        try:
            cls.resultado = subprocess.run(
                [sys.executable, str(RUTA_SCRIPT)],
                cwd=CARPETA_UNIDAD,
                capture_output=True,
                text=True,
                timeout=60,
            )
        except subprocess.TimeoutExpired as error:
            cls.resultado = error

    def _filas_o_skip(self):
        if self.resultado is None:
            self.skipTest("diferencias_finitas.py no existe todavía")
        if isinstance(self.resultado, subprocess.TimeoutExpired):
            self.skipTest("diferencias_finitas.py no terminó en 60s (¿un ciclo infinito?)")
        if self.resultado.returncode != 0:
            self.skipTest("el script no corrió; ver test_corre_sin_errores")
        if not RUTA_DATOS.exists():
            self.skipTest("datos/derivada_sin_x2.dat no existe todavía (Ejercicio 3)")

        filas = []
        with open(RUTA_DATOS) as archivo:
            for linea in archivo:
                if linea.startswith("#") or not linea.strip():
                    continue
                filas.append([float(c) for c in linea.split()])
        return filas

    def test_corre_sin_errores(self):
        if self.resultado is None:
            self.skipTest("diferencias_finitas.py no existe todavía")
        if isinstance(self.resultado, subprocess.TimeoutExpired):
            self.fail("diferencias_finitas.py no terminó en 60s (¿un ciclo infinito?)")
        self.assertEqual(
            self.resultado.returncode,
            0,
            msg=f"diferencias_finitas.py terminó con un error:\n{self.resultado.stderr}",
        )

    def test_genera_el_archivo_de_datos(self):
        if self.resultado is None or isinstance(self.resultado, subprocess.TimeoutExpired):
            self.skipTest("ver test_corre_sin_errores")
        if self.resultado.returncode != 0:
            self.skipTest("el script no corrió; ver test_corre_sin_errores")
        self.assertTrue(
            RUTA_DATOS.exists(),
            msg="no se generó datos/derivada_sin_x2.dat (Ejercicio 3)",
        )

    def test_formato_de_columnas(self):
        filas = self._filas_o_skip()
        self.assertGreaterEqual(
            len(filas),
            10,
            msg="se esperaban al menos 10 valores de h en el barrido",
        )
        for fila in filas[:5]:
            with self.subTest(fila=fila):
                self.assertEqual(
                    len(fila),
                    5,
                    msg=(
                        "cada línea debe tener 5 columnas: "
                        "h diff_forward error_forward diff_central error_central"
                    ),
                )

    def test_central_es_mas_precisa_que_adelante_para_h_grande(self):
        filas = self._filas_o_skip()
        if not filas:
            self.skipTest("datos/derivada_sin_x2.dat está vacío")

        for fila in filas[:3]:
            h, _, error_adelante, _, error_central = fila
            with self.subTest(h=h):
                self.assertLess(
                    error_central,
                    error_adelante,
                    msg=f"con h={h}, error_central debería ser menor que error_adelante",
                )

    def test_el_error_baja_al_reducir_h_al_principio(self):
        """Mientras h sea 'grande' (primeras filas, antes de que el
        error de redondeo tome el control), reducir h a la mitad debe
        bajar el error -- es la firma del error de truncamiento."""
        filas = self._filas_o_skip()
        if len(filas) < 5:
            self.skipTest("no hay suficientes filas para revisar la tendencia")

        errores_centrales = [fila[4] for fila in filas[:5]]
        self.assertEqual(
            errores_centrales,
            sorted(errores_centrales, reverse=True),
            msg="error_central debería ir bajando en las primeras filas del barrido (h grande)",
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
