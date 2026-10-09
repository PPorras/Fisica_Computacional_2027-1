# Solución de la Práctica 1

Una solución posible de la [Práctica 1](../practica_01.md). No es la
única: en los ejercicios 1 y 4 (reportes abiertos) cualquier
respuesta bien justificada es válida.

| Ejercicio | Archivo | Dónde va en la práctica |
|---|---|---|
| 1 (Millikan) | [`recoleccion_datos.py`](recoleccion_datos.py) | `practicas/practica1/` |
| 2 (pi con Leibniz) | [`calcular_pi.py`](calcular_pi.py) | `practicas/practica1/` |
| 2 (la constante guardada) | [`constantes.py`](constantes.py) | `fiscomp/` |
| 3 (`coseno`, `exponencial`, `ln`) | [`fiscomp/funciones_especiales.py`](../../../fiscomp/funciones_especiales.py) | ya está en el repositorio |
| 4 (errores) | [`comparar_errores.py`](comparar_errores.py), que genera [`reporte_ejercicio4.txt`](reporte_ejercicio4.txt) | `practicas/practica1/` |

La solución del Ejercicio 3 es la versión de `funciones_especiales.py`
que ya está en `fiscomp`. Esa versión es más completa que lo que pedía
la práctica: incluye la reducción de rango en `seno` y `coseno` (el
reto extra del Ejercicio 3), y calcula $e^x$ para $x < 0$ como
$`1/e^{|x|}`$. El reporte del Ejercicio 4 explica por qué eso último
importa.

## Cómo correrla

Desde la raíz del repositorio, con el entorno virtual activado:

```bash
python3 practicas/practica1/solucion_practica_01/calcular_pi.py        # ~4 s
python3 practicas/practica1/solucion_practica_01/comparar_errores.py
python3 practicas/practica1/solucion_practica_01/recoleccion_datos.py  # pide datos con input()
```

Para pasarla por las pruebas de la práctica, copien cada archivo al
lugar de la última columna de la tabla y corran
`python3 practicas/practica1/pruebas_practica_01.py`.
