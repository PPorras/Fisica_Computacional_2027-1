# Solución de la Práctica 3

Una solución posible de la [Práctica 3](../practica_03.md):
[`diferencias_finitas.py`](diferencias_finitas.py), la versión
completa de
[`unidades/08_diferencias_finitas/diferencias_finitas.py`](../../../unidades/08_diferencias_finitas/diferencias_finitas.py).

| Ejercicio | Dónde está en el archivo |
|---|---|
| 1 (`diff_backward`, `diff_central`) | sección "Las tres aproximaciones de diferencias finitas" |
| 2 (derivadas exactas y tabla de comparación) | secciones "Funciones de prueba y sus derivadas exactas" y "Comparar los tres métodos con una h fija" |
| 3 (barrido de `h`) | sección "Barrido de h: error de truncamiento vs. error de redondeo" |
| 3 (`datos/derivada_sin_x2.dat`) | sección "Guardar el barrido en un archivo, para graficarlo" |
| 4 (`h` óptimo medido contra la fórmula) | al final del barrido, con la respuesta en comentarios ("Respuesta del Ejercicio 4") |

Además trae cosas que la práctica no pedía:

- `h` óptimo por defecto en cada método (`h_optimo_adelante`,
  `h_optimo_central`);
- los métodos más precisos de la sección "Diferencias finitas más
  precisas" de `notas.md` (`diff_forward2`, `diff_central2`), que
  generan `datos/derivada_sin_x2_precisas.dat`;
- la respuesta del "Para pensar" de `notas.md` (extrapolación de
  Richardson).

## Cómo correrla

Copien el archivo encima de
`unidades/08_diferencias_finitas/diferencias_finitas.py` (guarden
antes el suyo) y córranlo desde la raíz del repositorio, con el
entorno virtual activado:

```bash
python3 unidades/08_diferencias_finitas/diferencias_finitas.py
python3 practicas/practica3/pruebas_practica_03.py
```

Así los datos quedan en `unidades/08_diferencias_finitas/datos/`, que
es donde los buscan `graficar_derivada.gp` y `graficar_derivada.py`
(Ejercicio 5).
