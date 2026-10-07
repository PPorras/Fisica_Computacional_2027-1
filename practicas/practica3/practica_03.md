# Práctica 3 — Diferencias finitas

Cubre lo visto en la unidad
[08 (diferencias finitas)](../../unidades/08_diferencias_finitas/).

## Contexto

En
[`unidades/08_diferencias_finitas/diferencias_finitas.py`](../../unidades/08_diferencias_finitas/diferencias_finitas.py)
ya está implementada `diff_forward(f, x0, h)`, la diferencia hacia
adelante:

```python
def diff_forward(f, x0, h):
    return (f(x0 + h) - f(x0)) / h
```

junto con las funciones de prueba (`const_5`, `ident`, `sqr`,
`sin_x2`) y un ejemplo de uso. El resto de la unidad — derivado a
partir de la serie de Taylor en
[`notas.md`](../../unidades/08_diferencias_finitas/notas.md) — lo
completan ustedes, en ese mismo archivo.

## Ejercicio 1 — `diff_backward` y `diff_central`

Agreguen, junto a `diff_forward`:

- `diff_backward(f, x0, h)`: diferencia hacia atrás,
  `(f(x0) - f(x0 - h)) / h`.
- `diff_central(f, x0, h)`: diferencia central,
  `(f(x0 + h/2) - f(x0 - h/2)) / h`.

Ambas fórmulas están deducidas paso a paso en `notas.md`, a partir de
la serie de Taylor de $f$ alrededor de $x_0$.

## Ejercicio 2 — derivadas exactas, para comparar

Para poder medir el error de sus tres aproximaciones necesitan la
derivada *exacta* (analítica) de cada función de prueba. Agreguen:

- `const_5_prima(x)` → `0`
- `ident_prima(x)` → `1`
- `sqr_prima(x)` → `2*x`
- `sin_x2_prima(x)` → `2*x*coseno(x**2)` (regla de la cadena; usen
  `coseno` de `fiscomp.funciones_especiales`, igual que `seno` ya se
  usa en `sin_x2`)

Con `error_relativo()` (de `fiscomp.precision_numerica`), impriman una
tabla que compare, para cada función de prueba y con una `h` fija
(por ejemplo `h=0.1` en `x0=6.0`), `diff_forward`, `diff_backward` y
`diff_central` contra la derivada exacta.

Revisen que para `const_5` y `ident` el error sea (esencialmente)
cero con los tres métodos: la serie de Taylor de una función
constante o lineal no tiene nada que truncar. Para `sqr` y
`sin_x2` sí debería aparecer un error, más chico en `diff_central`
que en `diff_forward`/`diff_backward`.

## Ejercicio 3 — barrido de h: truncamiento vs. redondeo

Para `f(x) = sin(x**2)` en `x0 = 1.0`, hagan un barrido de `h`
empezando en `h = 1.0` y dividiendo entre 2 repetidamente (unas 50
veces), calculando en cada paso `diff_forward`, `diff_central` y su
error relativo contra la derivada exacta.

Guarden el resultado en un archivo `datos/derivada_sin_x2.dat`
(creen la carpeta `datos/` junto al script si no existe — usen
`Path(__file__).resolve().parent / "datos"`, como en la unidad 04, en
vez de una ruta relativa a mano, para que funcione sin importar desde
dónde corran el script), con una línea de encabezado que empiece con
`#` y luego una línea por cada `h`, con **cinco columnas separadas
por espacios**:

```
h  diff_forward  error_forward  diff_central  error_central
```

Este es el formato que ya esperan
[`graficar_derivada.gp`](../../unidades/08_diferencias_finitas/graficar_derivada.gp)
y
[`graficar_derivada.py`](../../unidades/08_diferencias_finitas/graficar_derivada.py)
(Ejercicio 5), así que no hace falta tocar esos dos archivos.

## Ejercicio 4 — ¿dónde está el h óptimo?

Con los datos del Ejercicio 3, encuentren **con código** (no a ojo)
la `h` que da el menor error para cada método, y compárenla contra la
fórmula de `notas.md` (que asume que $f$ y sus derivadas son de orden
1):

```
h_opt_adelante ~ sqrt(4 * EPS)
h_opt_central  ~ (24 * EPS) ** (1/3)
```

(`EPS` es `fiscomp.precision_numerica.EPS`). Impriman ambos valores —
el que midieron y el de la fórmula — y respondan en un comentario:
¿coinciden exactamente? ¿Por qué no tendrían que coincidir
exactamente? (Piensen en que la fórmula asume $f(x_0)$ y sus
derivadas de orden 1 —¿lo son, para $\sin(x^2)$ en $x_0=1$?— y en que
el barrido solo prueba potencias de $1/2$, no cualquier valor de
$h$.)

## Ejercicio 5 (opcional) — graficar

Una vez generado el `.dat` del Ejercicio 3, corran, desde la carpeta
de la unidad:

```bash
gnuplot graficar_derivada.gp
```

o, si prefieren matplotlib (instrucciones de instalación al inicio de
ese archivo, con el entorno virtual activado):

```bash
python3 graficar_derivada.py
```

Ambos grafican el error contra `h` en escala log-log, con una línea
vertical punteada en la `h_opt` teórica de cada método. Comparen esa
línea contra el punto donde su curva realmente da vuelta (ver la
discusión de por qué no coinciden exactamente, del Ejercicio 4).

## Cómo revisar su trabajo

Hay un script de pruebas en
[`pruebas_practica_03.py`](pruebas_practica_03.py) que corre
`diferencias_finitas.py` completo y revisa el archivo
`datos/derivada_sin_x2.dat` que debe generar (Ejercicio 3): que
exista, que tenga el formato de columnas esperado, y que el error de
`diff_central` sea consistentemente más chico que el de
`diff_forward` para `h` grande — una consecuencia directa de que
`diff_central` esté bien implementada ($`O(h^2) < O(h)`$), no un juicio
de estilo. Córranlo desde la raíz del repositorio:

```bash
python3 practicas/practica3/pruebas_practica_03.py
```

No revisa `diff_backward` por separado, ni el contenido exacto de la
tabla del Ejercicio 2, ni la respuesta del Ejercicio 4: esas partes
son más abiertas y se revisan a mano.
