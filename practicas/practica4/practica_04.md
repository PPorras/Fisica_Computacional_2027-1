# Práctica 4 — Matrices triangulares, sistemas triangulares y complejidad

Cubre lo visto en la unidad
[09 (álgebra lineal)](../../unidades/09_algebra_lineal/), secciones
"Matrices triangulares" y "Conteo de operaciones" de
[`notas.md`](../../unidades/09_algebra_lineal/notas.md), y retoma la
clase `Matrix` de la [Práctica 2](../practica2/practica_02.md).

## Contexto

Hasta ahora tenemos dos piezas sueltas:

- En [`fiscomp/matrices.py`](../../fiscomp/matrices.py), la clase
  `Matrix`, que sirve para *cualquier* matriz rectangular.
- En [`fiscomp/algebra_lineal.py`](../../fiscomp/algebra_lineal.py),
  las funciones `sustitucion_adelante(L, b)` y `sustitucion_atras(U, b)`,
  que reciben una lista de listas y **confían** en que sea triangular y
  sin ceros en la diagonal. Si se les pasa otra cosa, o regresan
  basura sin avisar, o truenan con un `ZeroDivisionError` que no dice
  qué pasó.

En esta práctica vamos a juntar las dos piezas: una matriz triangular
va a ser un *tipo* de `Matrix` que **garantiza** ser triangular desde
que se crea (igual que `Intervalo` garantiza `lo <= hi` en la
unidad 07), y la solución de un sistema triangular va a ser una clase
que valida lo que recibe y explica con claridad cuando algo sale mal.
Al final, medimos cuánto cuesta cada operación.

No borren `sustitucion_adelante` ni `sustitucion_atras`: las usa
[`triangulares.py`](../../unidades/09_algebra_lineal/triangulares.py)
y las vamos a reutilizar en eliminación gaussiana. Lo que escriban en
esta práctica va **junto** a ellas.

## Ejercicio 1 — Clases hijas de `Matrix`: matrices triangulares

### Herencia y `super()`

En la unidad 07 mencionamos la **herencia**: `class Hija(Padre):`
define una clase que tiene automáticamente todos los atributos y
métodos de `Padre`, y que puede agregar métodos nuevos o
*sobreescribir* (redefinir) algunos. Una matriz triangular **es una**
matriz, así que no tiene sentido volver a escribir `shape()`,
`get_row()`, `+`, `*`, etc.: los heredamos de `Matrix`.

Lo que sí cambia es la validación al crearla. Para no repetir la que
ya hace `Matrix.__init__` (que no esté vacía y que sea rectangular),
desde el `__init__` de la clase hija se llama al de la clase padre con
`super()`:

```python
class LowerTriangularMatrix(Matrix):
    def __init__(self, data):
        super().__init__(data)   # corre Matrix.__init__: llena self.data, self.rows, self.cols
        # ... y aquí van las validaciones nuevas
```

### Lo que tienen que hacer

En `fiscomp/matrices.py`, después de la clase `Matrix`, definan dos
clases hijas (usamos nombres en inglés, como `Matrix` y sus métodos):

- `LowerTriangularMatrix(data)`: triangular inferior,
  `data[i][j] == 0` para `j > i`.
- `UpperTriangularMatrix(data)`: triangular superior,
  `data[i][j] == 0` para `j < i`.

Su `__init__` debe levantar `ValueError`, con un mensaje que diga qué
está mal, si:

- la matriz no es cuadrada;
- hay un elemento distinto de cero del lado equivocado de la diagonal
  (digan en el mensaje en qué posición está).

Los casos de matriz vacía o no rectangular ya los atrapa
`Matrix.__init__`, si llaman a `super().__init__(data)` **antes** de
sus propias validaciones (piensen por qué tiene que ser antes).

Además, ambas clases deben tener estos métodos:

| Método | Regresa |
|---|---|
| `diagonal()` | la lista `[A_00, A_11, ..., A_(n-1)(n-1)]` |
| `determinant()` | el determinante: para una matriz triangular, el producto de la diagonal |
| `is_singular()` | `True` si algún elemento de la diagonal es cero |
| `transpose()` | la transpuesta, **como matriz triangular del otro tipo** |

El último es una *sobreescritura*: `Matrix.transpose()` regresa una
`Matrix` común, pero la transpuesta de una triangular inferior es
triangular superior, y viceversa. Pueden reutilizar el método del padre
con `super().transpose()` y construir con su `.data` la clase que
corresponde.

**Sugerencia (no obligatoria).** `diagonal()`, `determinant()`,
`is_singular()` y la revisión de que sea cuadrada son *idénticos* en
las dos clases. En vez de escribirlos dos veces, pueden definir una
clase intermedia `TriangularMatrix(Matrix)` con lo común, y luego
`LowerTriangularMatrix(TriangularMatrix)` y
`UpperTriangularMatrix(TriangularMatrix)` con lo que cambia (cuál
lado de la diagonal revisar, y `transpose()`). Las pruebas solo piden
que ambas hereden, directa o indirectamente, de `Matrix`.

Escriban docstrings para cada clase y cada método, con el mismo
formato que los de `Matrix` (secciones `Parameters`, `Returns`,
`Raises`).

**Para probar a mano:** con

```python
L = LowerTriangularMatrix([[2, 0, 0], [1, 3, 0], [4, 5, 6]])
```

`L.diagonal()` es `[2, 3, 6]`, `L.determinant()` es `36`,
`isinstance(L, Matrix)` es `True`, `L.transpose()` es una
`UpperTriangularMatrix`, y `LowerTriangularMatrix([[1, 2], [0, 1]])`
debe levantar `ValueError` (el `2` está arriba de la diagonal).

## Ejercicio 2 — La clase `SistemaTriangular`

En `fiscomp/algebra_lineal.py`, al final del archivo, conviertan la
sustitución hacia adelante y hacia atrás en una clase que represente
el sistema $T\mathbf{x} = \mathbf{b}$, con $T$ triangular. (Necesitan
importar las dos clases del Ejercicio 1 al inicio del archivo:
`from fiscomp.matrices import LowerTriangularMatrix, UpperTriangularMatrix`.)

### Una excepción propia

Primero definan, como en la unidad 05, una excepción para el caso en
que la matriz es singular:

```python
class MatrizSingularError(ValueError):
    """..."""
```

Hereda de `ValueError` (no directamente de `Exception`) para que
quien ya atrape `ValueError` la atrape también: una matriz singular
*es* un valor inválido.

### La clase

`SistemaTriangular(T, b)` debe tener:

- **`__init__(self, T, b)`**: guarda `self.T`, `self.b` (una *copia*
  de la lista `b`, para que resolver el sistema nunca modifique la
  lista original de quien lo llama), `self.n` (el tamaño) y
  `self.flops = 0`. Levanta:
  - `TypeError` si `T` no es `LowerTriangularMatrix` ni
    `UpperTriangularMatrix` (por ejemplo, si es una `Matrix` común o
    una lista de listas). Usen `isinstance`.
  - `ValueError` si `b` no tiene tantas componentes como renglones
    tiene `T`.
- **`sustitucion_adelante(self)`**: la sustitución hacia adelante;
  regresa `x` como lista. Levanta `TypeError` si `T` no es triangular
  inferior.
- **`sustitucion_atras(self)`**: la sustitución hacia atrás; regresa
  `x`. Levanta `TypeError` si `T` no es triangular superior.
- **`resolver(self)`**: llama a `sustitucion_adelante` o a
  `sustitucion_atras` según el tipo de `T`. Quien usa la clase ya no
  tiene que acordarse de cuál de las dos toca.
- **`residuo(self, x)`**: el vector `b - T x`. Reutilicen la función
  `residuo` que ya está en el archivo (recuerden que recibe listas de
  listas: `self.T.data`).

Pueden partir del código de las funciones `sustitucion_adelante` y
`sustitucion_atras` que ya están en el archivo; lo nuevo es que ahora
trabajan con `self.T.data` y `self.b`.

### El cero en la diagonal: `try` / `except`

Si `T` tiene un cero en la diagonal, la división
`(b[i] - suma) / T[i][i]` levanta `ZeroDivisionError`. Ese mensaje
("float division by zero") no le dice al usuario qué hizo mal. Atrapen
esa excepción **en el lugar donde se divide** y levanten en su lugar
una `MatrizSingularError` que diga en qué renglón está el cero,
encadenándola con `from` (unidad 05):

```python
try:
    x[i] = (self.b[i] - suma) / T[i][i]
except ZeroDivisionError as error:
    raise MatrizSingularError(f"T[{i}][{i}] = 0: ...") from error
```

(Como la división aparece en los dos métodos, quizá les convenga
ponerla en un método auxiliar que usen ambos.)

¿Por qué no revisar `T.is_singular()` al principio y ya? También se
vale, y en `__init__` sería razonable. Pero aquí queremos practicar
`try`/`except`: en muchos problemas no se puede saber de antemano si
algo va a fallar, y lo que se hace es intentar y reaccionar.

### Contar operaciones: `self.flops`

Cada vez que se resuelve el sistema, `self.flops` debe quedar con el
número de operaciones de punto flotante que se hicieron (sumas,
restas, multiplicaciones y divisiones), contadas como en la sección
"Conteo de operaciones" de `notas.md`. Debe reiniciarse en cada
llamada: resolver dos veces el mismo sistema no duplica el conteo. Lo
van a usar en el Ejercicio 3; si lo hacen bien, debe dar exactamente
$n^2$.

### Docstrings

Documenten la clase (atributos y métodos, como el docstring de
`Matrix`) y cada método (`Parameters`, `Returns`, `Raises`). En
`Raises` pongan **todas** las excepciones que el método puede
levantar, incluyendo las que vienen de otros métodos que llama
(`resolver` puede levantar `MatrizSingularError`, aunque no tenga un
`raise` escrito directamente).

**Para probar a mano:** con la `L` del Ejercicio 1 y
`b = [2, 7, 26]`, la solución es `x = [1, 2, 2]`. Prueben también qué
pasa con `SistemaTriangular(Matrix([[1, 0], [0, 1]]), [1, 1])` y con
`SistemaTriangular(LowerTriangularMatrix([[1, 0], [1, 0]]), [1, 1]).resolver()`:
los mensajes de error deben bastar para entender qué salió mal, sin
leer el código.

## Ejercicio 3 — Complejidad: contar y medir

Creen un script `unidades/09_algebra_lineal/complejidad.py`. Como en
la Práctica 3, guarden los datos en la carpeta `datos/` junto al
script, usando `Path(__file__).resolve().parent / "datos"`.

### 3a. Contar a mano

En `notas.md` contamos los flops del producto matriz-vector
($`2n^2 - n`$) y de la sustitución hacia adelante ($`n^2`$), y dejamos
pendientes otros dos. En un comentario al inicio del script (o en el
docstring del módulo), cuenten, con el mismo razonamiento que en
`notas.md`, los flops de:

- el producto punto de dos vectores de $n$ componentes;
- el producto de dos matrices $n \times n$.

Den el número **exacto**, con todo y términos de grado menor, y luego
su orden $O(\cdot)$. Usen la misma convención que en las notas: sumar
$k$ números cuesta $k - 1$ sumas. (Pista: el producto de matrices es
$n^2$ productos punto.)

### 3b. Contar con código

Escriban versiones "contando" de tres operaciones, como
`sustitucion_adelante_contando` en
[`triangulares.py`](../../unidades/09_algebra_lineal/triangulares.py):
cada una regresa `(resultado, flops)`.

- `producto_punto_contando(x, y)`
- `mat_vec_contando(A, x)`
- `mat_mat_contando(A, B)`

(Pista: si escriben primero `producto_punto_contando`, las otras dos
se pueden escribir llamándola, y sumando los flops que regresa.) Para
la sustitución no hace falta una función nueva: usen
`SistemaTriangular` y lean su atributo `flops`.

Para `n` en `[1, 2, 4, 8, 16, 32, 64]` (o más valores), calculen los
cuatro conteos e impriman una tabla que los compare contra sus
fórmulas del 3a y las de `notas.md`. Escriban los resultados en
`datos/conteo_flops.dat`, con una línea de encabezado que empiece con
`#` y luego **cinco columnas separadas por espacios**:

```
n  producto_punto  mat_vec  mat_mat  sustitucion
```

Los conteos deben coincidir **exactamente** con las fórmulas. Si no
coinciden, o el código cuenta mal, o la fórmula está mal: averigüen
cuál.

### 3c. Medir el tiempo

Contar operaciones es teoría; ahora veamos si la computadora le hace
caso. Si un método cuesta $C n^p$ operaciones y cada operación tarda
más o menos lo mismo, el tiempo es $t(n) \approx C' n^p$, y al
**duplicar** $n$:

$$
\frac{t(2n)}{t(n)} \approx 2^p
\qquad\Longrightarrow\qquad
p \approx \log_2 \frac{t(2n)}{t(n)}.
$$

Midan, con `time.perf_counter()` como en `triangulares.py`, el tiempo
de:

- `SistemaTriangular(...).resolver()` (esperamos $`p \approx 2`$), y
- el producto `A * B` de dos `Matrix` $n \times n$ (esperamos
  $`p \approx 3`$),

para al menos 4 valores de `n`, duplicando cada vez (por ejemplo
`[25, 50, 100, 200]`; no se pasen de 200 para el producto de
matrices: con listas de Python ya tarda alrededor de un segundo).
Repitan cada medición unas 3 veces y quédense con la más rápida.

Escriban `datos/tiempos.dat`, con encabezado `#` y **tres columnas**:

```
n  t_sustitucion  t_mat_mat
```

e impriman, **calculado con código**, el exponente $p$ estimado entre
cada par de valores consecutivos de `n`.

### 3d. Preguntas

Respondan en comentarios al final del script:

1. ¿Qué valores de $p$ obtuvieron? ¿Se acercan más a 2 y a 3 para `n`
   chica o para `n` grande? ¿Por qué creen que para `n` chica el
   tiempo no sigue tan bien a la fórmula? (Piensen en qué más hace
   Python además de las operaciones de punto flotante: crear listas,
   llamar funciones, recorrer ciclos...)
2. La sustitución y el producto matriz-vector son ambos $O(n^2)$, pero
   con prefactores distintos ($`n^2`$ contra $`2n^2 - n`$). ¿Cuál esperan
   que sea más rápido para la misma `n`, y por cuánto, aproximadamente?
3. Si resolver un sistema triangular con $n = 1000$ tarda $t$, ¿cuánto
   tardaría con $n = 10\,000$? ¿Y un producto de matrices, si con
   $n = 1000$ tarda $t'$? En unas líneas: ¿por qué en álgebra lineal
   importa tanto si un método es $O(n^2)$ o $`O(n^3)`$?

### 3e. (Opcional) Graficar

Una vez generado `datos/tiempos.dat`, corran, desde la carpeta de la
unidad:

```bash
gnuplot graficar_complejidad.gp
```

Grafica los tiempos contra `n` en escala log-log, junto con rectas
proporcionales a $n^2$ y a $n^3$. En log-log, $t = C n^p$ es una recta
de pendiente $p$: comparen las pendientes de sus datos con las de las
rectas de referencia.

## Cómo revisar su trabajo

Hay un script de pruebas en
[`pruebas_practica_04.py`](pruebas_practica_04.py). Córranlo desde la
raíz del repositorio, con el entorno virtual activado:

```bash
python3 practicas/practica4/pruebas_practica_04.py
```

Revisa:

- **Ejercicio 1:** que las dos clases hereden de `Matrix`; que
  rechacen matrices no cuadradas, vacías, no rectangulares o con
  elementos del lado equivocado; y `diagonal()`, `determinant()`,
  `is_singular()` y `transpose()`.
- **Ejercicio 2:** que `SistemaTriangular` resuelva sistemas con
  solución conocida, con `resolver()` y con cada método; que levante
  `TypeError`, `ValueError` y `MatrizSingularError` en los casos
  descritos arriba; que no modifique la lista `b`; y que `flops`
  valga exactamente $n^2$.
- **Ejercicio 3:** corre `complejidad.py` completo y revisa que
  `datos/conteo_flops.dat` coincida exactamente con las fórmulas, y el
  formato de `datos/tiempos.dat` (que `n` se duplique y que el tiempo
  crezca).

Mientras una clase o el script no existan, sus pruebas aparecen como
`skipped`, con un mensaje que dice qué falta. Como en las prácticas
anteriores, los mensajes no dicen el valor esperado ni el obtenido.

No revisa los docstrings ni las respuestas del 3a y el 3d: esas partes
se revisan a mano.
