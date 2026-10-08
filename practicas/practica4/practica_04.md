# Práctica 4 — Matrices triangulares, sistemas lineales con objetos y complejidad

Cubre lo visto en la unidad
[09 (álgebra lineal)](../../unidades/09_algebra_lineal/), secciones
"Matrices triangulares", "Eliminación gaussiana" y "Descomposición
LU" de [`notas.md`](../../unidades/09_algebra_lineal/notas.md), y retoma la
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
Luego medimos cuánto cuesta cada operación, y al final usamos todas
esas piezas para resolver sistemas lineales generales con objetos
(eliminación gaussiana y LU, con la inversa y el determinante).

No borren `sustitucion_adelante` ni `sustitucion_atras`: las usa
[`triangulares.py`](../../unidades/09_algebra_lineal/triangulares.py)
y las reutilizan `eliminacion_gaussiana` y `resolver_lu`, en el mismo
archivo. Lo que escriban en
esta práctica va **junto** a ellas.

## Ejercicio 1 — Clases hijas de `Matrix`: matrices triangulares

### Herencia y `super()`

La **herencia** se explica con detalle en la sección "Herencia" de
[`poo.md`](../../unidades/07_programacion_orientada_a_objetos/poo.md),
en la unidad 07; léanla antes de empezar. En resumen:
`class Hija(Madre):` define una clase que tiene automáticamente todos
los atributos y métodos de `Madre`, y que puede agregar métodos nuevos
o *sobreescribir* (redefinir) algunos. Una matriz triangular **es
una** matriz, así que no tiene sentido volver a escribir `shape()`,
`get_row()`, `+`, `*`, etc.: los heredamos de `Matrix`.

Lo que sí cambia es la validación al crearla. Para no repetir la que
ya hace `Matrix.__init__` (que no esté vacía y que sea rectangular),
desde el `__init__` de la clase hija se llama al de la clase madre con
`super()`:

```python
class LowerTriangularMatrix(Matrix):
    def __init__(self, data):
        super().__init__(data)   # corre Matrix.__init__: llena self.data, self.rows, self.cols
        # ... y aquí van las validaciones nuevas
```

En [`fiscomp/matrices.py`](../../fiscomp/matrices.py) ya hay una
clase hija de `Matrix` escrita así, `Identity`, que les puede servir
de modelo (está explicada en `poo.md`, en "Un ejemplo del curso").

Ojo con lo que se hereda tal cual (ver "Cuidado con lo que se hereda
tal cual" en `poo.md`): los métodos de `Matrix` construyen una
`Matrix(...)` nueva, así que `L + L` o `2 * L` regresan una `Matrix`
común, no una `LowerTriangularMatrix`, aunque el resultado siga
siendo triangular. En esta práctica eso está bien, y no hace falta
cambiarlo; el único método heredado que sí tienen que sobreescribir
es `transpose()` (abajo).

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
triangular superior, y viceversa. Pueden reutilizar el método de la madre
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

```math
\frac{t(2n)}{t(n)} \approx 2^p
\qquad\Longrightarrow\qquad
p \approx \log_2 \frac{t(2n)}{t(n)}.
```

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
   tardaría con $`n = 10\,000`$? ¿Y un producto de matrices, si con
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

## Ejercicio 4 — Sistemas lineales con objetos: eliminación gaussiana y LU

En [`fiscomp/algebra_lineal.py`](../../fiscomp/algebra_lineal.py) ya
están, como funciones sobre listas de listas, la eliminación
gaussiana, la descomposición LU, la inversa y el determinante
(secciones "Eliminación gaussiana" y "Descomposición LU" de
[`notas.md`](../../unidades/09_algebra_lineal/notas.md)). Ahora las
vamos a construir con objetos, **por piezas y en orden**, de forma
que cada pieza nueva use las anteriores en vez de repetir su código:

1. Sustitución hacia adelante y hacia atrás: ya la tienen, es
   `SistemaTriangular` (Ejercicio 2).
2. Eliminación gaussiana: la clase `EliminacionGaussiana` (4a).
3. Resolver un sistema por eliminación gaussiana: la clase
   `SistemaLineal`, que usa la eliminación y la sustitución (4b).
4. La descomposición LU: la clase `FactorizacionLU`, que usa la misma
   eliminación (4c).
5. La inversa y el determinante, con la LU (4d).

Todo va en `fiscomp/algebra_lineal.py`, después de
`SistemaTriangular`. Agreguen a los imports del inicio lo que
necesiten de `fiscomp.matrices` (`Matrix`, `Identity`, ...).

### Diseño: SOLID

**SOLID** es un conjunto de cinco principios para diseñar programas
con objetos (cada letra es uno). No son reglas del lenguaje: Python
no los revisa. Son guías para que el código sea fácil de entender, de
corregir y de extender. Así se ven en este ejercicio:

- **S, responsabilidad única** (*single responsibility*): cada clase
  hace **una** cosa. `SistemaTriangular` sustituye,
  `EliminacionGaussiana` elimina, `FactorizacionLU` factoriza y usa la
  factorización, y `SistemaLineal` representa un sistema y lo
  resuelve. Por eso en todo el archivo debe haber **un solo** ciclo de
  eliminación (en `EliminacionGaussiana`) y **un solo** par de ciclos
  de sustitución (en `SistemaTriangular`): si hay un error en la
  eliminación, se corrige en un lugar y queda corregido para Gauss y
  para LU.
- **O, abierto/cerrado** (*open/closed*): una clase debe estar
  abierta a extensiones y cerrada a modificaciones. Por ejemplo, el
  pivoteo se puede agregar con una clase hija de
  `EliminacionGaussiana` (ver el opcional, al final), sin tocar las
  demás clases.
- **L, sustitución de Liskov** (*Liskov substitution*): un objeto de
  una clase hija debe servir en cualquier lugar donde se espera uno de
  la madre. Una `LowerTriangularMatrix` es una `Matrix`, así que
  `SistemaLineal` debe aceptarla igual que a cualquier `Matrix`.
- **I, segregación de interfaces** (*interface segregation*): clases
  con pocos métodos, que tengan sentido para quien las usa.
  `SistemaTriangular` no tiene `inversa()`; quien solo quiere resolver
  un sistema no tiene que saber nada de inversas.
- **D, inversión de dependencias** (*dependency inversion*): una pieza
  depende de **lo que hace** otra (sus métodos), no de **cómo** lo
  hace. `FactorizacionLU` le pide a `SistemaTriangular` que resuelva,
  sin saber cómo está escrita la sustitución: si mañana la cambian,
  `FactorizacionLU` no se toca.

En la práctica, lo que más van a usar es la S: antes de escribir un
ciclo, pregúntense si ya existe una pieza que lo hace.

### 4a. La eliminación gaussiana: `EliminacionGaussiana`

Primero, una excepción para cuando la eliminación encuentra un pivote
cero. Sin pivoteo, la eliminación falla si un pivote $A_{jj}$ vale
cero, **aunque la matriz no sea singular** (sección "Un pivote cero"
de `notas.md`: el ejemplo $`\begin{pmatrix} 0 & -1 \\ 1 & 1 \end{pmatrix}`$
tiene determinante 1). Ese error no es "la matriz es singular" sino
"este método no sirve para esta matriz", así que merece su propia
excepción:

```python
class PivoteCeroError(ZeroDivisionError):
    """..."""
```

Hereda de `ZeroDivisionError`, que es lo que de verdad pasó: quien ya
atrape `ZeroDivisionError` la atrapa también. Como en el Ejercicio 2,
atrapen el `ZeroDivisionError` **donde se divide entre el pivote** y
levanten en su lugar una `PivoteCeroError`, encadenada con `from`, con
un mensaje que diga en qué renglón quedó el pivote cero y que sugiera
el pivoteo.

La clase `EliminacionGaussiana(A, b=None, coeficientes=False)` hace
**solo la fase de eliminación** (sección "Caso general" de
`notas.md`), al crearse, sobre **copias** de `A.data` y de `b`:

- `A` debe ser una `Matrix` cuadrada; `b` es opcional.
- Si se pasa `b`, se le aplican las mismas operaciones de renglón que
  a `A`.
- `coeficientes` es una bandera (*flag*): si es `True`, además se
  guardan los coeficientes que se usaron para eliminar, en la matriz
  triangular inferior que será la $L$ de la LU (sección "La
  descomposición de Doolittle": $L_{ij}$ es el coeficiente que anuló
  a $A_{ij}$, con unos en la diagonal). Si es `False`, no se construye.

Así la misma eliminación sirve para los dos usos: para resolver un
sistema se necesita la triangular superior y el `b` transformado; para
la LU se necesita la triangular superior y los coeficientes.

Al terminar, el objeto tiene los atributos:

| Atributo | Valor |
|---|---|
| `U` | la triangular superior que quedó, como `UpperTriangularMatrix` |
| `b_nuevo` | el `b` transformado (lista), o `None` si no se pasó `b` |
| `L` | la matriz de coeficientes como `LowerTriangularMatrix`, con unos en la diagonal, o `None` si `coeficientes=False` |
| `n` | el tamaño de `A` |

No debe modificar ni `A` ni `b`. Levanta:

- `TypeError` si `A` no es una `Matrix` (por ejemplo, si es una lista
  de listas).
- `ValueError` si `A` no es cuadrada, o si se pasó un `b` que no tiene
  $n$ componentes.
- `PivoteCeroError` si aparece un pivote cero.

**Un detalle de punto flotante: ceros que no son cero.** El elemento
$A_{ij}$ de abajo de la diagonal se "anula" con

```python
coeficiente = A[i][j] / A[j][j]
A[i][j] -= coeficiente * A[j][j]
```

En aritmética exacta eso da 0, pero en punto flotante (unidad 06) no
siempre: con $A_{jj} = 0.3$ y $A_{ij} = 0.7$ queda
$`-1.1 \times 10^{-16}`$. Y entonces `UpperTriangularMatrix` rechaza
la $U$, porque tiene un elemento distinto de cero abajo de la
diagonal. Como sabemos que ese elemento *debe* ser cero, la solución
es no calcularlo: asígnenle `0.0` directamente (o hagan que el ciclo
sobre las columnas empiece en `j + 1`). Las pruebas incluyen una
matriz con ese problema.

**Para probar a mano:** con el ejemplo $3\times 3$ de `notas.md`,

```python
A = Matrix([[2, 1, 1], [1, 1, -2], [5, 10, 5]])
e = EliminacionGaussiana(A, [8, -2, 10], coeficientes=True)
```

deben obtener `e.U.data == [[2, 1, 1], [0, 0.5, -2.5], [0, 0, 40]]`,
`e.b_nuevo == [8, -6, 80]` y
`e.L.data == [[1, 0, 0], [0.5, 1, 0], [2.5, 15, 1]]` (compárenlo con
las cuentas a mano de `notas.md`).

### 4b. Resolver por eliminación gaussiana: `SistemaLineal`

`SistemaLineal(A, b)` representa el sistema $A\mathbf{x} = \mathbf{b}$
con $A$ cuadrada cualquiera, como `SistemaTriangular` lo hace para $T$
triangular:

- **`__init__(self, A, b)`**: guarda `self.A`, `self.b` (una copia),
  `self.n`. Levanta `TypeError` si `A` no es una `Matrix`, y
  `ValueError` si no es cuadrada o si `b` no tiene $n$ componentes.
  (Las mismas revisiones que en `EliminacionGaussiana`: escríbanlas
  una vez, en una función auxiliar que usen las dos clases.)
- **`resolver_gauss(self)`**: resuelve el sistema por eliminación
  gaussiana, **sin escribir ningún ciclo**: crea una
  `EliminacionGaussiana(self.A, self.b)` y resuelve el sistema
  triangular que quedó, `U x = b_nuevo`, con un `SistemaTriangular`.
  Regresa `x`.
- **`residuo(self, x)`**: el vector $\mathbf{b} - A\mathbf{x}$, con la
  función `residuo` de arriba.

Revisen qué pasa con una matriz singular, por ejemplo
$`\begin{pmatrix} 1 & 2 \\ 2 & 4 \end{pmatrix}`$: la eliminación sí
termina (el pivote cero aparece hasta el último renglón, donde ya no se
divide entre él), y es la sustitución la que levanta la
`MatrizSingularError` del Ejercicio 2, sin que ustedes escriban ningún
`raise` nuevo. En `Raises` de los docstrings, documenten también las
excepciones que vienen de las piezas que usan.

### 4c. La descomposición LU: `FactorizacionLU`

¿Por qué otra clase, si `SistemaLineal` ya resuelve? Porque la
eliminación gaussiana resuelve **un** sistema, pero la LU es una
propiedad de la matriz $A$ sola: se calcula **una vez** y luego sirve
para resolver con tantos $\mathbf{b}$ como se quiera, cada uno con
solo dos sustituciones ($`2n^2`$ flops, contra $`2n^3/3`$ de
factorizar). La inversa aprovecha eso ($n$ sistemas con la misma
$A$), y en el tema de eigenvalores vamos a resolver docenas de
sistemas con la misma matriz (la iteración inversa): ahí un objeto que
guarda su $L$ y su $U$ va a ser la pieza central.

`FactorizacionLU(A)` debe tener:

- **`__init__(self, A)`**: calcula $A = LU$ **con
  `EliminacionGaussiana(A, coeficientes=True)`**, sin `b` y sin
  escribir otra vez la eliminación, y guarda `self.A`, `self.n`,
  `self.L` y `self.U`. Levanta las mismas excepciones que
  `EliminacionGaussiana` (que vienen de ella).
- **`resolver(self, b)`**: resuelve $A\mathbf{x} = \mathbf{b}$ en dos
  pasos, $L\mathbf{y} = \mathbf{b}$ y luego $U\mathbf{x} = \mathbf{y}$,
  cada uno con un `SistemaTriangular`. Regresa `x`.

Y agreguen a `SistemaLineal` el método **`resolver_lu(self)`**, que
resuelve con `FactorizacionLU(self.A).resolver(self.b)`.

### 4d. La inversa y el determinante

Agreguen a `FactorizacionLU`:

- **`determinante(self)`**: $\det A = \det L \det U = \det U$, porque
  $L$ tiene unos en la diagonal. Usen `determinant()` del
  Ejercicio 1.
- **`inversa(self)`**: la inversa como `Matrix`. Su columna $k$ es la
  solución de $A\mathbf{x}_k = \mathbf{e}_k$, con $`\mathbf{e}_k`$ la
  columna $k$ de la identidad (sección "La matriz inversa" de
  `notas.md`): saquen $`\mathbf{e}_k`$ de `Identity(n).get_col(k)` y
  resuelvan con `self.resolver`, sin volver a factorizar.

Con la matriz singular de arriba, `determinante()` da 0 e `inversa()`
levanta `MatrizSingularError`.

**Para probar a mano:** con el ejemplo $3\times 3$ y `b = [8, -2, 10]`,
la solución es `x = [4, -2, 2]` con `resolver_gauss()` y con
`resolver_lu()`, `FactorizacionLU(A).determinante()` es `40`, y
`A * FactorizacionLU(A).inversa()` es, salvo redondeo, la identidad.
Prueben también `FactorizacionLU(Matrix([[0, -1], [1, 1]]))` y
`SistemaLineal(Matrix([[1, 2], [2, 4]]), [1, 1]).resolver_gauss()`:
como en el Ejercicio 2, los mensajes de error deben bastar para
entender qué salió mal.

### 4e. (Opcional) Pivoteo, sin modificar lo anterior

No lo revisan las pruebas. Escriban una clase hija
`EliminacionConPivoteo(EliminacionGaussiana)` que haga pivoteo parcial
(sección "Pivoteo" de `notas.md`): antes de eliminar la columna $j$,
intercambia el renglón $j$ con el que tenga el $|A_{ij}|$ más grande
(en `A`, en `b` y, si `coeficientes=True`, en los coeficientes ya
guardados). Comprueben que resuelve el ejemplo
$`\begin{pmatrix} 0 & -1 \\ 1 & 1 \end{pmatrix}`$ que
`EliminacionGaussiana` no puede.

Para que la hija solo tenga que sobreescribir un método pequeño (por
ejemplo, `elegir_pivote(self, j)`), quizá les convenga separar ese
paso en su propio método dentro de `EliminacionGaussiana`. Eso es el
principio O: la clase madre queda abierta a extensiones sin
modificarla después. Piensen: ¿por qué `FactorizacionLU` **no**
puede usar esta clase tal cual, sin guardar también los intercambios?
(Pista: con los renglones intercambiados, $LU$ ya no es $A$.)

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
- **Ejercicio 4:** que `EliminacionGaussiana` deje la $U$ y el $b$
  transformado correctos, y la $L$ solo si `coeficientes=True`
  (también cuando la eliminación deja "ceros" de redondeo); que
  `SistemaLineal` resuelva con `resolver_gauss()` y `resolver_lu()`
  sin modificar `A` ni `b`; que `FactorizacionLU` tenga $LU = A$,
  resuelva varios sistemas con la misma factorización y calcule el
  determinante y la inversa; y que se levanten `TypeError`,
  `ValueError`, `PivoteCeroError` y `MatrizSingularError` en los casos
  descritos.

Mientras una clase o el script no existan, sus pruebas aparecen como
`skipped`, con un mensaje que dice qué falta. Como en las prácticas
anteriores, los mensajes no dicen el valor esperado ni el obtenido.

No revisa los docstrings, las respuestas del 3a y el 3d, ni el diseño
del Ejercicio 4 (que cada pieza use las anteriores en vez de repetir
sus ciclos): esas partes se revisan a mano.
