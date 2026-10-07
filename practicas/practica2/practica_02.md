# Práctica 2 — Programación orientada a objetos: `VectorND` y `Matrix`

Cubre lo visto en la unidad
[07 (programación orientada a objetos)](../../unidades/07_programacion_orientada_a_objetos/).
Es el primer paso hacia la unidad de **álgebra lineal numérica**.

## Parte 1 — `VectorND` (entregada en clase)

Esta parte ya se hizo en clase y se entregó por correo; queda
documentada aquí para que el registro de la práctica esté completo.

La entrega consiste en subir la versión realizada en clase como parte de esta practica. 
Recuerden que la versión que suban debe coincidir con la que me mandaron por correo.

## Parte 2 — `Matrix`

La clase `Matrix` que realicen aquí es la que van a seguir usando en
la unidad de álgebra lineal numérica.

### Contexto

En [`fiscomp/matrices.py`](../../fiscomp/matrices.py) ya está la clase
`Matrix`, con:

- El constructor `Matrix(data)`, que guarda una lista de listas
  rectangular en `self.data` y sus dimensiones en `self.rows`,
  `self.cols`, validando que no esté vacía y que todos los renglones
  tengan la misma longitud.
- `__str__`, para que `print(matriz)` se vea bien.
- `shape()`, `get_row(i)`, `get_col(j)`: para consultar la forma y
  acceder a renglones/columnas.
- `transpose()` y `copy()`.

Lo que falta — marcado con `TODO` en el archivo — son los tres
operadores aritméticos: `__add__`, `__sub__` y `__mul__`.

### Ejercicio 1 — Suma y resta (`__add__`, `__sub__`)

Completen `__add__` para que `matriz1 + matriz2` sume elemento a
elemento, y solo cuando ambas tengan la misma forma (`shape()`); si no
coinciden, levanten `ValueError`.

Para `__sub__`, no hace falta repetir el recorrido elemento por
elemento: `A - B` es lo mismo que `A + (B * -1)`, así que, una vez que
tengan `__add__` y `__mul__`, `__sub__` se puede escribir en una línea
reutilizándolos (revisen también ahí que las dimensiones coincidan).

### Ejercicio 2 — Multiplicación (`__mul__`)

`__mul__` tiene que distinguir dos casos, usando `isinstance`:

- Si `other` es un número (`int` o `float`): multiplicación por
  escalar, cada entrada de la matriz se multiplica por `other`.
- Si `other` es otra `Matrix`: multiplicación matricial. La entrada
  `(i, j)` del resultado es
  $\sum_k \text{self.data}[i][k] \cdot \text{other.data}[k][j]$.
  Solo está definida si el número de columnas de `self` coincide con
  el número de renglones de `other` (`self.cols == other.rows`) — si
  no, levanten `ValueError`.
- Si `other` no es ninguno de los dos casos anteriores: levanten
  `TypeError`.

Ya en el archivo está la línea `__rmul__ = __mul__`, así que en cuanto
completen `__mul__`, `escalar * matriz` va a funcionar igual que
`matriz * escalar` sin que tengan que escribir nada más (mismo truco
que usamos con `VectorND` en la unidad 07).

**Pista para probar la multiplicación matricial a mano:** con
$`A = \begin{pmatrix}1&2\\3&4\end{pmatrix}`$ y
$`B = \begin{pmatrix}5&6\\7&8\end{pmatrix}`$,
$`A \times B = \begin{pmatrix}19&22\\43&50\end{pmatrix}`$.

### Manejo de errores — resumen

Al terminar, `Matrix` debe levantar una excepción clara en estos
casos:

- Construir con renglones de longitudes distintas, o sin renglones
  (`ValueError`, ya implementado en `__init__`).
- Sumar o restar matrices de dimensiones distintas (`ValueError`).
- Multiplicar dos matrices cuyas dimensiones internas no coinciden
  (`ValueError`).
- Multiplicar una matriz por algo que no es ni número ni `Matrix`
  (`TypeError`).

### Cómo revisar su trabajo

Hay un script de pruebas en
[`pruebas_practica_02.py`](pruebas_practica_02.py) que revisa, sobre
matrices pequeñas de prueba: `shape()`/`get_row()`/`get_col()`, que
`print()` muestre las entradas, `+`, `-`, `*` por escalar (de los dos
lados), multiplicación de matrices (cuadradas y no cuadradas),
`transpose()`, y los cuatro casos de error de la sección anterior.

Como en la Práctica 1, las pruebas usan `subTest` (si una falla, las
demás se siguen corriendo) y no les dicen el valor esperado ni el
obtenido, solo qué caso falló. Córranlo desde la raíz del repositorio:

```bash
python3 practicas/practica2/pruebas_practica_02.py
```

Mientras `__add__`, `__sub__` o `__mul__` sigan con su
`NotImplementedError` del `TODO`, las pruebas que dependen de ellos
van a aparecer como `ERROR` (no `FAILED`) — es la señal de que a esa
parte todavía le falta código, no de que esté mal.
