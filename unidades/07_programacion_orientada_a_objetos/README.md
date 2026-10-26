# Unidad 07 — Programación orientada a objetos

## Temas
- **Python (Capítulo 9 — Clases):**
  - Definición de una clase, atributos de instancia, `__init__` y
    `self`
  - Métodos de instancia, `@property`
  - Métodos especiales (*dunder methods*): `__repr__`, `__eq__`,
    sobrecarga de operadores aritméticos (`__add__`, `__sub__`,
    `__mul__`, `__truediv__`, `__neg__`), `__contains__`
  - Métodos "reflejados" (`__radd__`, `__rmul__`, ...) para que la
    clase también funcione del lado derecho de un operador con un
    `int`/`float`
  - Herencia: clases hijas, `super()`, sobreescritura de métodos,
    búsqueda de métodos (`__mro__`), `isinstance`/`issubclass`, la
    jerarquía de excepciones, y cuándo conviene la composición en su
    lugar; polimorfismo
- **Física computacional:** aritmética de intervalos como caso de uso
  de POO — representar un intervalo $[a,b]$ como un objeto con sus
  propios operadores, y usarlo para acotar rigurosamente el rango de
  una función (en vez de solo aproximarlo con muestreo).

## Definiciones

Un **intervalo** $X = [a, b]$, con $a \le b$, representa el conjunto
de reales $`\{x \in \mathbb{R} : a \le x \le b\}`$. La **aritmética de
intervalos** extiende las operaciones aritméticas a intervalos, de
forma que el resultado siempre contenga todos los resultados posibles
de operar cualquier par de reales tomados de los intervalos
originales.

Para $X = [a,b]$ y $Y = [c,d]$:

```math
X + Y = [a+c,\; b+d]
```

```math
X - Y = [a-d,\; b-c]
```

```math
X \times Y = [\min(ac,ad,bc,bd),\; \max(ac,ad,bc,bd)]
```

```math
X \div Y = X \times \left[\tfrac{1}{d}, \tfrac{1}{c}\right],
\qquad \text{si } 0 \notin Y
```

(la división no está definida si $0 \in Y$, porque $1/y$ no está
acotado cerca de $y=0$).

### Por qué importa para cómputo numérico

Si a cada operación de punto flotante se le aplica **redondeo
dirigido hacia afuera** (el extremo inferior se redondea hacia $-\infty$
y el superior hacia $+\infty$), el intervalo resultante contiene
garantizadamente el resultado exacto, a pesar del error de redondeo de
la unidad 06. Esto convierte a la aritmética de intervalos en una
herramienta para obtener cotas de error *rigurosas*, en vez de solo
estimarlas. En este curso no implementamos
el redondeo dirigido (Python no lo expone directamente), pero las
fórmulas de arriba ya bastan para ilustrar la idea central.

### El "problema de dependencia" (*dependency problem*)

Una propiedad contraintuitiva: si la misma variable aparece más de una
vez en una expresión, la aritmética de intervalos la trata como si
fueran variables *independientes* en cada aparición. Por ejemplo, con
$X = [1,2]$:

```math
X - X = [1,2] - [1,2] = [1-2,\; 2-1] = [-1, 1] \neq [0,0].
```

Por eso, la misma función matemática puede dar cotas más anchas de lo
necesario según cómo esté escrita la expresión (p. ej. $X^2$ da una
cota más angosta que $X \times X$ cuando $X$ contiene negativos, ya
que `__mul__` no sabe que ambos factores son "el mismo" intervalo).

## Contenido
- [`poo.md`](poo.md): teoría general de programación orientada a
  objetos (clases, objetos, `__init__`, métodos, `@property`, métodos
  especiales, sobrecarga de operadores, encapsulación, herencia y
  polimorfismo), con una clase `Punto` de ejemplo.
- [`punto.py`](punto.py): versión completa y ejecutable de la clase
  `Punto` que se va construyendo por partes en `poo.md`.
- [`herencia.py`](herencia.py): la clase hija `PuntoConMasa(Punto)`
  (con `super()`, un método nuevo y `__repr__` sobreescrito), los
  métodos heredados que sorprenden, `__mro__`, `isinstance`,
  polimorfismo y la jerarquía de excepciones.
- [`aritmetica_intervalos.py`](aritmetica_intervalos.py): clase
  `Intervalo` con los operadores aritméticos sobrecargados, y ejemplos
  que ilustran la contención garantizada, el problema de dependencia
  ($X - X \neq [0,0]$), y cómo evaluar una función sobre un intervalo
  para acotar su rango.
