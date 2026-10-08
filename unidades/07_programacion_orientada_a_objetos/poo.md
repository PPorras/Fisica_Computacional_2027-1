# Programación orientada a objetos (POO)

Notas de teoría sobre POO en general. Para la aritmética de intervalos
como aplicación concreta, ver [`notas.md`](notas.md); el código vive en
[`aritmetica_intervalos.py`](aritmetica_intervalos.py). Los ejemplos
de aquí usan una clase `Punto` más simple, solo para ilustrar cada
concepto por separado antes de verlos combinados en `Intervalo`; la
versión completa y ejecutable de `Punto` está en
[`punto.py`](punto.py).

## ¿Qué es?

Hasta ahora, el código del curso ha sido principalmente
**procedural**: datos por un lado (números, listas, diccionarios,
tuplas) y funciones sueltas por otro que reciben esos datos y los
transforman. La **programación orientada a objetos** (*object-oriented
programming*, OOP) propone empaquetar juntos los datos y las funciones
que operan sobre ellos en una sola unidad: el **objeto**.

No es un reemplazo de lo anterior — dentro de cada método de una clase
se sigue escribiendo código procedural normal (`if`, `for`, funciones
auxiliares, etc.) — sino una forma de organizar el código cuando hay
un concepto con **estado propio** (datos que le pertenecen) y
**comportamiento propio** (operaciones que solo tienen sentido para
ese estado). `Intervalo` es un buen ejemplo: sus datos son `lo` y
`hi`, y su comportamiento (`+`, `-`, `*`, `/`, "¿contiene a x?") solo
tiene sentido en términos de esos dos datos.

## Clases y objetos

Una **clase** (*class*) es una plantilla: describe qué atributos y
métodos va a tener cada ejemplar, pero no es, por sí misma, ningún
ejemplar concreto. Un **objeto** (*object*) o **instancia**
(*instance*) es un ejemplar concreto construido a partir de una clase,
con sus propios valores para esos atributos.

```python
class Punto:
    pass

p1 = Punto()  # p1 es una instancia de Punto
p2 = Punto()  # p2 es otra instancia, distinta de p1

print(type(p1))       # <class '__main__.Punto'>
print(p1 is p2)        # False -- son objetos distintos
print(isinstance(p1, Punto))  # True
```

La relación clase/objeto es la misma que la de tipo/valor que ya
conocen: `int` es un tipo y `3` es un valor de ese tipo; `Punto` es
una clase y `p1` es un objeto de esa clase. De hecho, en Python,
`int`, `list`, `dict`, etc. son clases igual que las que uno define, y
`3`, `[1, 2]`, `{"a": 1}` son instancias de ellas.

## Atributos e `__init__`

Un objeto vacío como `Punto()` no es muy útil todavía. El método
especial `__init__` (el **constructor**, aunque técnicamente
inicializa, no construye) se ejecuta automáticamente al crear un
objeto, y es donde normalmente se asignan sus **atributos de
instancia**:

```python
class Punto:
    def __init__(self, x, y):
        self.x = x
        self.y = y

p = Punto(3, 4)
print(p.x, p.y)  # 3 4
```

`self` es el primer parámetro de todo método de instancia, y es la
forma en que el método recibe una referencia al objeto sobre el que
fue llamado (`p.__init__(p, 3, 4)` es, conceptualmente, lo que pasa
por debajo cuando se escribe `Punto(3, 4)`). No es una palabra
reservada de Python — es solo una convención universal, casi nadie usa
otro nombre — pero **sí** hay que escribirlo explícitamente como
primer parámetro de cada método.

`self.x = x` crea un atributo `x` en *esa* instancia particular.
Instancias distintas tienen copias independientes:

```python
p1 = Punto(3, 4)
p2 = Punto(0, 0)
p1.x = 100
print(p1.x, p2.x)  # 100 0 -- no se afectan entre sí
```

## Métodos de instancia

Un **método** es una función definida dentro de una clase, que recibe
`self` como primer parámetro y típicamente opera sobre los atributos
de esa instancia:

```python
class Punto:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def distancia_al_origen(self):
        return (self.x**2 + self.y**2) ** 0.5

p = Punto(3, 4)
print(p.distancia_al_origen())  # 5.0
```

`p.distancia_al_origen()` es azúcar sintáctica por
`Punto.distancia_al_origen(p)`: Python pasa `p` automáticamente como
`self`. Esto es justo lo que ya hacían, sin saberlo, con métodos de
tipos incorporados como `lista.append(x)` (equivalente a
`list.append(lista, x)`) o `"abc".upper()`.

### `@property`: métodos que se usan como atributos

A veces conviene que algo calculado se vea, desde afuera, como un
atributo normal (sin paréntesis) en vez de como un método. El
decorador `@property` hace justo eso — ya lo usa `Intervalo.ancho` y
`Intervalo.centro`:

```python
class Punto:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    @property
    def distancia_al_origen(self):
        return (self.x**2 + self.y**2) ** 0.5

p = Punto(3, 4)
print(p.distancia_al_origen)  # 5.0 -- sin paréntesis
```

Es puramente una cuestión de qué tan natural se ve la sintaxis del
lado de quien usa la clase (`intervalo.ancho` se lee mejor que
`intervalo.ancho()`), no cambia nada sobre cómo se escribe el cuerpo
del método.

## Métodos especiales (*dunder methods*)

Los métodos cuyo nombre empieza y termina con doble guion bajo
(*dunder*, de *d*ouble *under*score) son "ganchos" que Python llama
automáticamente en situaciones específicas, en vez de tener que
llamarlos por su nombre. Ya vimos `__init__` (se llama al construir el
objeto); estos son los que más se usan:

| Método | Se llama cuando... |
|---|---|
| `__repr__(self)` | se pide una representación del objeto como texto (`repr(x)`, la consola interactiva, o dentro de un `!r` en un f-string) |
| `__str__(self)` | `str(x)` o `print(x)` (si no está definido, Python usa `__repr__` en su lugar) |
| `__eq__(self, otro)` | se evalúa `x == otro` |
| `__len__(self)` | se evalúa `len(x)` |
| `__contains__(self, item)` | se evalúa `item in x` |
| `__add__(self, otro)` | se evalúa `x + otro` |

```python
class Punto:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __repr__(self):
        return f"Punto({self.x!r}, {self.y!r})"

    def __eq__(self, otro):
        return self.x == otro.x and self.y == otro.y

p = Punto(3, 4)
print(p)               # Punto(3, 4) -- usa __repr__
print(p == Punto(3, 4))  # True -- usa __eq__
```

Sin `__repr__`, `print(p)` mostraría algo inútil como
`<__main__.Punto object at 0x7f...>` (la dirección de memoria del
objeto); sin `__eq__`, `p == Punto(3, 4)` daría `False`, porque el
`==` por default compara identidad (¿son el *mismo* objeto?), no
igualdad de contenido.

### Sobrecarga de operadores (*operator overloading*)

`__add__`, `__sub__`, `__mul__`, `__truediv__`, `__neg__`, etc. son la
razón por la que `Intervalo(1, 2) + Intervalo(3, 5)` funciona: `x + y`
en Python **siempre** se traduce, por debajo, a un llamado a método.
Concretamente, para `x + y`:

1. Python intenta `x.__add__(y)`.
2. Si ese método no existe, o regresa la constante especial
   `NotImplemented` (por ejemplo porque `y` es de un tipo que `x` no
   sabe cómo sumar), Python intenta `y.__radd__(x)` (el método
   "reflejado", de *reflected*).
3. Si ninguno de los dos funciona, se lanza `TypeError`.

Siguiendo con `Punto`, basta con definir `__add__` para que el `+`
funcione entre dos puntos (sumándolos como vectores, coordenada a
coordenada):

```python
class Punto:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __repr__(self):
        return f"Punto({self.x!r}, {self.y!r})"

    def __add__(self, otro):
        return Punto(self.x + otro.x, self.y + otro.y)

p1 = Punto(3, 4)
p2 = Punto(1, 1)
print(p1 + p2)  # Punto(4, 5) -- Python traduce esto a p1.__add__(p2)
```

Con solo esto, `p1 + p2` funciona, pero `p1 + 5` no: adentro de
`__add__`, `otro.x` fallaría con `AttributeError` porque un `int` no
tiene atributo `x`. Si quisiéramos que `Punto` también aceptara sumar
un escalar (interpretándolo como sumarlo a ambas coordenadas, digamos)
habría que revisar el tipo de `otro` dentro de `__add__` — exactamente
el papel que cumple `_coaccionar` en `Intervalo` (ver
`aritmetica_intervalos.py`), que convierte cualquier `int`/`float`
suelto en un `Intervalo` de ancho 0 antes de operar.

`__sub__` sigue exactamente el mismo patrón que `__add__` (resta
coordenada a coordenada), y `__mul__` aprovecha que, para un punto
visto como vector, "multiplicar" sí tiene un significado natural con
un escalar suelto (no con otro `Punto`) — escalar cada coordenada:

```python
class Punto:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __repr__(self):
        return f"Punto({self.x!r}, {self.y!r})"

    def __add__(self, otro):
        return Punto(self.x + otro.x, self.y + otro.y)

    def __sub__(self, otro):
        return Punto(self.x - otro.x, self.y - otro.y)

    def __mul__(self, escalar):
        return Punto(self.x * escalar, self.y * escalar)

    def __truediv__(self, escalar):
        return Punto(self.x / escalar, self.y / escalar)

    def __neg__(self):
        return Punto(-self.x, -self.y)

p1 = Punto(3, 4)
p2 = Punto(1, 1)
print(p1 - p2)  # Punto(2, 3) -- usa __sub__
print(p1 * 2)   # Punto(6, 8) -- usa __mul__, aquí sí con un int suelto
print(p1 / 2)   # Punto(1.5, 2.0) -- usa __truediv__, mismo patrón que __mul__
print(-p1)      # Punto(-3, -4) -- usa __neg__
```

Nótese la asimetría entre los cuatro: `__add__` y `__sub__` esperan
otro `Punto` del lado derecho (`otro.x`, `otro.y`), `__mul__` y
`__truediv__` esperan un número suelto (`escalar`) — no hay una única
forma "correcta" de definir estos métodos, cada uno recibe el tipo que
tenga sentido para la operación que están modelando. `__neg__` es
distinto de todos los anteriores porque es un **operador unario**: no
recibe ningún segundo operando (`-p1` solo necesita `self`), a
diferencia de `+`, `-`, `*`, `/`, que son binarios y siempre reciben
`otro`.

Y, a diferencia de la suma, `p1 - p2` y `p2 - p1` dan resultados
distintos (la resta no es conmutativa): eso es justo lo que hace que,
más abajo, `Intervalo` necesite un `__rsub__` separado (no un simple
alias a `__sub__`) para que `escalar - Intervalo(...)` reste en el
orden correcto — el mismo razonamiento aplica a `__truediv__` /
`__rtruediv__`, que tampoco es conmutativa.

Esto explica por qué `Intervalo` define tanto `__add__` como
`__radd__` (aliasado al mismo método, porque la suma es conmutativa):
`Intervalo(1,2) + 3` usa `__add__` normalmente, pero `3 + Intervalo(1,2)`
necesita `__radd__`, porque `int.__add__(3, Intervalo(1,2))` no sabe
qué hacer con un `Intervalo` y regresa `NotImplemented`. Para la resta
y la división, que **no** son conmutativas, hace falta escribir
`__rsub__`/`__rtruediv__` por separado (no como simple alias), como en
`aritmetica_intervalos.py`.

## Encapsulación

**Encapsulación** (*encapsulation*) es la idea de que un objeto expone
una interfaz (sus métodos) y oculta los detalles de cómo logra su
comportamiento (sus atributos internos), de forma que quien usa la
clase no necesita saber cómo está implementada por dentro, solo qué
métodos puede llamar y qué hacen.

Python no tiene atributos verdaderamente privados (a diferencia de
otros lenguajes de POO). Las convenciones son:

- `_atributo` (un guion bajo): "privado por convención" — indica "no
  lo uses desde afuera de la clase", pero Python no lo impide.
- `__atributo` (doble guion bajo): activa *name mangling* (Python lo
  renombra internamente a `_NombreDeClase__atributo`), lo que hace más
  difícil (no imposible) acceder desde afuera por accidente.

En este curso, con clases pequeñas como `Intervalo`, no hace falta
tanta ceremonia: sus atributos `lo`/`hi` son públicos a propósito
(leerlos desde afuera es razonable), y la "protección" real está en
que el `__init__` valida `lo <= hi` al crear el objeto.

## Herencia

El código completo de esta sección, listo para correr, está en
[`herencia.py`](herencia.py).

### La idea: "es un"

Muchas veces una clase nueva es un caso particular de una que ya
existe:

- una partícula puntual **es un** punto que además tiene masa;
- la matriz identidad **es una** matriz (con unos en la diagonal y
  ceros fuera);
- `TemperaturaImposibleError`, de la
  [unidad 05](../05_manejo_de_excepciones/), **es una** excepción.

Cuando la frase "B es un A" tiene sentido, conviene la **herencia**
(*inheritance*): la clase nueva, llamada **clase hija** (*subclass* o
*child class*), se define a partir de la **clase madre** (*superclass*,
*parent class* o *base class*) y recibe automáticamente todos sus
atributos y métodos. Solo hay que escribir lo que es nuevo o lo que
cambia. La sintaxis es poner la madre entre paréntesis:

```python
class Hija(Madre):
    ...
```

### Una clase hija vacía ya funciona

Para ver que la herencia de verdad trae todo, empecemos con una clase
hija que no agrega nada:

```python
class PuntoConMasa(Punto):
    pass

p = PuntoConMasa(3, 4)
print(p)                      # Punto(3, 4) -- __repr__ heredado
print(p.distancia_al_origen)  # 5.0 -- heredado
print(p + Punto(1, 1))        # Punto(4, 5) -- __add__ heredado
```

`PuntoConMasa(3, 4)` funciona aunque `PuntoConMasa` no define
`__init__`: Python usa el de `Punto`.

### Cómo encuentra Python un método

Cuando se escribe `p.distancia_al_origen`, Python busca ese nombre:

1. en el objeto mismo (sus atributos de instancia, como `x` y `y`);
2. si no está, en la clase del objeto, `PuntoConMasa`;
3. si no está, en la clase madre, `Punto`;
4. y así hacia arriba, hasta llegar a `object`.

Usa **el primero que encuentra**. Esa lista ordenada de clases se
llama el *method resolution order* y se puede ver con `__mro__`:

```python
print(PuntoConMasa.__mro__)
# (PuntoConMasa, Punto, object)  -- simplificado; Python imprime <class '...'>
```

`object` es la madre de todas las clases: escribir `class Punto:` es
lo mismo que `class Punto(object):`. De ahí salen los métodos "por
default" de la sección de métodos especiales: el `__repr__` que
imprime `<__main__.Punto object at 0x7f...>` y el `__eq__` que compara
identidad son los de `object`, y `Punto` los reemplaza con los suyos.

Como un `PuntoConMasa` es un `Punto`, `isinstance` lo reconoce como
las dos cosas. Al revés no: un `Punto` común no es un `PuntoConMasa`.
`issubclass` hace la misma pregunta, pero entre clases:

```python
p = PuntoConMasa(3, 4)
print(isinstance(p, PuntoConMasa))           # True
print(isinstance(p, Punto))                  # True
print(isinstance(Punto(0, 0), PuntoConMasa)) # False
print(issubclass(PuntoConMasa, Punto))       # True
```

### Agregar atributos: `__init__` y `super()`

Para que `PuntoConMasa` tenga masa, necesita su propio `__init__`.
Pero al definirlo, **reemplaza** al de `Punto`: Python ya no llama al
`__init__` de la madre por su cuenta. Para no repetir el código que
guarda `x` y `y`, se llama explícitamente con `super()`:

```python
class PuntoConMasa(Punto):
    def __init__(self, x, y, masa):
        super().__init__(x, y)  # Punto.__init__ guarda self.x y self.y
        if masa <= 0:
            raise ValueError("La masa debe ser positiva.")
        self.masa = masa

p = PuntoConMasa(3, 4, masa=2.0)
print(p.x, p.y, p.masa)  # 3 4 2.0
```

`super()` quiere decir "la clase madre, actuando sobre este mismo
`self`". Aquí, `super().__init__(x, y)` hace lo mismo que
`Punto.__init__(self, x, y)`, con la ventaja de que no hay que
escribir el nombre de la madre (ni pasar `self` a mano).

Un error común es olvidar esa línea. Python no se queja al crear el
objeto, pero `Punto.__init__` nunca se ejecutó, así que el objeto no
tiene `x` ni `y`:

```python
class PuntoSinSuper(Punto):
    def __init__(self, x, y, masa):
        self.masa = masa  # falta super().__init__(x, y)

roto = PuntoSinSuper(3, 4, masa=2.0)
print(roto.x)
# AttributeError: 'PuntoSinSuper' object has no attribute 'x'
```

### Agregar métodos

Un método nuevo se escribe igual que en cualquier clase. Dentro de
él, `self` tiene todo lo heredado, así que se puede usar sin más:

```python
class PuntoConMasa(Punto):
    # ... __init__ como arriba ...

    def momento_de_inercia(self):
        """Momento de inercia respecto al origen, m r^2."""
        return self.masa * self.distancia_al_origen**2

p = PuntoConMasa(3, 4, masa=2.0)
print(p.momento_de_inercia())  # 50.0, es decir 2 * 5^2
```

`distancia_al_origen` no está escrito en `PuntoConMasa`: Python lo
encuentra en `Punto`, siguiendo la búsqueda de arriba.

### Sobreescribir métodos (*overriding*)

Con lo que llevamos, `print(p)` todavía imprime `Punto(3, 4)`: usa el
`__repr__` de `Punto`, que no sabe nada de la masa. Si la hija define
un método **con el mismo nombre** que uno de la madre, la búsqueda
encuentra primero el de la hija, y ese es el que se usa. A eso se le
llama **sobreescribir** el método:

```python
class PuntoConMasa(Punto):
    # ... __init__ y momento_de_inercia como arriba ...

    def __repr__(self):
        return f"PuntoConMasa({self.x!r}, {self.y!r}, masa={self.masa!r})"

print(PuntoConMasa(3, 4, masa=2.0))  # PuntoConMasa(3, 4, masa=2.0)
```

Al sobreescribir no hay que reescribir todo desde cero: se puede usar
`super().metodo(...)` para que la madre haga su parte y luego agregar
lo nuevo, como hicimos con `__init__`.

### Cuidado con lo que se hereda tal cual

La hija hereda **todos** los métodos de la madre, y cada uno hay que
revisarlo: ¿sigue teniendo sentido para la hija? En `PuntoConMasa` hay
dos sorpresas:

```python
p = PuntoConMasa(3, 4, masa=2.0)
q = PuntoConMasa(1, 1, masa=5.0)

print(p == PuntoConMasa(3, 4, masa=9.0))  # True -- ¡no compara la masa!
print(p + q)                              # Punto(4, 5) -- se perdió la masa
```

- `__eq__` es el de `Punto`, que solo compara `x` y `y`.
- `__add__` es el de `Punto`, que construye un `Punto(...)` nuevo. El
  resultado es un `Punto` común, sin masa.

¿Es eso un error? Depende de la física. No es obvio qué debería ser
"la suma de dos partículas", así que regresar un `Punto` es una opción
razonable. Si se quiere otra cosa, hay que sobreescribir el método en
la hija.

### Un ejemplo del curso: `Identity(Matrix)`

En [`fiscomp/matrices.py`](../../fiscomp/matrices.py), la clase
`Identity` es una clase hija de `Matrix`:

```python
class Identity(Matrix):
    def __init__(self, n: int):
        if n < 1:
            raise ValueError("La dimensión n debe ser al menos 1.")
        data = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
        super().__init__(data)

    def __mul__(self, other):
        if isinstance(other, Matrix):
            if self.cols != other.rows:
                raise ValueError("Dimensiones incompatibles para multiplicación.")
            return other.copy()
        return super().__mul__(other)
```

Ahí aparece todo lo de esta sección:

- **Un constructor distinto:** `Identity(3)` recibe la dimensión en
  lugar de una lista de listas. Arma esa lista y se la pasa a
  `Matrix.__init__` con `super()`, que valida y guarda `data`, `rows`
  y `cols`.
- **Métodos heredados:** `shape()`, `transpose()`, `+` y `-` no se
  escriben de nuevo.
- **Un método sobreescrito:** como $I A = A$, `I * A` regresa una
  copia de `A` sin hacer las $n^3$ multiplicaciones. Para un escalar,
  `I * 2`, le deja el trabajo a la madre con `super().__mul__(other)`.
- **El tipo del resultado:** igual que `p + q` arriba, `2 * I` o
  `I + A` regresan una `Matrix` común. Aquí eso es justo lo correcto,
  porque ya no son la identidad.

En la [práctica 4](../../practicas/practica4/practica_04.md) ustedes
escriben otras dos clases hijas de `Matrix`: las matrices
triangulares.

### Otro ejemplo que ya usaron: las excepciones

`class TemperaturaImposibleError(Exception)` es herencia, y las
excepciones de Python forman una jerarquía completa:

```python
print(ZeroDivisionError.__mro__)
# (ZeroDivisionError, ArithmeticError, Exception, BaseException, object)  -- simplificado
```

Eso explica cómo funciona `except`: **atrapa la clase que se nombra y
todas sus hijas**, igual que `isinstance`. `except ArithmeticError`
atrapa un `ZeroDivisionError` (o un `OverflowError`), pero no un
`ValueError`, que no está en esa rama. Y `except Exception` atrapa
casi todo, porque casi todas las excepciones son hijas de
`Exception`.

### Cuándo no usar herencia: "es un" contra "tiene un"

Heredar solo para reutilizar código de una clase que "se parece" es un
error común. `Intervalo` y `Punto` guardan dos números cada uno, pero
un intervalo **no es** un punto. Si escribiéramos
`class Intervalo(Punto)`, heredaría:

- `distancia_al_origen`, que no significa nada para un intervalo;
- un `__sub__` que calcula $`[a-c,\; b-d]`$, cuando la resta de
  intervalos es $`[a-d,\; b-c]`$.

Cuando un objeto **tiene** otro adentro, en vez de **ser** uno, lo
correcto es guardarlo como atributo. A eso se le llama **composición**
(*composition*). Por ejemplo, una partícula que se mueve *tiene* una
posición y una velocidad:

```python
class Particula:
    def __init__(self, posicion, velocidad, masa):
        self.posicion = posicion    # un Punto
        self.velocidad = velocidad  # otro Punto
        self.masa = masa
```

Aquí no hay herencia: `Particula` usa dos objetos `Punto`, pero no es
uno.

(Python también permite heredar de varias madres a la vez,
`class C(A, B):`, lo que se llama herencia múltiple. No la usamos en
el curso.)

### Para pensar

1. ¿Qué imprime `PuntoConMasa(3, 4, masa=2.0) * 2`? ¿De qué tipo es el
   resultado? Sobreescriban `__mul__` en `PuntoConMasa` para que el
   resultado conserve la masa.
2. Sobreescriban `__eq__` en `PuntoConMasa` para que también compare
   la masa, usando `super().__eq__(otro)` para la parte de la
   posición.
3. En `fiscomp/matrices.py`, `I * A` usa el atajo de `Identity`, pero
   `A * I` hace todas las multiplicaciones. ¿Por qué? (Pista: ¿el
   `__mul__` de cuál de los dos objetos llama Python?)

## Polimorfismo

**Polimorfismo** (*polymorphism*) quiere decir que el mismo código
funciona con objetos de clases distintas, siempre que tengan los
métodos que ese código usa.

Con herencia es inmediato: una función escrita para `Punto` funciona
con un `PuntoConMasa`, porque tiene todo lo que tiene un `Punto`. Y
cuando la hija sobreescribe un método, cada objeto usa **el de su
propia clase**:

```python
objetos = [Punto(1, 2), PuntoConMasa(3, 4, masa=2.0), Punto(0, 5)]
for objeto in objetos:
    print(objeto)
# Punto(1, 2)
# PuntoConMasa(3, 4, masa=2.0)
# Punto(0, 5)
```

La línea `print(objeto)` es la misma, pero el segundo objeto usa el
`__repr__` de `PuntoConMasa`.

Pero el polimorfismo no necesita herencia. Código que llama a un
método (p. ej. `f(x + y)`) funciona igual sin importar la clase
concreta de `x`, mientras esa clase implemente `__add__`. Esto es lo
que permite que `x**2 - 2*x` en `aritmetica_intervalos.py` funcione
idéntico ya sea que `x` sea un `float` o un `Intervalo`, aunque
`Intervalo` no hereda de `float`: la función `f` no sabe, ni le
importa, cuál de los dos es.

## Por qué usar una clase aquí (y no, digamos, una tupla `(lo, hi)`)

Con funciones sueltas, `sumar_intervalos((1,2), (3,5))` también
resolvería el problema numérico. La ventaja de la clase es la sintaxis
y las garantías:

- `X + Y` en vez de `sumar_intervalos(X, Y)` — las expresiones se leen
  como las fórmulas matemáticas que representan.
- `__init__` valida `lo <= hi` **una sola vez**, al crear el objeto —
  con tuplas sueltas, cada función tendría que validar sus argumentos
  por separado (o confiar en que ya vienen bien formados).
- Un `Intervalo(1, 2)` se puede pasar a cualquier función que espere
  un número (gracias a `__add__`, `__pow__`, etc.) sin que esa función
  tenga que enterarse de que está tratando con un intervalo — es el
  polimorfismo de la sección anterior, aplicado directamente al
  ejemplo de acotar el rango de `f(x) = x**2 - 2*x` sobre un intervalo
  en vez de un float.
