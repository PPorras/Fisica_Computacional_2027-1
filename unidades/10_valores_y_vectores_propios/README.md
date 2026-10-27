# Unidad 10 — Valores y vectores propios

## Temas
- **Física computacional:** el problema de eigenvalores (valores
  propios) y eigenvectores (vectores propios)
  $A\mathbf{v} = \lambda\mathbf{v}$; análisis de error: número de
  condición para eigenvalores simples y sensibilidad de los
  eigenvectores; el problema de eigenvalores numéricamente
  (diagonalización, por qué no usar el polinomio característico); el
  método de la potencia, el cociente de Rayleigh y el método de la
  potencia inversa con desplazamiento; la descomposición QR
  (Gram-Schmidt), las transformaciones de similitud, la iteración
  simultánea y el método QR para todos los eigenvalores; todos los
  eigenpares combinando QR con la potencia inversa; la descomposición
  en valores singulares (SVD): definición, prueba de existencia, norma
  2 y número de condición, y su cálculo como problema de eigenvalores.
- La mayor parte de este material se ha tomado de las secciones 4.4
  (hasta la 4.4.4) y 4.5 de [1], y de la parte de eigenvalores de la
  4.2.
- Se apoya en la [unidad 09](../09_algebra_lineal/) (sistemas de
  ecuaciones lineales): normas, número de condición, descomposición LU
  y el criterio de paro de Jacobi.

## Notación

Los índices van de $0$ a $n-1$, como en Python. Las
matrices se escriben con mayúsculas ($A$, $L$, $U$) y los vectores en
negritas minúsculas ($\mathbf{x}$, $\mathbf{b}$). Cuando hace falta,
los dos índices de un elemento se separan con una coma: $A_{1,n-1}$.
Un superíndice entre paréntesis es el número de iteración
($`A^{(k)}`$ es la matriz de la iteración $k$) y uno sin paréntesis es
una potencia ($`A^k`$ es $A$ multiplicada por sí misma $k$ veces).
Usamos "eigenvalor" y "valor propio" como sinónimos (igual
"eigenvector" y "vector propio"); "eigenpar" es un eigenvalor junto
con su eigenvector.

En el código, como en la [unidad 09](../09_algebra_lineal/), una
matriz es un objeto de la clase `Matrix` de la práctica 2
([`fiscomp/matrices.py`](../../fiscomp/matrices.py)), con `A[i][j]`
para el elemento del renglón `i`, columna `j`, y un vector es una
lista simple (`x[i]`). Los productos, restas y transpuestas de
matrices son las operaciones de `Matrix` (`A * B`, `A - B`,
`A.transpose()`), y la identidad es `Identity(n)`. En este tema **no
usamos NumPy**: todos los métodos
se programan desde cero.

## Motivación: ejemplos de física

**Momentos principales de inercia.** Para un cuerpo rígido que rota
alrededor de un eje arbitrario, el tensor de inercia es

```math
I_{\alpha\beta} = \int \rho(\mathbf{r})\left(\delta_{\alpha\beta}\, r^2 - r_\alpha r_\beta\right) d^3r,
\qquad
I = \begin{pmatrix}
I_{xx} & I_{xy} & I_{xz} \\
I_{yx} & I_{yy} & I_{yz} \\
I_{zx} & I_{zy} & I_{zz}
\end{pmatrix},
```

con $\rho$ la densidad de masa y $\alpha, \beta$ componentes
cartesianas. Los elementos de la diagonal son los momentos de inercia
y los de fuera, los productos de inercia; por su definición, esta
matriz sí es simétrica ($I_{xy} = I_{yx}$). Existe un sistema de
coordenadas, los **ejes principales**, en el que los productos de
inercia se anulan y el tensor es diagonal,
$I_P = \mathrm{diag}(I_0, I_1, I_2)$, con $I_0, I_1, I_2$ los
**momentos principales de inercia**. Encontrar los ejes principales
equivale a **diagonalizar** una matriz $3\times 3$: es un problema de
eigenvalores.

**Mecánica cuántica.** La ecuación de Schrödinger independiente del
tiempo, $H\psi = E\psi$, es un problema de eigenvalores: en una base
finita, $H$ es una matriz y las energías son sus eigenvalores. El
proyecto de espines ([`simulaciones/espines/`](../../simulaciones/espines/))
es un ejemplo.

## El problema de eigenvalores

La forma estándar es

```math
A\mathbf{v} = \lambda\mathbf{v}.
```

La diferencia crucial con $A\mathbf{x} = \mathbf{b}$ es que aquí
**tanto** el número $\lambda$ (el **eigenvalor**) **como** el vector
$\mathbf{v}$ (el **eigenvector**) son incógnitas. Pasando todo al lado
izquierdo,

```math
(A - \lambda I)\mathbf{v} = \mathbf{0},
```

con $I$ la matriz identidad $n\times n$. Es un sistema lineal con
coeficientes $A - \lambda I$ y lado derecho cero, pero con $n+1$
incógnitas ($\lambda, v_0, \ldots, v_{n-1}$), así que no podemos
esperar una solución única. La solución trivial $\mathbf{v} = \mathbf{0}$
siempre existe; para que haya una no trivial, la matriz
$A - \lambda I$ tiene que ser **singular**:

```math
|A - \lambda I| = 0.
```

Desarrollando el determinante se obtiene la **ecuación
característica**, un polinomio de grado $n$ en $\lambda$:

```math
(-1)^n\lambda^n + c_{n-1}\lambda^{n-1} + \cdots + c_1\lambda + c_0 = 0.
```

Así que una matriz $n\times n$ tiene a lo más $n$ eigenvalores
distintos: las raíces del polinomio característico. Si una raíz
aparece dos veces decimos que tiene **multiplicidad** (algebraica) 2;
si aparece una sola vez, es un eigenvalor **simple**. (Dos datos
curiosos: la matriz $A$ satisface su propia ecuación característica,
teorema de Cayley–Hamilton; y el producto de los eigenvalores es el
determinante de $A$.)

Conocido un eigenvalor $\lambda_i$, el eigenvector correspondiente se
obtiene resolviendo $(A - \lambda_i I)\mathbf{v}_i = \mathbf{0}$. Como
la matriz es singular, $\mathbf{v}_i$ no queda determinado de forma
única: solo se pueden calcular los valores *relativos* de sus
componentes (por eso normalmente se normaliza,
$`\|\mathbf{v}_i\| = 1`$).

**Ojo con la notación:** $\mathbf{v}_0$ es un vector columna de $n$
elementos (el eigenvector número 0), no una componente. Sus
componentes se escriben $(\mathbf{v}_0)_0, (\mathbf{v}_0)_1, \ldots$:
el primer índice dice de qué eigenvector se trata, el segundo de qué
componente.
## Análisis de error

En la [unidad 09](../09_algebra_lineal/README.md) (sección "Análisis de
error") vimos qué tan sensible es un sistema lineal a cambios en sus
datos, con el número de condición $\kappa(A)$ y las normas de
matrices. Aquí hacemos lo mismo para los eigenvalores y los
eigenvectores.

### Número de condición para eigenvalores simples

Ahora nos interesa el problema $A\mathbf{v} = \lambda\mathbf{v}$. (Uno
podría pensar que $\kappa(A)$ también sirve aquí, pero no: veremos que
la sensibilidad de los eigenvalores se mide con otra cosa.) Con
índices explícitos, para distinguir los eigenvalores:

```math
A\mathbf{v}_i = \lambda_i\mathbf{v}_i.
```

No vamos a suponer que $A$ es simétrica (en física hay casos
importantes que no lo son). Entonces conviene distinguir entre los
**eigenvectores derechos** $\mathbf{v}_i$ de arriba y los
**eigenvectores izquierdos** $\mathbf{u}_i$, definidos por

```math
\mathbf{u}_i^T A = \lambda_i\,\mathbf{u}_i^T.
```

($\mathbf{u}_i^T$ es un vector renglón $1\times n$.) Transponiendo,
$A^T\mathbf{u}_i = \lambda_i\mathbf{u}_i$: **los eigenvectores
izquierdos de $A$ son los eigenvectores derechos de $A^T$**, con los
mismos eigenvalores. Si $A$ no es simétrica, en general
$\mathbf{u}_i \ne \mathbf{v}_i$.

Un dato que usaremos abajo: para eigenvalores distintos,
**$\mathbf{u}_k^T\mathbf{v}_i = 0$ si $k \ne i$**. Demostración:
$\mathbf{u}_k^T A \mathbf{v}_i$ se puede evaluar de dos formas,
$`\lambda_i\,\mathbf{u}_k^T\mathbf{v}_i`$ (usando $A\mathbf{v}_i$) o
$`\lambda_k\,\mathbf{u}_k^T\mathbf{v}_i`$ (usando $\mathbf{u}_k^T A$);
restando, $`(\lambda_k - \lambda_i)\,\mathbf{u}_k^T\mathbf{v}_i = 0`$.

Perturbamos $A$, y con ella cambian eigenvalores y eigenvectores:

```math
(A + \Delta A)(\mathbf{v}_i + \Delta\mathbf{v}_i) = (\lambda_i + \Delta\lambda_i)(\mathbf{v}_i + \Delta\mathbf{v}_i).
```

Desarrollando, cancelando $A\mathbf{v}_i = \lambda_i\mathbf{v}_i$ y
despreciando los términos de segundo orden ($\Delta\times\Delta$):

```math
A\,\Delta\mathbf{v}_i + \Delta A\,\mathbf{v}_i = \lambda_i\,\Delta\mathbf{v}_i + \Delta\lambda_i\,\mathbf{v}_i.
```

Multiplicando por la izquierda por $\mathbf{u}_i^T$, el primer término
de cada lado se cancela (porque
$\mathbf{u}_i^T A = \lambda_i\mathbf{u}_i^T$), y queda

```math
\mathbf{u}_i^T\,\Delta A\,\mathbf{v}_i = \Delta\lambda_i\,\mathbf{u}_i^T\mathbf{v}_i
\quad\Longrightarrow\quad
|\Delta\lambda_i| = \frac{|\mathbf{u}_i^T\,\Delta A\,\mathbf{v}_i|}{|\mathbf{u}_i^T\mathbf{v}_i|}.
```

Por la desigualdad de Cauchy–Schwarz,
$`|\mathbf{u}_i^T\,\Delta A\,\mathbf{v}_i| \le \|\mathbf{u}_i\|\,\|\Delta A\|\,\|\mathbf{v}_i\|`$,
y si normalizamos los eigenvectores, $`\|\mathbf{u}_i\| = \|\mathbf{v}_i\| = 1`$
(como hacen las bibliotecas estándar):

```math
|\Delta\lambda_i| \le \frac{1}{|\mathbf{u}_i^T\mathbf{v}_i|}\,\|\Delta A\|.
```

(Aquí $`\|\mathbf{u}\|`$ es la norma euclídea; $`\|\Delta A\|`$ es la norma
de matrices asociada a ella, que siempre es menor o igual que la de
Frobenius, así que la cota vale también con $`\|\Delta A\|_F`$.)

El factor que amplifica la perturbación es el **número de condición
para el eigenvalor simple** $\lambda_i$:

```math
\kappa^{ev}_{\lambda_i}(A) = \frac{1}{|\mathbf{u}_i^T\mathbf{v}_i|}.
```

El subíndice $ev$ recuerda que es para eigenvalores; el superíndice,
de cuál eigenvalor se trata (cada uno tiene el suyo). Si
$\kappa^{ev} \approx 1$ el eigenvalor está bien condicionado; si
$\kappa^{ev} \gg 1$, mal condicionado. Tres observaciones:

- Aquí la cota es para el error **absoluto** $|\Delta\lambda_i|$, no el
  relativo: tiene sentido, porque un eigenvalor puede valer cero (si la
  matriz es singular).
- No hizo falta $A^{-1}$ en ningún momento.
- Para una matriz **simétrica**, $A^T = A$, así que
  $\mathbf{u}_i = \mathbf{v}_i$ y $\kappa^{ev} = 1$: **los eigenvalores
  de matrices simétricas siempre están bien condicionados.** Los
  problemas aparecen con matrices no simétricas, cuyos eigenvectores
  izquierdo y derecho pueden ser casi perpendiculares.

Ejemplo (en `condicion_eigenvalores.py`): para
$`A = \begin{pmatrix}1 & 1000\\ 0 & 2\end{pmatrix}`$, los eigenvalores
son 1 y 2, y ambos tienen $\kappa^{ev} \approx 1000$. Sumar
$10^{-6}$ a $A_{10}$ mueve cada eigenvalor $\approx 10^{-3}$.

### Sensibilidad de los eigenvectores

¿Y los eigenvectores? (Tampoco basta $\kappa^{ev}$.) Suponemos
eigenvalores distintos, así que los eigenvectores son linealmente
independientes. Partimos de la misma ecuación de primer orden:

```math
A\,\Delta\mathbf{v}_i + \Delta A\,\mathbf{v}_i = \lambda_i\,\Delta\mathbf{v}_i + \Delta\lambda_i\,\mathbf{v}_i,
```

y desarrollamos la perturbación del eigenvector en términos de los
*demás* eigenvectores:

```math
\Delta\mathbf{v}_i = \sum_{j\ne i} t_{ji}\,\mathbf{v}_j,
```

con coeficientes $t_{ji}$ por determinar. (No hace falta un término
$j = i$: una componente de $\Delta\mathbf{v}_i$ a lo largo de
$\mathbf{v}_i$ solo cambia su normalización.) Sustituyendo, y usando
$A\mathbf{v}_j = \lambda_j\mathbf{v}_j$:

```math
\sum_{j\ne i}(\lambda_j - \lambda_i)\,t_{ji}\,\mathbf{v}_j + \Delta A\,\mathbf{v}_i = \Delta\lambda_i\,\mathbf{v}_i.
```

Multiplicamos por la izquierda por $\mathbf{u}_k^T$ con $k \ne i$.
Como $\mathbf{u}_k^T\mathbf{v}_j = 0$ para $j \ne k$, de la suma solo
sobrevive el término $j = k$, y el lado derecho se anula:

```math
(\lambda_k - \lambda_i)\,t_{ki}\,\mathbf{u}_k^T\mathbf{v}_k + \mathbf{u}_k^T\,\Delta A\,\mathbf{v}_i = 0.
```

Despejando $t_{ki}$ y sustituyendo en el desarrollo de
$\Delta\mathbf{v}_i$:

```math
\Delta\mathbf{v}_i = \sum_{j\ne i}\frac{\mathbf{u}_j^T\,\Delta A\,\mathbf{v}_i}{(\lambda_i - \lambda_j)\,\mathbf{u}_j^T\mathbf{v}_j}\,\mathbf{v}_j.
```

Este es el resultado principal. Nótese que:

- A diferencia de los casos anteriores, es una **suma**: la
  perturbación de un eigenvector tiene componentes a lo largo de todos
  los demás.
- El numerador contiene la perturbación $\Delta A$.
- El denominador tiene **dos** fuentes de amplificación:
  1. $\mathbf{u}_j^T\mathbf{v}_j$, lo mismo que en $\kappa^{ev}$: si
     algún eigenvalor está mal condicionado, los eigenvectores también.
  2. $\lambda_i - \lambda_j$, la **separación entre eigenvalores**: si
     dos eigenvalores están muy cerca, sus eigenvectores son muy
     sensibles, **aunque la matriz sea simétrica**.

Lo segundo generaliza algo que ya sabíamos: si dos eigenvalores
coinciden, los eigenvectores no quedan determinados de forma única
(cualquier combinación de ellos sirve). Si solo están cerca, quedan
determinados, pero muy mal.

Ejemplo (en `condicion_eigenvalores.py`): la matriz simétrica
$`\begin{pmatrix}1 & \epsilon\\ \epsilon & 1+\delta\end{pmatrix}`$ con
$\epsilon = 10^{-7}$. Si $\delta = 1$, el eigenvector de $\lambda \approx 1$
gira unos $10^{-5}$ grados respecto a $(1,0)$; si $\delta = 10^{-9}$,
gira más de 40 grados, mientras que el eigenvalor se mueve
apenas $10^{-7}$.

## El problema de eigenvalores, numéricamente

Pasamos a la "segunda mitad" del álgebra lineal:

```math
A\mathbf{v}_i = \lambda_i\mathbf{v}_i .
```

El truco principal de la unidad 09 ya no sirve: restarle
a un renglón un múltiplo de otro (la eliminación gaussiana) **cambia
los eigenvalores** de la matriz. Así que aquí no vamos a eliminar.

Nos restringimos al caso en el que la matriz $A$ de $n\times n$ tiene
$n$ eigenvalores **distintos**. Entonces sus $n$ eigenvectores son
linealmente independientes, y si los ponemos como columnas de una
matriz $V$ (la *matriz de eigenvectores*) se cumple

```math
V^{-1} A V = \Lambda,
\qquad
\Lambda = \begin{pmatrix}
\lambda_0 & 0 & \cdots & 0\\
0 & \lambda_1 & \cdots & 0\\
\vdots & \vdots & \ddots & \vdots\\
0 & 0 & \cdots & \lambda_{n-1}
\end{pmatrix},
\qquad
V = \begin{pmatrix} \mathbf{v}_0 & \mathbf{v}_1 & \cdots & \mathbf{v}_{n-1}\end{pmatrix}.
```

(La prueba es el problema 4.21 de [1]: escriban $AV$ columna por
columna.) Por eso a resolver el problema de eigenvalores también se le
llama **diagonalizar** la matriz, y a $A = V\Lambda V^{-1}$ se le
llama la **eigendescomposición**.

No vamos a suponer que $A$ es simétrica: muchos problemas de física
dan matrices simétricas, pero no todos. Eso sí, una matriz no
simétrica puede tener eigenvalores complejos; aquí solo veremos
matrices no simétricas con eigenvalores reales (y distintos). La
matriz de prueba de [1], `crear_prueba(4, 21)`, es de ese tipo.

**¿Por qué no usar el polinomio característico?** En la sección "El
problema de eigenvalores" vimos que los eigenvalores son las raíces de
$|A - \lambda I| = 0$. Pero encontrar las raíces de un polinomio suele
ser un problema **mal condicionado** (lo veremos en el capítulo 5),
aunque el problema de eigenvalores original esté bien condicionado. El
proyecto de espines ya se topó con otra limitación: barriendo
$\lambda$ y buscando cambios de signo del determinante, un eigenvalor
repetido no cambia el signo y se pierde. Es más sensato transformar la
matriz en otra de la que se puedan leer los eigenvalores, cuidando que
sean los mismos.

Todos los métodos para eigenvalores son **iterativos**. El que se usa
en la práctica cuando se quieren *todos* los eigenvalores es el
**método QR**. Antes de llegar a él veremos dos métodos más sencillos,
el de la **potencia** y el de la **potencia inversa**, que dan un solo
eigenpar (un eigenvalor y su eigenvector) cada uno. Sirven
pedagógicamente, se parecen conceptualmente a QR, y la potencia
inversa será una pieza del programa final.

## Método de la potencia

Como los $n$ eigenvalores son distintos, los podemos ordenar por
magnitud:

```math
|\lambda_0| > |\lambda_1| > |\lambda_2| > \cdots > |\lambda_{n-1}|.
```

El método de la potencia encuentra $\lambda_0$, el eigenvalor
**dominante**, y su eigenvector $\mathbf{v}_0$.

### Primera versión

La receta es muy simple: se escoge un vector inicial
$`\mathbf{z}^{(0)}`$ y se multiplica por $A$ una y otra vez (el
superíndice entre paréntesis es el número de iteración, como en
Jacobi, unidad 09):

```math
\mathbf{z}^{(k)} = A\mathbf{z}^{(k-1)}, \qquad k = 1, 2, \ldots
```

Entonces $`\mathbf{z}^{(1)} = A\mathbf{z}^{(0)}`$,
$`\mathbf{z}^{(2)} = A^2\mathbf{z}^{(0)}`$, y en general
$`\mathbf{z}^{(k)} = A^k\mathbf{z}^{(0)}`$: aparecen las *potencias*
de $A$, de ahí el nombre.

¿Qué tiene que ver esto con los eigenvalores? Escribamos el vector
inicial como combinación lineal de los eigenvectores (que no
conocemos; los coeficientes $c_i$ tampoco):

```math
\mathbf{z}^{(0)} = \sum_{i=0}^{n-1} c_i\mathbf{v}_i .
```

Multiplicando por $A^k$ y usando $k$ veces que
$A\mathbf{v}_i = \lambda_i\mathbf{v}_i$:

```math
\mathbf{z}^{(k)} = A^k\mathbf{z}^{(0)}
= \sum_{i=0}^{n-1} c_i\lambda_i^k\mathbf{v}_i
= \lambda_0^k\left[c_0\mathbf{v}_0 + \sum_{i=1}^{n-1} c_i\left(\frac{\lambda_i}{\lambda_0}\right)^k\mathbf{v}_i\right].
```

En el último paso sacamos el término $i = 0$ de la suma y factorizamos
$\lambda_0^k$. Como $\lambda_0$ es el más grande en magnitud, cada
$`(\lambda_i/\lambda_0)^k \to 0`$ cuando $k \to \infty$. Así que, si
$c_0 \neq 0$ (es decir, si el vector inicial tiene alguna componente
en la dirección de $\mathbf{v}_0$), $`\mathbf{z}^{(k)}`$ apunta cada vez
más en la dirección de $\mathbf{v}_0$.

¿Qué tan rápido? El término que tarda más en desaparecer es el de
$\lambda_1$, el segundo más grande. El error se reduce en cada
iteración más o menos por un factor

```math
\left|\frac{\lambda_1}{\lambda_0}\right| .
```

Si $\lambda_1$ y $\lambda_0$ son parecidos, la convergencia es lenta.

**El eigenvalor.** Sabemos que $`\mathbf{z}^{(k)}`$ tiende a
$c_0\lambda_0^k\mathbf{v}_0$, pero no conocemos ni $c_0$, ni
$\lambda_0$, ni $\mathbf{v}_0$. Aun así basta para sacar el eigenvalor
con el **cociente de Rayleigh** de un vector $\mathbf{x}$:

```math
\mu(\mathbf{x}) = \frac{\mathbf{x}^T A\mathbf{x}}{\mathbf{x}^T\mathbf{x}} .
```

Si $\mathbf{x}$ es un eigenvector, $`A\mathbf{x} = \lambda\mathbf{x}`$
y $`\mu(\mathbf{x}) = \lambda`$ exactamente (cualquier factor que
multiplique a $\mathbf{x}$ se cancela arriba y abajo). Si
$\mathbf{x}$ no es un eigenvector, $`\mu(\mathbf{x})`$ es lo más
parecido a un eigenvalor que se le puede asociar (en el sentido de
mínimos cuadrados, capítulo 6). Así que $`\mu(\mathbf{z}^{(k)})`$ es
nuestra mejor estimación de $\lambda_0$ en la iteración $k$.

### Normalizando

Hay un problema que ignoramos: el factor $\lambda_0^k$. Si
$|\lambda_0| > 1$ crece sin límite (desbordamiento), y si
$|\lambda_0| < 1$ se va a cero (subdesbordamiento). La solución es
reescalar el vector en cada paso para que tenga norma 1 (norma
euclídea):

```math
\mathbf{z}^{(k)} = A\mathbf{q}^{(k-1)},
\qquad
\mathbf{q}^{(k)} = \frac{\mathbf{z}^{(k)}}{\|\mathbf{z}^{(k)}\|},
\qquad
\mu(\mathbf{q}^{(k)}) = [\mathbf{q}^{(k)}]^T A\mathbf{q}^{(k)},
\qquad k = 1, 2, \ldots
```

empezando con un $`\mathbf{q}^{(0)}`$ de norma 1. Como
$`\|\mathbf{q}^{(k)}\| = 1`$, el denominador del cociente de Rayleigh
vale 1. Ahora $`\mathbf{q}^{(k)}`$ ya no es *igual* a
$`A^k\mathbf{q}^{(0)}`$ sino *proporcional*: seguimos multiplicando por
$A$ en cada paso, pero quitamos el factor que sobra. De paso, el
eigenvector ya sale normalizado: es el propio $`\mathbf{q}^{(k)}`$.
(Los nombres $\mathbf{q}$ y $\mathbf{z}$ no son casualidad: los
volveremos a usar en el método QR.)

### Implementación

En [`fiscomp/algebra_lineal.py`](../../fiscomp/algebra_lineal.py)
agregamos tres funciones auxiliares, `producto_punto`, `normalizar`
(divide entre `norma_euclidea`) y `cociente_rayleigh`, y el método:

```python
def potencia(A, kmax=6):
    q = normalizar([1.0] * A.rows)
    for k in range(1, kmax):
        z = mat_vec(A, q)
        q = normalizar(z)
    return cociente_rayleigh(A, q), q
```

El vector inicial es $`(1, 1, \ldots, 1)`$
normalizado, se hace un número fijo de iteraciones, y el cociente de
Rayleigh se calcula **una sola vez**, al final: calcularlo en cada
iteración sería desperdiciar un producto matriz-vector por vuelta.

La primera parte de [`potencia.py`](potencia.py) lo prueba con la
matriz de prueba y reproduce el resultado de [1],
$\lambda_0 = 21.3166626635$ y
$`\mathbf{v}_0 = (0.44439562, 0.48218122, 0.51720906, 0.55000781)`$.
Como aquí no usamos NumPy (y por lo tanto tampoco `numpy.linalg.eig`),
comprobamos con el residuo
$`\|A\mathbf{v} - \lambda\mathbf{v}\|/\|\mathbf{v}\|`$, que sale
$\sim 5\times 10^{-13}$. El script también imprime el error relativo
de $`\mu(\mathbf{q}^{(k)})`$ en cada iteración: el cociente entre dos
errores seguidos es $\approx 3.15\times 10^{-3}$, que es exactamente
$|\lambda_1/\lambda_0| = 0.0672/21.3$. Por eso converge tan rápido:
cada iteración gana unos dos dígitos y medio, y en 5 iteraciones ya
estamos en el épsilon de la máquina.

En vez de escoger el número de iteraciones a ojo, sería mejor un
criterio de paro como el de Jacobi, con `cambio_relativo` entre dos
$\mathbf{q}$ seguidos (eso lo hace la potencia inversa, abajo).

### Conteo de operaciones

El número total de operaciones de un método iterativo depende de
cuántas iteraciones $m$ hagan falta, que no se sabe de antemano
(depende del vector inicial, de la tolerancia, del criterio de paro y,
sobre todo, de $|\lambda_1/\lambda_0|$). Por iteración, el producto
$A\mathbf{q}$ cuesta $\sim 2n^2$ y la norma $\sim 2n$; al final, el
cociente de Rayleigh cuesta otro $\sim 2n^2$. En total, el método de
la potencia cuesta $`\sim 2(m+1)n^2`$ operaciones.

## Método de la potencia inversa con desplazamiento

Ahora queremos el eigenvalor de **menor** magnitud, $`\lambda_{n-1}`$.

### La potencia inversa

Multiplicando $A\mathbf{v}_i = \lambda_i\mathbf{v}_i$ por $A^{-1}$ por
la izquierda y dividiendo entre $\lambda_i$:

```math
A^{-1}\mathbf{v}_i = \frac{1}{\lambda_i}\mathbf{v}_i .
```

Es decir, $A^{-1}$ tiene **los mismos eigenvectores** que $A$, y sus
eigenvalores son los **recíprocos**. (Si algún $\lambda_i = 0$, $A$ es
singular y no hay $A^{-1}$; recuerden que el determinante es el
producto de los eigenvalores.) El eigenvalor más grande de $A^{-1}$ es
$`1/\lambda_{n-1}`$, así que aplicar el método de la potencia a
$A^{-1}$,

```math
\mathbf{z}^{(k)} = A^{-1}\mathbf{q}^{(k-1)},
```

hace que $`\mathbf{q}^{(k)}`$ tienda a $`\mathbf{v}_{n-1}`$. Y como los
eigenvectores son los mismos, al final podemos usar el cociente de
Rayleigh con la $A$ original, $`[\mathbf{q}^{(k)}]^T A\mathbf{q}^{(k)}`$,
que da directamente $`\lambda_{n-1}`$ (con $A^{-1}$ daría
$`1/\lambda_{n-1}`$).

Pero calcular $A^{-1}$ es caro, y no hace falta. Multiplicando por
$A$ por la izquierda, cada paso es un **sistema lineal**:

```math
A\mathbf{z}^{(k)} = \mathbf{q}^{(k-1)} .
```

Es siempre la misma matriz $A$ con distintos lados derechos: el caso
ideal para LU (unidad 09, sección "Resolver un sistema con LU"). Se factoriza
$A = LU$ **una vez**, fuera del ciclo, y en cada iteración solo se
hacen las dos sustituciones con `resolver_lu`. El método queda casi
idéntico al de la potencia; la única diferencia (crucial) es que donde
había un producto matriz-vector ahora hay un sistema que resolver:

```math
A\mathbf{z}^{(k)} = \mathbf{q}^{(k-1)},
\qquad
\mathbf{q}^{(k)} = \frac{\mathbf{z}^{(k)}}{\|\mathbf{z}^{(k)}\|},
\qquad
\mu(\mathbf{q}^{(k)}) = [\mathbf{q}^{(k)}]^T A\mathbf{q}^{(k)} .
```

### Desplazamiento de eigenvalores

Un refinamiento más. Restemos $s\mathbf{v}_i$ de los dos lados de
$A\mathbf{v}_i = \lambda_i\mathbf{v}_i$, con $s$ un número cualquiera
(el **desplazamiento**, *shift*):

```math
(A - sI)\mathbf{v}_i = (\lambda_i - s)\mathbf{v}_i .
```

La matriz $`A^* = A - sI`$ tiene los mismos eigenvectores que $A$, y sus
eigenvalores son los de $A$ recorridos: $`\lambda^*_i = \lambda_i - s`$.

Si aplicamos la potencia inversa a $`A^*`$, encontramos su eigenvalor
de menor magnitud, el $`\lambda^*_i = \lambda_i - s`$ más cercano a
cero. Eso es lo mismo que encontrar el eigenvalor $\lambda_i$ de $A$
**más cercano a** $s$. Esa es la gracia de combinar las dos ideas: la
parte "inversa" encuentra el más chico, y el desplazamiento decide qué
quiere decir "más chico". Con $s = 0$ recuperamos la potencia inversa
simple.

Cambiando $s$ podemos ir encontrando distintos eigenpares, aunque sin
información previa no es una manera eficiente de encontrarlos todos
(¿qué valores de $s$ probamos? ¿una malla?). El uso más importante es
otro: si ya tenemos una **estimación** de un eigenvalor, la potencia
inversa con ese desplazamiento la refina. Y si ya conocemos un
eigenvalor con mucha precisión pero no su eigenvector, unas cuantas
iteraciones con $s \approx \lambda_i$ dan $\mathbf{v}_i$. Así lo
usaremos después de QR, que da los eigenvalores pero no los
eigenvectores. (La primera versión del proyecto de espines ya lo hacía,
para refinar los eigenvalores que encontraba con el determinante.)

**Rapidez.** Por el mismo argumento que en la potencia directa, la
potencia inversa converge según el cociente entre el eigenvalor de
$`A^*`$ que buscamos y el siguiente más chico. Para $s = 0$ es
$`|\lambda_{n-1}/\lambda_{n-2}|`$. Con un desplazamiento $s$ muy cerca
del eigenvalor buscado, $`\lambda^*_i = \lambda_i - s`$ es diminuto y el
cociente también: el desplazamiento **acelera** la convergencia. (No
hay problema en que $A - sI$ sea casi singular: el error de redondeo
al resolver el sistema cae, justamente, en la dirección del
eigenvector que buscamos.)

### Implementación

```python
def potencia_inversa(A, desplazamiento=20.0, kmax=200, tol=1e-8):
    A_desplazada = A - Identity(A.rows) * desplazamiento
    L, U = descomposicion_lu(A_desplazada)
    q = normalizar([1.0] * A.rows)
    for k in range(1, kmax):
        z = resolver_lu(L, U, q)
        q_nuevo = normalizar(z)
        if producto_punto(q, q_nuevo) < 0:
            q_nuevo = [-q_i for q_i in q_nuevo]
        error = cambio_relativo(q, q_nuevo)
        q = q_nuevo
        if error < tol:
            break
    else:
        return None, None, kmax
    return cociente_rayleigh(A, q), q, k
```

Comparada con `potencia`:

- Al principio se construye $A - sI$ y se factoriza con
  `descomposicion_lu`, **una sola vez**. Dentro del ciclo solo se usa
  `resolver_lu` (las dos sustituciones), que es mucho más barato.
- Tiene criterio de paro: el mismo `cambio_relativo` y el mismo
  `for`-`else` de `jacobi`. Si converge regresa
  `(eigenvalor, eigenvector, iteraciones)`; si no,
  `(None, None, kmax)`. Igual que en Jacobi, `q` y `q_nuevo` son
  listas distintas (`normalizar` crea una lista nueva); con NumPy habría
  que usar `np.copy` para lograrlo.
- **El signo.** Si el eigenvalor de $A - sI$ que domina es negativo,
  cada iteración multiplica por un número negativo y $\mathbf{q}$
  cambia de signo en cada vuelta. Un vector y su opuesto son el mismo
  eigenvector, pero para `cambio_relativo` son completamente
  distintos (cada componente cambia en $2$, la suma da $2n$) y el
  criterio de paro nunca se cumpliría. Por eso, si
  $`\mathbf{q}^{(k-1)}\cdot\mathbf{q}^{(k)} < 0`$ (apuntan en
  sentidos opuestos), se le cambia el signo a $`\mathbf{q}^{(k)}`$.
  `potencia` no lo necesita porque no tiene criterio de paro.
- El eigenvalor se calcula con el cociente de Rayleigh de la $A$
  **original**, así que sale directamente $\lambda_i$ (con $`A^*`$
  saldría $\lambda_i - s$ y habría que sumarle $s$).
- Un detalle de `cambio_relativo`: divide entre cada componente, y
  algunas componentes de un eigenvector pueden ser cero. Para los
  ejemplos de aquí no pasa, pero es un criterio frágil.

[`potencia.py`](potencia.py) lo prueba con la matriz de prueba:

- Con $s = 20$ encuentra el mismo eigenpar que `potencia`, en 8
  iteraciones (no las 200 de `kmax`).
- Probando desplazamientos encuentra los **cuatro** eigenvalores:
  $21.3$, $-0.0672$, $`-6.93\times 10^{-5}`$ y $`-7.50\times 10^{-8}`$.
  Se comprueban sin NumPy con dos datos de la sección "El problema de
  eigenvalores": su suma coincide con la traza de $A$ y su producto
  con el determinante (calculado con LU), $`-7.45\times 10^{-12}`$.
  (Tres eigenvalores diminutos: la matriz está muy cerca de ser
  singular, como ya nos decía $`\kappa \approx 4\times 10^{8}`$ en la unidad 09.)
- Con $s = 1$ **no converge** en 200 iteraciones: los dos eigenvalores
  más cercanos a 1 están casi a la misma distancia de él, y el
  cociente que controla la rapidez es $0.99993$.
- Buscando $\lambda = -0.0672$ con desplazamientos cada vez más
  cercanos, las iteraciones bajan de 104 ($s = -0.5$) a 4
  ($s = -0.0672$).
- Con $-A$ (eigenvalor dominante negativo) muestra el cambio de signo
  de $\mathbf{q}$ en cada iteración.

### Conteo de operaciones

La descomposición LU de $`A^*`$ cuesta $\sim 2n^3/3$, una sola vez
(construir $A - sI$ solo toca la diagonal: $`O(n)`$). Por iteración, las
dos sustituciones cuestan $n^2 + n^2$ y la norma $\sim 2n$; al final,
el cociente de Rayleigh $\sim 2n^2$. En total,

```math
\sim 2(m+1)n^2 + \frac{2n^3}{3} .
```

El $n^3$ viene únicamente de LU: si la matriz ya estaba factorizada,
la potencia inversa cuesta lo mismo que la directa. (Y si $m \approx n$,
cosa común, los dos métodos son $`O(n^3)`$.)

**Para pensar.**

- Agreguen a `potencia` un criterio de paro con `cambio_relativo`,
  como en `potencia_inversa` (problema 4.23 de [1]). Prueben con
  `crear_prueba(20, 21)`: según [1], bastan 6 iteraciones.
  ¿Qué pasa con $-A$ si no corrigen el signo?
- Repitan el desarrollo de $`\mathbf{z}^{(k)} = A^k\mathbf{z}^{(0)}`$
  para la versión normalizada y muestren que $`\mathbf{q}^{(k)}`$ es
  proporcional a $`A^k\mathbf{q}^{(0)}`$ (problema 4.22).
- ¿Qué pasa con el método de la potencia si el vector inicial no
  tiene componente en $\mathbf{v}_0$ ($c_0 = 0$)? Construyan un
  ejemplo $2\times 2$ y córranlo. (Pista: el redondeo ayuda.)

## Método QR

La potencia y la potencia inversa dan **un** eigenpar a la vez. Con
desplazamientos podemos ir buscando los demás, pero hay que saber
dónde buscar. El **método QR** (J. Francis, principios de los 60)
calcula **todos** los eigenvalores de una sola vez, y es la base de lo
que hacen las bibliotecas profesionales.

Ojo con los nombres, que se parecen mucho:

- La **descomposición QR** escribe una matriz como producto de otras
  dos, $A = QR$. Es un método **directo**, como LU: un número fijo de
  pasos que depende solo de $n$.
- El **método QR** calcula los eigenvalores. Es **iterativo** (como
  todos los métodos de eigenvalores) y en cada iteración hace una
  descomposición QR.

### Descomposición QR

Toda matriz cuadrada real (con columnas linealmente independientes)
se puede escribir como

```math
A = QR ,
```

con $Q$ **ortogonal** y $R$ **triangular superior** (se llama $R$ y no
$U$ por razones históricas: *right*). Ortogonal quiere decir que
$`Q^{-1} = Q^T`$, o sea $`Q^TQ = I`$: las columnas de $Q$ son vectores
**ortonormales** (ortogonales entre sí y de norma 1).

La prueba es **constructiva**: construimos $Q$ a partir de $A$, y en el
camino nos sale $R$ gratis.

**Construyendo Q: Gram-Schmidt.** Escribimos $A$ por columnas,

```math
A = \begin{pmatrix} \mathbf{a}_0 & \mathbf{a}_1 & \cdots & \mathbf{a}_{n-1}\end{pmatrix},
\qquad
Q = \begin{pmatrix} \mathbf{q}_0 & \mathbf{q}_1 & \cdots & \mathbf{q}_{n-1}\end{pmatrix},
```

y buscamos vectores $\mathbf{q}_j$ ortonormales. El primero es fácil:
la dirección de $\mathbf{a}_0$, normalizada,

```math
\mathbf{q}_0 = \frac{\mathbf{a}_0}{\|\mathbf{a}_0\|} .
```

Para el segundo no basta normalizar $\mathbf{a}_1$: casi seguro no es
perpendicular a $\mathbf{a}_0$. La receta de **Gram-Schmidt** es
quitarle a $\mathbf{a}_1$ su componente en la dirección de
$\mathbf{q}_0$ y luego normalizar:

```math
\mathbf{a}'_1 = \mathbf{a}_1 - (\mathbf{q}_0^T\mathbf{a}_1)\,\mathbf{q}_0,
\qquad
\mathbf{q}_1 = \frac{\mathbf{a}'_1}{\|\mathbf{a}'_1\|} .
```

(Comprueben que $`\mathbf{q}_0^T\mathbf{a}'_1 = 0`$.) Para el tercero
quitamos las componentes en las direcciones de $\mathbf{q}_0$ **y** de
$\mathbf{q}_1$, y así sucesivamente. En general, para
$j = 0, 1, \ldots, n-1$:

```math
\mathbf{a}'_j = \mathbf{a}_j - \sum_{i=0}^{j-1} (\mathbf{q}_i^T\mathbf{a}_j)\,\mathbf{q}_i,
\qquad
\mathbf{q}_j = \frac{\mathbf{a}'_j}{\|\mathbf{a}'_j\|} .
```

Para $j = 0$ la suma no tiene términos y recuperamos la fórmula de
$\mathbf{q}_0$.

**Encontrando R.** Supongamos que $A = QR$ y escribámoslo columna por
columna. Como $R$ es triangular superior, la columna $j$ de $A$ solo
usa las primeras $j+1$ columnas de $Q$:

```math
\begin{aligned}
\mathbf{a}_0 &= R_{00}\,\mathbf{q}_0 \\
\mathbf{a}_1 &= R_{01}\,\mathbf{q}_0 + R_{11}\,\mathbf{q}_1 \\
&\;\;\vdots \\
\mathbf{a}_j &= R_{0j}\,\mathbf{q}_0 + R_{1j}\,\mathbf{q}_1 + \cdots + R_{jj}\,\mathbf{q}_j .
\end{aligned}
```

Despejando $\mathbf{q}_j$ de la última,

```math
\mathbf{q}_j = \frac{\mathbf{a}_j - \sum_{i=0}^{j-1} R_{ij}\,\mathbf{q}_i}{R_{jj}},
```

y comparando con Gram-Schmidt se lee qué tiene que ser $R$:

```math
R_{ij} = \mathbf{q}_i^T\mathbf{a}_j \quad (i < j),
\qquad
R_{jj} = \|\mathbf{a}'_j\| ,
\qquad
R_{ij} = 0 \quad (i > j).
```

Lo importante: **los dos ya los calculamos** al construir $Q$. Así que
$Q$ y $R$ se van llenando juntos. (Escogimos $R_{jj} > 0$; con esa
convención la descomposición QR es única.)

**De paso.** Gram-Schmidt convierte $n$ vectores linealmente
independientes en $n$ vectores ortonormales. La misma idea funciona
con *funciones* si el producto punto se cambia por la integral del
producto de dos funciones (el espacio de Hilbert de mecánica
cuántica). Ortonormalizando los monomios $1, x, x^2, \ldots$ en
$[-1, 1]$ salen (múltiplos de) los **polinomios de Legendre**
(problema 4.27 de [1]).

### Implementación

Una `Matrix` guarda sus elementos por **renglones** y Gram-Schmidt
trabaja con **columnas**, así que usamos la transpuesta: las columnas
de $A$ son los renglones de $A^T$ (`A.transpose().data`). En
[`fiscomp/algebra_lineal.py`](../../fiscomp/algebra_lineal.py):

```python
def descomposicion_qr(A):
    n = A.rows
    columnas = A.transpose().data
    q = []  # columnas de Q, ya ortonormales
    R = [[0.0] * n for _ in range(n)]
    for j in range(n):
        a_prima = columnas[j][:]
        for i in range(j):
            R[i][j] = producto_punto(q[i], columnas[j])
            a_prima = [a - R[i][j] * q_i for a, q_i in zip(a_prima, q[i])]
        R[j][j] = norma_euclidea(a_prima)
        q.append([a / R[j][j] for a in a_prima])
    return Matrix(q).transpose(), Matrix(R)
```

Es una traducción directa de las fórmulas: `a_prima` empieza siendo
una copia de $\mathbf{a}_j$ y se le van restando las componentes; el
ciclo interno no se ejecuta para $j = 0$. Al final, `q` es la lista de
columnas de $Q$: armamos con ella una `Matrix` y la transponemos para
que las columnas queden como columnas. `R` se llena como lista de
listas y al final se convierte en `Matrix`. No hay `kmax` ni tolerancia: es un método directo.

¿Cómo sabemos si funciona? Con dos pruebas. La obvia,
$`\|A - QR\|`$, que debe ser cero. Y otra: en aritmética exacta $Q$ es
ortogonal por construcción, pero con redondeo hay que comprobarlo con
$`\|Q^TQ - I\|`$. La parte 1 de [`qr.py`](qr.py) lo hace con la matriz
de prueba (norma de Frobenius):

| $n$ | $`\|A - QR\|`$ | $`\|Q^TQ - I\|`$ |
|---|---|---|
| 4 | $`8.9\times 10^{-16}`$ | $`7.3\times 10^{-3}`$ |
| 6 | $0$ | $2.4$ |
| 8 | $`1.5\times 10^{-15}`$ | $4.5$ |

$QR$ reconstruye $A$ al nivel del épsilon de la máquina, pero $Q$
**deja de ser ortogonal** muy rápido: ya para $n = 6$ no se parece en
nada a una matriz ortogonal. Gram-Schmidt tal como lo escribimos (el
**Gram-Schmidt clásico**) se porta mal con el redondeo, sobre todo con
una matriz casi singular como esta: si dos columnas son casi
paralelas, $`\mathbf{a}'_j`$ es la resta de dos vectores casi iguales
(cancelación catastrófica, unidad 06). Hay una variante, el
**Gram-Schmidt modificado**, que se porta mucho mejor (problema 4.25),
y otra manera de hacer QR con reflexiones de Householder (problema
4.28). Para el método QR esta falta de ortogonalidad resulta no ser
tan grave.

### Transformaciones de similitud

Diagonalizar es encontrar $V$ con $V^{-1}AV = \Lambda$. Hagamos algo
más modesto: tomemos cualquier matriz invertible $S$ y definamos

```math
A' = S^{-1} A S .
```

A esto se le llama una **transformación de similitud** y se dice que
$A$ y $`A'`$ son **similares**. No pedimos que $`A'`$ sea diagonal. ¿Para
qué sirve? Sustituyendo $`A = SA'S^{-1}`$ en
$A\mathbf{v}_i = \lambda_i\mathbf{v}_i$ y multiplicando por $S^{-1}$
por la izquierda:

```math
A'\,(S^{-1}\mathbf{v}_i) = \lambda_i\,(S^{-1}\mathbf{v}_i) .
```

Es decir, **$A$ y $`A'`$ tienen los mismos eigenvalores**. Los
eigenvectores no son los mismos: los de $`A'`$ son
$`\mathbf{v}'_i = S^{-1}\mathbf{v}_i`$. (También tienen el mismo
determinante, problema 4.29.)

Si $S$ es una matriz ortogonal $Q$, la inversa es gratis
($`Q^{-1} = Q^T`$, no hay que resolver $n$ sistemas):

```math
A' = Q^T A Q .
```

La idea del método QR: si encontramos una $Q$ que haga a $`A'`$
**triangular**, los eigenvalores de $A$ se leen de la diagonal de
$`A'`$ (los eigenvalores de una matriz triangular son los elementos de
su diagonal). A $`A = QA'Q^T`$ con $`A'`$ triangular se le llama la
**descomposición de Schur**. Falta ver cómo encontrar esa $Q$ (la $Q$
de $A = QR$ no sirve, pruébenlo).

**¿Y LU?** $U$ también es triangular, y el producto de su diagonal da
el determinante. ¿Por qué no leer los eigenvalores de la diagonal de
$U$? Porque $A = LU$ **no** es una transformación de similitud: los
eigenvalores de $U$ no son los de $A$. La parte 2 de
[`qr.py`](qr.py) lo muestra: la diagonal de $U$ es
$`(4.58, -0.0186, -8.4\times 10^{-5}, -1.0\times 10^{-6})`$ y los
eigenvalores son $`(21.3, -0.0672, -6.9\times 10^{-5}, -7.5\times 10^{-8})`$.
Solo coinciden los productos.

### Iteración simultánea: primer intento

El método de la potencia sigue **un** vector: lo multiplica por $A$
una y otra vez, y el vector se va alineando con el eigenvector
dominante $\mathbf{v}_0$. La idea ahora es seguir $n$ vectores a la
vez, con la esperanza de sacar los $n$ eigenvectores.

En vez de un vector inicial $`\mathbf{z}^{(0)}`$ usamos $n$, y los
acomodamos como columnas de una matriz:

```math
Z^{(0)} = \begin{pmatrix} \mathbf{z}_0^{(0)} & \mathbf{z}_1^{(0)} & \cdots & \mathbf{z}_{n-1}^{(0)}\end{pmatrix} .
```

La elección natural es la identidad, $`Z^{(0)} = I`$: sus columnas
$`\mathbf{e}_0, \ldots, \mathbf{e}_{n-1}`$ son linealmente
independientes. Multiplicar $A$ por una matriz es multiplicar $A$ por
cada una de sus columnas (la columna $j$ de $AZ$ es
$`A\mathbf{z}_j`$), así que hacer la potencia con las $n$ columnas a
la vez es un solo producto de matrices:

```math
Z^{(k)} = A Z^{(k-1)}, \qquad k = 1, 2, \ldots
```

Aplicándolo $k$ veces, $`Z^{(k)} = A A \cdots A\, Z^{(0)} = A^k Z^{(0)}`$,
y como $`Z^{(0)} = I`$, simplemente $`Z^{(k)} = A^k`$.

(Cuidado con la notación: un superíndice **entre paréntesis** es el
número de iteración, y $`Z^{(k)}`$ es "la $Z$ de la iteración $k$". Un
superíndice **sin paréntesis** es una potencia: $`A^k`$ es $A$
multiplicada por sí misma $k$ veces.)

¿Qué le pasa a cada columna? Como en la potencia, expandimos en
eigenvectores, ahora con coeficientes de dos índices: $i$ dice de qué
eigenvector se trata y $j$ de qué columna:

```math
\mathbf{z}_j^{(0)} = \sum_{i=0}^{n-1} c_{ij}\mathbf{v}_i
\quad\Longrightarrow\quad
\mathbf{z}_j^{(k)} = \sum_{i=0}^{n-1} c_{ij}\lambda_i^k\mathbf{v}_i
= \lambda_0^k\left[c_{0j}\mathbf{v}_0 + \sum_{i=1}^{n-1} c_{ij}\left(\frac{\lambda_i}{\lambda_0}\right)^k\mathbf{v}_i\right] .
```

Uno esperaría que cada columna se quedara con "su" eigenvector. No
pasa: si $`|\lambda_0| > |\lambda_i|`$, los cocientes
$`(\lambda_i/\lambda_0)^k`$ se van a cero, y en **todas** las columnas
(las que tengan $`c_{0j} \neq 0`$) domina el término de
$\mathbf{v}_0$. Cada columna está haciendo, por su cuenta, el método de
la potencia, y todas llegan al mismo lugar.

**Un ejemplo a mano.** Tomemos

```math
A = \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix},
\qquad
\lambda_0 = 3,\ \ \mathbf{v}_0 = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 \\ 1 \end{pmatrix},
\qquad
\lambda_1 = 1,\ \ \mathbf{v}_1 = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 \\ -1 \end{pmatrix} .
```

Como $`\mathbf{e}_0 = (\mathbf{v}_0 + \mathbf{v}_1)/\sqrt{2}`$ y
$`\mathbf{e}_1 = (\mathbf{v}_0 - \mathbf{v}_1)/\sqrt{2}`$, al multiplicar
$k$ veces por $A$ cada $\mathbf{v}_i$ se multiplica por $`\lambda_i^k`$:

```math
Z^{(k)} = A^k = \frac{1}{2}\begin{pmatrix} 3^k + 1 & 3^k - 1 \\ 3^k - 1 & 3^k + 1 \end{pmatrix} .
```

(Para $k = 1$ sale $A$, compruébenlo.) Con $k = 10$, $`3^{10} = 59049`$
y las columnas son $(29525, 29524)$ y $(29524, 29525)$: las dos casi
paralelas a $(1, 1)$, o sea a $\mathbf{v}_0$. La información sobre
$\mathbf{v}_1$ no ha desaparecido del todo: la resta de las dos columnas
es $(1, -1)$ para todo $k$. Pero es la diferencia de dos números
enormes y casi iguales. Con números de punto flotante, en cuanto
$`3^k`$ pasa de unos $`10^{16}`$ el redondeo se la come (cancelación
catastrófica, unidad 06): con $k = 35$ la resta ya da $(4, -4)$, y
desde $k = 37$ da $(0, 0)$.

Con la matriz de prueba pasa lo mismo. La parte 3a de [`qr.py`](qr.py)
calcula $`A^{10}`$ y normaliza sus columnas: las cuatro son
$`(0.4444, 0.4822, 0.5172, 0.5500)`$, el eigenvector dominante, repetido
cuatro veces. (Además está el problema del desbordamiento, como en la
potencia sin normalizar.)

La lección: no basta con normalizar cada columna. Eso arregla el
desbordamiento, pero no que todas apunten al mismo lado. Hay que
impedir que las columnas se alineen.

### Iteración simultánea: ortonormalizando

Con un solo vector bastaba normalizar. Con varios hay que, además,
**quitarle a cada columna lo que tenga de las anteriores**:
ortogonalizar. Eso es justo lo que hace Gram-Schmidt, o sea, una
descomposición QR: las columnas de la $Q$ son ortonormales, y las
normas y componentes que se quitaron quedan guardadas en la $R$.

Para no confundir estas matrices con las del método QR (siguiente
sección), las de la iteración simultánea llevan una **raya abajo**:
$`\underline{Q}^{(k)}`$ y $`\underline{R}^{(k)}`$. La raya es solo parte
del nombre, no una operación.

La **iteración simultánea** (o *iteración ortogonal*) empieza con
$`\underline{Q}^{(0)} = I`$ y, para $k = 1, 2, \ldots$, hace dos pasos:

1. **Multiplicar** por $A$ las columnas ortonormales que ya tenemos:
   $`Z^{(k)} = A\,\underline{Q}^{(k-1)}`$.
2. **Ortonormalizar** el resultado con una descomposición QR:
   $`Z^{(k)} = \underline{Q}^{(k)}\,\underline{R}^{(k)}`$.

En fórmulas:

```math
Z^{(k)} = A\,\underline{Q}^{(k-1)},
\qquad
Z^{(k)} = \underline{Q}^{(k)}\,\underline{R}^{(k)} .
```

Para la siguiente iteración solo se usa $`\underline{Q}^{(k)}`$;
$`\underline{R}^{(k)}`$ no hace falta para el algoritmo, pero la vamos a
necesitar para entender qué calcula.

Compárenlo con la potencia normalizada,
$`\mathbf{z}^{(k)} = A\mathbf{q}^{(k-1)}`$,
$`\mathbf{q}^{(k)} = \mathbf{z}^{(k)}/\|\mathbf{z}^{(k)}\|`$: es lo
mismo con matrices en vez de vectores, y "ortonormalizar" en vez de
"normalizar". De hecho, la descomposición QR de una matriz con una
sola columna **es** normalizar: $`\mathbf{z} = \mathbf{q}\,R_{00}`$ con
$`R_{00} = \|\mathbf{z}\|`$.

**Por qué funciona.** La primera columna de $`\underline{Q}^{(k)}`$ no
cambia con Gram-Schmidt más que por la normalización: hace el método de
la potencia normalizado y tiende a $\mathbf{v}_0$. A la segunda columna
Gram-Schmidt le quita su componente en la dirección de la primera (que
ya es casi $\mathbf{v}_0$), así que en ella ya no puede dominar
$\mathbf{v}_0$: domina lo que sigue, $\mathbf{v}_1$. A la tercera se le
quitan las componentes en las direcciones de las dos primeras, y
domina $\mathbf{v}_2$. Y así sucesivamente: cada columna hace la
potencia en el espacio que le dejan las anteriores.

**La relación con las potencias de A.** Juntemos las $`\underline{R}`$
de todas las iteraciones en un producto, la más nueva a la izquierda:
$`\underline{R}^{(k)}\underline{R}^{(k-1)}\cdots\underline{R}^{(2)}\underline{R}^{(1)}`$.
Se cumple

```math
A^k = \underline{Q}^{(k)}\;\underline{R}^{(k)}\underline{R}^{(k-1)}\cdots\underline{R}^{(2)}\underline{R}^{(1)} .
```

Veamos las primeras $k$ una por una.

**Para $k = 1$:** como $`\underline{Q}^{(0)} = I`$, $`Z^{(1)} = A`$, y el paso
2 dice $`A = \underline{Q}^{(1)}\underline{R}^{(1)}`$.

**Para $k = 2$:** usamos el caso $k = 1$ para la segunda $A$, luego el paso 1
y luego el paso 2:

```math
\begin{aligned}
A^2 &= A\,A \\
&= A\,\underline{Q}^{(1)}\underline{R}^{(1)} \\
&= Z^{(2)}\,\underline{R}^{(1)} \\
&= \underline{Q}^{(2)}\underline{R}^{(2)}\underline{R}^{(1)} .
\end{aligned}
```

**En general,** si ya vale para $k - 1$, lo mismo:

```math
\begin{aligned}
A^k &= A\,A^{k-1} \\
&= A\,\underline{Q}^{(k-1)}\;\underline{R}^{(k-1)}\cdots\underline{R}^{(1)} \\
&= Z^{(k)}\;\underline{R}^{(k-1)}\cdots\underline{R}^{(1)} \\
&= \underline{Q}^{(k)}\underline{R}^{(k)}\;\underline{R}^{(k-1)}\cdots\underline{R}^{(1)} .
\end{aligned}
```

El producto de matrices triangulares superiores es triangular superior
(problema 4.2 de [1]), así que
$`\underline{R}^{(k)}\cdots\underline{R}^{(1)}`$ lo es. Entonces la
fórmula dice que $`\underline{Q}^{(k)}`$ es la $Q$ de la
**descomposición QR de $`A^k`$**. Ortonormalizar en cada paso da lo
mismo que calcular $`A^k`$ y ortonormalizar sus columnas al final (como
en la potencia, donde normalizar en cada paso o solo al final da la
misma dirección), pero sin desbordamiento y sin la cancelación
catastrófica del primer intento.

**El ejemplo, ortonormalizado.** Con la matriz de $2\times 2$ de
arriba, $`\underline{Q}^{(k)}`$ son las columnas de $`A^k`$
ortonormalizadas. La primera es la primera columna de $`A^k`$
normalizada,

```math
\mathbf{q}_0 = \frac{1}{\sqrt{(3^k+1)^2 + (3^k-1)^2}}\begin{pmatrix} 3^k + 1 \\ 3^k - 1 \end{pmatrix}
\;\longrightarrow\; \frac{1}{\sqrt{2}}\begin{pmatrix} 1 \\ 1 \end{pmatrix} = \mathbf{v}_0 .
```

La segunda tiene que ser perpendicular a la primera y, en dos
dimensiones, eso deja una sola dirección (Gram-Schmidt escoge el signo):

```math
\mathbf{q}_1 = \frac{1}{\sqrt{(3^k+1)^2 + (3^k-1)^2}}\begin{pmatrix} -(3^k - 1) \\ 3^k + 1 \end{pmatrix}
\;\longrightarrow\; \frac{1}{\sqrt{2}}\begin{pmatrix} -1 \\ 1 \end{pmatrix} = -\mathbf{v}_1 .
```

La segunda columna encontró a $\mathbf{v}_1$ (el signo no importa: un
eigenvector multiplicado por $-1$ sigue siendo eigenvector). Lo que en
el primer intento era la resta de dos números enormes, aquí se hace en
cada iteración, cuando las columnas todavía no son casi paralelas.

**Los eigenvectores.** Si $A$ es simétrica, sus eigenvectores son
ortonormales y $`\underline{Q}^{(k)}`$ tiende (salvo signos) a la matriz
de eigenvectores, como en el ejemplo. Si no es simétrica (el caso de la
matriz de prueba), los eigenvectores no son ortogonales y no pueden ser
las columnas de una matriz ortogonal. Lo que pasa entonces es que
$`\underline{Q}^{(k)}`$ tiende al "factor ortogonal" $\tilde{Q}$ de la
matriz de eigenvectores $V$: el de su descomposición QR,
$V = \tilde{Q}U$, con $U$ triangular superior (problema 4.32).

**Los eigenvalores.** En la potencia, el eigenvalor salía del cociente
de Rayleigh $`[\mathbf{q}^{(k)}]^T A\mathbf{q}^{(k)}`$. La
generalización a $n$ vectores es

```math
A^{(k)} = [\underline{Q}^{(k)}]^T A\,\underline{Q}^{(k)} .
```

(El elemento $j$ de su diagonal es $`\mathbf{q}_j^T A\mathbf{q}_j`$, el
cociente de Rayleigh de la columna $j$.) Es una **transformación de
similitud ortogonal**: $`A^{(k)}`$ tiene los mismos eigenvalores que
$A$. Y cuando $`\underline{Q}^{(k)} \to \tilde{Q}`$:

```math
\begin{aligned}
A^{(k)} &\to \tilde{Q}^T A\tilde{Q} \\
&= \tilde{Q}^T A V U^{-1} \\
&= \tilde{Q}^T V \Lambda U^{-1} \\
&= \tilde{Q}^T \tilde{Q} U \Lambda U^{-1} \\
&= U \Lambda U^{-1} .
\end{aligned}
```

(Se usó $`\tilde{Q} = VU^{-1}`$, luego $AV = V\Lambda$, luego otra vez
$V = \tilde{Q}U$, y por último $`\tilde{Q}^T\tilde{Q} = I`$.) El
resultado es un producto de tres matrices triangulares superiores (la
inversa de una triangular superior también lo es), así que es
**triangular superior**, con $\Lambda$ en la diagonal: los
eigenvalores se leen de la diagonal de $`A^{(k)}`$. Esta es la
descomposición de Schur que buscábamos. Si $A$ es simétrica, $`A^{(k)}`$
tiende a una matriz **diagonal**. (En el ejemplo de $2\times 2$,
$`A^{(k)}`$ tiende a la matriz diagonal con $3$ y $1$.)

Resumiendo, la iteración simultánea da

```math
A^k = \underline{Q}^{(k)}\;\underline{R}^{(k)}\underline{R}^{(k-1)}\cdots\underline{R}^{(1)},
\qquad
A^{(k)} = [\underline{Q}^{(k)}]^T A\,\underline{Q}^{(k)} .
```

La primera tiene que ver con los eigenvectores y la segunda con los
eigenvalores. Como criterio de paro se puede usar un número fijo de
iteraciones o vigilar cuánto cambia $`A^{(k)}`$.

### El método QR: algoritmo

El método QR es sorprendentemente simple. Empezando con
$`A^{(0)} = A`$, para $k = 1, 2, \ldots$ hace dos pasos:

1. **Descomponer** la matriz de la iteración anterior:
   $`A^{(k-1)} = Q^{(k)} R^{(k)}`$.
2. **Multiplicar en orden inverso**: $`A^{(k)} = R^{(k)} Q^{(k)}`$.

```math
A^{(k-1)} = Q^{(k)} R^{(k)},
\qquad
A^{(k)} = R^{(k)} Q^{(k)} .
```

Y otra vez, y otra vez. Ojo: estas $`Q^{(k)}`$ y $`R^{(k)}`$, **sin
raya**, no son las de la iteración simultánea. Aquí se descompone
$`A^{(k-1)}`$, allá se descomponía $`Z^{(k)}`$. Vamos a ver que están
relacionadas, y que al final los dos métodos calculan la misma
$`A^{(k)}`$.

**Primero: no hay magia.** Como $`Q^{(k)}`$ es ortogonal,
$`[Q^{(k)}]^TQ^{(k)} = I`$. Multiplicando el paso 1 por $`[Q^{(k)}]^T`$
por la izquierda, $`R^{(k)} = [Q^{(k)}]^T A^{(k-1)}`$, y sustituyendo
en el paso 2:

```math
A^{(k)} = [Q^{(k)}]^T A^{(k-1)}\, Q^{(k)} .
```

Cada paso del método QR es una **transformación de similitud
ortogonal**: $`A^{(k)}`$ tiene los mismos eigenvalores que
$`A^{(k-1)}`$, que tiene los mismos que $`A^{(k-2)}`$, y así hasta $A$.

**Las primeras iteraciones, una por una.** Usando lo anterior y que la
transpuesta de un producto es el producto de las transpuestas en orden
inverso, $`(BC)^T = C^TB^T`$:

**Para $k = 1$:** $`A^{(1)} = [Q^{(1)}]^T A\, Q^{(1)}`$.

**Para $k = 2$:** sustituimos $`A^{(1)}`$,

```math
\begin{aligned}
A^{(2)} &= [Q^{(2)}]^T A^{(1)}\, Q^{(2)} \\
&= [Q^{(2)}]^T [Q^{(1)}]^T A\, Q^{(1)} Q^{(2)} \\
&= [Q^{(1)} Q^{(2)}]^T A\, Q^{(1)} Q^{(2)} .
\end{aligned}
```

**Para $k = 3$:** igual,

```math
\begin{aligned}
A^{(3)} &= [Q^{(3)}]^T A^{(2)}\, Q^{(3)} \\
&= [Q^{(1)} Q^{(2)} Q^{(3)}]^T A\, Q^{(1)} Q^{(2)} Q^{(3)} .
\end{aligned}
```

En cada iteración, la $Q$ nueva se agrega **a la derecha** del
producto. En general:

```math
A^{(k)} = [Q^{(1)} Q^{(2)} \cdots Q^{(k-1)} Q^{(k)}]^T\, A\; Q^{(1)} Q^{(2)} \cdots Q^{(k-1)} Q^{(k)} .
```

El producto de matrices ortogonales es ortogonal (problema 4.2), así
que $`Q^{(1)}Q^{(2)}\cdots Q^{(k)}`$ lo es. Comparen con la iteración
simultánea, $`A^{(k)} = [\underline{Q}^{(k)}]^T A\,\underline{Q}^{(k)}`$:
es la misma fórmula si $`\underline{Q}^{(k)} = Q^{(1)}Q^{(2)}\cdots Q^{(k)}`$.
Eso es lo que vamos a probar.

**Los dos métodos son el mismo.** Queremos ver que, para todo $k$,

```math
\underline{Q}^{(k)} = Q^{(1)} Q^{(2)} \cdots Q^{(k-1)} Q^{(k)},
\qquad
\underline{R}^{(k)} = R^{(k)} .
```

La herramienta es que la descomposición QR es **única** (con la
diagonal de $R$ positiva, que es lo que da Gram-Schmidt): si logramos
escribir $`Z^{(k)} = PR`$ con $P$ ortogonal y $R$ triangular superior
con diagonal positiva, entonces forzosamente $`\underline{Q}^{(k)} = P`$
y $`\underline{R}^{(k)} = R`$.

**Para $k = 1$:** $`Z^{(1)} = A\,\underline{Q}^{(0)} = A`$ y $`A^{(0)} = A`$.
Los dos métodos descomponen la misma matriz, así que
$`\underline{Q}^{(1)} = Q^{(1)}`$ y $`\underline{R}^{(1)} = R^{(1)}`$.

**Para $k = 2$:** $`Z^{(2)} = A\,\underline{Q}^{(1)} = A\,Q^{(1)}`$. Ahora el
truco: escribimos la $A$ como $`Q^{(1)}R^{(1)}`$ y reconocemos
$`R^{(1)}Q^{(1)} = A^{(1)}`$, que el método QR descompone como
$`Q^{(2)}R^{(2)}`$:

```math
\begin{aligned}
Z^{(2)} &= A\,Q^{(1)} \\
&= Q^{(1)} R^{(1)} Q^{(1)} \\
&= Q^{(1)} A^{(1)} \\
&= Q^{(1)} Q^{(2)} R^{(2)} .
\end{aligned}
```

$`Q^{(1)}Q^{(2)}`$ es ortogonal y $`R^{(2)}`$ es triangular superior:
es la descomposición QR de $`Z^{(2)}`$. Así que
$`\underline{Q}^{(2)} = Q^{(1)}Q^{(2)}`$ y $`\underline{R}^{(2)} = R^{(2)}`$.

**En general,** si ya vale para $k - 1$, entonces
$`Z^{(k)} = A\,Q^{(1)}\cdots Q^{(k-1)}`$. Multiplicando la fórmula
de $`A^{(k-1)}`$ por $`Q^{(1)}\cdots Q^{(k-1)}`$ por la izquierda
(con una matriz ortogonal $P$, $`PP^T = I`$) queda
$`A\,Q^{(1)}\cdots Q^{(k-1)} = Q^{(1)}\cdots Q^{(k-1)}\,A^{(k-1)}`$, y
con el paso 1 del método QR:

```math
\begin{aligned}
Z^{(k)} &= A\,Q^{(1)} \cdots Q^{(k-1)} \\
&= Q^{(1)} \cdots Q^{(k-1)}\, A^{(k-1)} \\
&= Q^{(1)} \cdots Q^{(k-1)}\, Q^{(k)} R^{(k)} ,
\end{aligned}
```

que otra vez es una descomposición QR de $`Z^{(k)}`$.

Con esto, la relación con las potencias de $A$ de la iteración
simultánea se escribe con las matrices del método QR:

```math
A^k = Q^{(1)} Q^{(2)} \cdots Q^{(k-1)} Q^{(k)}\; R^{(k)} R^{(k-1)} \cdots R^{(2)} R^{(1)} .
```

Fíjense en el orden: las $Q$ van de la más vieja a la más nueva y las
$R$ al revés, y se encuentran en el centro, en $`Q^{(k)}R^{(k)}`$. Para
$k = 2$ se puede ver directamente:
$`A^2 = Q^{(1)}R^{(1)}\,Q^{(1)}R^{(1)} = Q^{(1)}A^{(1)}R^{(1)} = Q^{(1)}Q^{(2)}R^{(2)}R^{(1)}`$.

**¿Por qué preferir el método QR?** Los dos calculan lo mismo, pero
cada uno da directo una cosa distinta:

| | Iteración simultánea | Método QR |
|---|---|---|
| Se descompone | $`Z^{(k)} = A\,\underline{Q}^{(k-1)}`$ | $`A^{(k-1)}`$ |
| Sale directo | $`\underline{Q}^{(k)}`$ (eigenvectores) | $`A^{(k)} = R^{(k)}Q^{(k)}`$ (eigenvalores) |
| Lo otro cuesta | $`A^{(k)} = [\underline{Q}^{(k)}]^TA\,\underline{Q}^{(k)}`$: dos productos más | ir acumulando $`Q^{(1)}Q^{(2)}\cdots Q^{(k)}`$ (problema 4.35) |

Si solo queremos los eigenvalores, el método QR hace un producto de
matrices por iteración y la iteración simultánea tres
($`A\,\underline{Q}^{(k-1)}`$ y los dos de $`A^{(k)}`$). Y ya no hay
magia: el segundo paso, $`A^{(k)} = R^{(k)}Q^{(k)}`$, **es** una
transformación de similitud ortogonal.

Tal como está, el método QR es poco eficiente. En la práctica se
combina con desplazamientos (como en la potencia inversa) y con otros
trucos.

La parte 3b de [`qr.py`](qr.py) corre los dos métodos lado a lado y
compara las diagonales de $`A^{(k)}`$: coinciden hasta
$`\sim 10^{-13}`$. (Los elementos de arriba de la diagonal difieren
más, $`\sim 10^{-4}`$: la iteración simultánea usa
$`[\underline{Q}^{(k)}]^T`$ como si fuera la inversa de
$`\underline{Q}^{(k)}`$, y con Gram-Schmidt clásico la $Q$ no es del
todo ortogonal.)

### Implementación

Con `descomposicion_qr` y el producto de matrices, el método QR es
casi trivial:

```python
def metodo_qr(A, kmax=100):
    A_k = A
    for k in range(1, kmax):
        Q, R = descomposicion_qr(A_k)
        A_k = R * Q
    return [A_k[i][i] for i in range(A.rows)]
```

Se hace un número fijo de iteraciones, y por eso no
hace falta guardar los superíndices: `A_k`, `Q` y `R` se van
reemplazando. `A_k = A` no copia la matriz, pero no hace falta: nunca
se cambia un elemento de `A_k`, solo se reemplaza `A_k` completa por la
`Matrix` nueva que regresa `R * Q`, y la `A` de quien llama a la
función queda intacta.

La parte 4 de [`qr.py`](qr.py) lo corre paso a paso con la matriz de
prueba. La diagonal converge a los cuatro eigenvalores que encontramos
con la potencia inversa ($21.3166626635$, $-0.0672337148$,
$`-6.93\times 10^{-5}`$, $`-7.50\times 10^{-8}`$), y la parte de abajo
de la diagonal se va a cero, unas 300 veces más chica en cada
iteración: en 8 iteraciones pasa de $0.1$ a $`3\times 10^{-19}`$. El
factor 300 no es casualidad: el elemento $`A^{(k)}_{i+1,i}`$ decrece
como $`|\lambda_{i+1}/\lambda_i|^k`$, y el más lento de esos cocientes
es $`|\lambda_1/\lambda_0| = 3.15\times 10^{-3}`$, el mismo que
controlaba la potencia.

Con la matriz simétrica $`S = A + A^T`$, la parte de abajo también se
va a cero, pero la de arriba se queda en $`\sim 10^{-7}`$, aunque en
aritmética exacta $`A^{(k)}`$ seguiría siendo simétrica (y por lo tanto
diagonal). Otra vez es culpa de Gram-Schmidt clásico: si $Q$ no es
ortogonal, $RQ$ no es exactamente una transformación de similitud
ortogonal. Los eigenvalores, en la diagonal, sí salen bien (coinciden
con `numpy.linalg.eigvalsh`).

### Conteo de operaciones

En la descomposición QR, lo que más se repite es el ciclo interno:
para cada par $i < j$, el producto punto $`\mathbf{q}_i^T\mathbf{a}_j`$
($n$ multiplicaciones y $n-1$ sumas) y la resta
$`\mathbf{a}'_j - R_{ij}\mathbf{q}_i`$ ($n$ multiplicaciones y $n$
restas): $\sim 4n$ operaciones. Sumando sobre todos los pares:

```math
\sum_{j=0}^{n-1}\sum_{i=0}^{j-1} 4n = \sum_{j=0}^{n-1} 4nj = 4n\,\frac{(n-1)n}{2} \sim 2n^3 .
```

Una descomposición QR cuesta $\sim 2n^3$, **tres veces** más que LU
($\sim 2n^3/3$). El problema 4.34 pide la cuenta completa.

Cada iteración del método QR hace una descomposición ($\sim 2n^3$) y
un producto de matrices ($\sim 2n^3$). Con $m$ iteraciones, el total es
$\sim 4mn^3$. ¿Cuánto vale $m$? Depende de qué tan separados estén los
eigenvalores, pero no es raro que $m$ sea del orden de $n$, y entonces
el método QR es $`O(n^4)`$: demasiado para matrices grandes. Lo que se
hace en la práctica es llevar primero $A$ a la **forma de Hessenberg**
(ceros abajo de la subdiagonal, con reflexiones de Householder); eso
baja el costo de cada iteración de $`O(n^3)`$ a $`O(n^2)`$.

## Todos los eigenvalores y eigenvectores

`metodo_qr` da todos los eigenvalores. Para tener también los
eigenvectores podríamos acumular $`Q^{(1)}Q^{(2)}\cdots Q^{(k)}`$, pero para una
matriz no simétrica eso da el factor ortogonal $\tilde{Q}$, no los
eigenvectores. Mejor usamos lo que ya tenemos: con una buena
estimación del eigenvalor, la **potencia inversa con desplazamiento**
converge en muy pocas iteraciones y da el eigenvector. Esto era lo que
anunciamos en la sección de la potencia inversa:

```python
def eigenpares(A, eps=1e-12):
    eigenvalores = []
    eigenvectores = []
    for lam_qr in metodo_qr(A):
        lam, v, _ = potencia_inversa(A, desplazamiento=lam_qr + eps)
        eigenvalores.append(lam)
        eigenvectores.append(v)
    return eigenvalores, eigenvectores
```

Detalles:

- Al desplazamiento se le suma `eps`. Si el desplazamiento fuera
  *exactamente* un eigenvalor, $A - sI$ sería singular y LU
  dividiría entre cero; con el valor que da QR, que es casi exacto, el
  sistema queda muy mal condicionado. Lo que queremos es estar muy
  cerca pero no encima. (Recuerden que estar cerca **no** es problema
  para la potencia inversa: el error cae en la dirección del
  eigenvector.)
- Los eigenvalores que regresa son los de la potencia inversa
  (cociente de Rayleigh), no los de QR: QR solo da el punto de partida.
- Regresa una **lista de eigenvectores** (`eigenvectores[i]` es el de
  `eigenvalores[i]`), porque nuestros vectores son listas. Con NumPy
  sería más natural regresarlos como columnas de una matriz.
- `eigenpares` es tan corta porque usa todo lo anterior:
  `metodo_qr` usa `descomposicion_qr`, y `potencia_inversa` usa
  `descomposicion_lu` y las sustituciones.

La parte 5 de [`qr.py`](qr.py) la prueba. Los residuos
$`\|A\mathbf{v} - \lambda\mathbf{v}\|/\|\mathbf{v}\|`$ son todos
$`\sim 10^{-15}`$, la suma de los eigenvalores da la traza y su producto
el determinante, y el primer eigenvector es el de la potencia. Como se
reporta en [1], la potencia inversa cambia un poco el eigenvalor más
chico: QR da $`-7.501283325\times 10^{-8}`$ y la potencia inversa
$`-7.501283353\times 10^{-8}`$. (¿Cuál es mejor? Para un eigenvalor tan
chico comparado con $`\|A\| \approx 21`$, la diferencia es del orden
del redondeo.)

**Para pensar.**

- Programen el **Gram-Schmidt modificado** (problema 4.25): en cuanto
  tienen $\mathbf{q}_j$, réstenle su componente a **todas** las
  columnas que faltan, $`\mathbf{a}'_k \leftarrow \mathbf{a}'_k - (\mathbf{q}_j^T\mathbf{a}'_k)\,\mathbf{q}_j`$
  para $k > j$. En aritmética exacta es lo mismo; comparen
  $`\|Q^TQ - I\|`$ para $n = 4, 6, 8$. ¿Se arregla lo de la matriz
  simétrica?
- Programen la iteración simultánea como función (problema 4.33),
  parecida a `potencia`. Prueben primero la versión sin ortonormalizar
  y vean cómo todas las columnas acaban en $\mathbf{v}_0$.
- Modifiquen `metodo_qr` para que también acumule el producto
  $`Q^{(1)}Q^{(2)}\cdots Q^{(k)}`$ (problema 4.35). Con la matriz simétrica $S$,
  comparen sus columnas con los eigenvectores de `eigenpares`.
- Agreguen a `metodo_qr` un criterio de paro: por ejemplo, que la
  parte de abajo de la diagonal sea menor que una tolerancia.
- Prueben `eigenpares` con `eps = 0`.

## Descomposición en valores singulares (SVD)

Ya llevamos varias factorizaciones de una matriz: la **LU** (para
resolver sistemas lineales), la **QR** (para ortonormalizar, y como
paso del método QR) y la **eigendescomposición**
$`A = V\Lambda V^{-1}`$. Falta una quinta, igual de importante: la
**descomposición en valores singulares** (SVD, por *singular-value
decomposition*). En los cursos de álgebra lineal suele aparecer como un
resultado formal, y para calcularla se usa una biblioteca. Aquí, como
con todo lo demás, vamos a ver qué significa, probar que existe y
programarla.

### Definición

La SVD factoriza **cualquier** matriz $A$ como

```math
A = U\Sigma V^T ,
```

donde:

- $U$ es una matriz **ortogonal** ($`U^TU = I`$). Sus columnas
  $`\mathbf{u}_i`$ son los **vectores singulares izquierdos**.
- $\Sigma$ es **diagonal**, con los **valores singulares**
  $`\sigma_i \ge 0`$ en la diagonal, ordenados de mayor a menor:
  $`\sigma_0 \ge \sigma_1 \ge \cdots \ge \sigma_{n-1} \ge 0`$.
- $V$ también es **ortogonal**. Sus columnas $`\mathbf{v}_i`$ son los
  **vectores singulares derechos**.

Con ese orden, los valores singulares son únicos. Aquí nos quedamos
con matrices cuadradas de $n\times n$ (y entonces
$U$, $\Sigma$ y $V$ también lo son), pero toda matriz, cuadrada o no,
tiene una SVD (problema 4.45 de [1]).

Multiplicando $`A = U\Sigma V^T`$ por $V$ por la derecha queda
$`AV = U\Sigma`$, que columna por columna es

```math
A\mathbf{v}_i = \sigma_i\mathbf{u}_i .
```

Se parece al problema de eigenvalores, pero los vectores de la
izquierda y de la derecha son **distintos**.

**Toda matriz es diagonal en las bases adecuadas.** Tomemos el sistema
$`A\mathbf{x} = \mathbf{b}`$ de la unidad 09. La eliminación gaussiana
lo volvía **triangular**; la SVD lo vuelve **diagonal**. Sustituyendo
la SVD y multiplicando por $`U^T`$ por la izquierda:

```math
U^TU\Sigma V^T\mathbf{x} = U^T\mathbf{b}
\quad\Longrightarrow\quad
\Sigma\mathbf{x}' = \mathbf{b}',
\qquad
\mathbf{x}' = V^T\mathbf{x},\quad \mathbf{b}' = U^T\mathbf{b} .
```

$`\mathbf{x}'`$ es $\mathbf{x}$ escrito en la base de los vectores
singulares derechos (las columnas de $V$), y $`\mathbf{b}'`$ es
$\mathbf{b}$ en la base de los izquierdos (las columnas de $U$). En
esas coordenadas cada ecuación es independiente de las demás:
$`\sigma_i x'_i = b'_i`$.

**¿No es lo mismo que diagonalizar?** En esta sección llamamos $W$ a
la matriz de eigenvectores (la $V$ de antes, porque ahora $V$ es otra
cosa): $`A = W\Lambda W^{-1}`$. Con $`\mathbf{x}' = W^{-1}\mathbf{x}`$ y
$`\mathbf{b}' = W^{-1}\mathbf{b}`$ el sistema también queda diagonal,
$`\Lambda\mathbf{x}' = \mathbf{b}'`$. Pero hay diferencias de fondo:

- La SVD existe para **cualquier** matriz. La eigendescomposición solo
  para las diagonalizables (no para una matriz defectuosa, por
  ejemplo, ni para una que no sea cuadrada).
- La SVD usa **dos** bases, $U$ y $V$; la eigendescomposición, una
  sola ($W$).
- Las bases de la SVD son **ortonormales**. Los eigenvectores de una
  matriz no simétrica son linealmente independientes, pero no
  ortogonales.
- Los valores singulares son **reales y no negativos**; los
  eigenvalores pueden ser negativos, o complejos.

La SVD tiene muchísimas aplicaciones. La más cercana a lo que ya vimos:
la **norma 2** de una matriz es $`\sigma_0`$, y el número de condición
en esa norma es $`\kappa_2(A) = \sigma_0/\sigma_{n-1}`$ (lo vemos abajo).
Otras: trabajar con matrices no invertibles (resolver sistemas
homogéneos $`A\mathbf{x} = \mathbf{0}`$), **aproximar** una matriz por
otra más sencilla que conserve lo esencial, y el ajuste por mínimos
cuadrados (capítulo 6 de [1]).

### La SVD existe

Antes de calcularla, probemos que toda matriz $A$ de $n\times n$ tiene
una SVD, sin suponer nada de $A$ (ni que sea simétrica, ni
diagonalizable). Es una prueba por inducción, como las del método QR.

**Caso base, $n = 1$.** $A$ es un número $a$. Si $a \neq 0$, basta
$`U = a/|a|`$ (que es $\pm 1$), $`\Sigma = |a|`$ y $V = 1$. (Si $a = 0$,
$U = V = 1$ y $\Sigma = 0$.)

**La norma 2.** Necesitamos una norma de matrices nueva, la **norma 2**
o **norma espectral**, definida a partir de la norma euclídea de
vectores:

```math
\|A\|_2 = \max_{\mathbf{x}\neq\mathbf{0}} \frac{\|A\mathbf{x}\|}{\|\mathbf{x}\|}
= \max_{\|\mathbf{x}\| = 1} \|A\mathbf{x}\| .
```

Es lo más que $A$ puede "estirar" a un vector. (Las dos formas son
iguales porque $`\|A(c\mathbf{x})\|/\|c\mathbf{x}\|`$ no depende de $c$.)
Una propiedad que vamos a usar: las matrices ortogonales **no cambian**
la norma euclídea, porque
$`\|Q\mathbf{x}\|^2 = \mathbf{x}^TQ^TQ\mathbf{x} = \mathbf{x}^T\mathbf{x} = \|\mathbf{x}\|^2`$.
Por eso, si $`Q_1`$ y $`Q_2`$ son ortogonales,
$`\|Q_1AQ_2\|_2 = \|A\|_2`$ (problema 4.40).

**Preparación.** Sea $\mathbf{v}$ el vector de norma 1 en el que se
alcanza el máximo, y

```math
\sigma = \|A\|_2 = \|A\mathbf{v}\|,
\qquad
\mathbf{u} = \frac{A\mathbf{v}}{\|A\mathbf{v}\|} = \frac{A\mathbf{v}}{\sigma}
```

(si $A = 0$ no hay nada que probar, así que $\sigma > 0$). Ahora
completamos $\mathbf{v}$ a una base ortonormal (por ejemplo con
Gram-Schmidt): $`V = (\mathbf{v}\ \ \tilde{V})`$ es ortogonal de
$n\times n$, y $`\tilde{V}`$ tiene las otras $n-1$ columnas. Igual con
$\mathbf{u}$: $`U = (\mathbf{u}\ \ \tilde{U})`$. Multiplicando por bloques,

```math
U^TAV =
\begin{pmatrix} \mathbf{u}^T \\ \tilde{U}^T \end{pmatrix}
A
\begin{pmatrix} \mathbf{v} & \tilde{V} \end{pmatrix}
=
\begin{pmatrix}
\mathbf{u}^TA\mathbf{v} & \mathbf{u}^TA\tilde{V} \\
\tilde{U}^TA\mathbf{v} & \tilde{U}^TA\tilde{V}
\end{pmatrix} .
```

Veamos cada bloque:

- Arriba a la izquierda (un número):
  $`\mathbf{u}^TA\mathbf{v} = (A\mathbf{v})^T(A\mathbf{v})/\sigma = \sigma^2/\sigma = \sigma`$.
- Abajo a la izquierda (una columna de $n-1$ elementos):
  $`\tilde{U}^TA\mathbf{v} = \sigma\,\tilde{U}^T\mathbf{u} = \mathbf{0}`$,
  porque las columnas de $`\tilde{U}`$ son ortogonales a $\mathbf{u}$.
- Arriba a la derecha (un renglón de $n-1$ elementos): le ponemos
  nombre, $`\mathbf{w}^T = \mathbf{u}^TA\tilde{V}`$.
- Abajo a la derecha (una matriz de $`(n-1)\times(n-1)`$):
  $`\tilde{A} = \tilde{U}^TA\tilde{V}`$.

Así que

```math
U^TAV = \begin{pmatrix} \sigma & \mathbf{w}^T \\ \mathbf{0} & \tilde{A} \end{pmatrix} .
```

**El bloque de arriba a la derecha también es cero.** Multipliquemos
esa matriz por el vector que tiene $\sigma$ arriba y $\mathbf{w}$
abajo:

```math
\left\| \begin{pmatrix} \sigma & \mathbf{w}^T \\ \mathbf{0} & \tilde{A} \end{pmatrix}
\begin{pmatrix} \sigma \\ \mathbf{w} \end{pmatrix} \right\|
=
\left\| \begin{pmatrix} \sigma^2 + \|\mathbf{w}\|^2 \\ \tilde{A}\mathbf{w} \end{pmatrix} \right\|
=
\sqrt{\left(\sigma^2 + \|\mathbf{w}\|^2\right)^2 + \|\tilde{A}\mathbf{w}\|^2}
\ \ge\ \sigma^2 + \|\mathbf{w}\|^2 .
```

Entonces

```math
\sigma = \|A\|_2 = \|U^TAV\|_2
\ \ge\ \frac{\left\| U^TAV \begin{pmatrix} \sigma \\ \mathbf{w} \end{pmatrix}\right\|}
{\left\| \begin{pmatrix} \sigma \\ \mathbf{w} \end{pmatrix}\right\|}
\ \ge\ \frac{\sigma^2 + \|\mathbf{w}\|^2}{\sqrt{\sigma^2 + \|\mathbf{w}\|^2}}
= \sqrt{\sigma^2 + \|\mathbf{w}\|^2} .
```

En el segundo paso se usó que $U$ y $V$ son ortogonales; en el
tercero, que la norma 2 es un **máximo**, así que es mayor o igual que
el cociente para un vector en particular; en el cuarto, la desigualdad
anterior. Comparando el principio con el final,
$`\sigma \ge \sqrt{\sigma^2 + \|\mathbf{w}\|^2}`$, y eso solo puede
pasar si $`\mathbf{w} = \mathbf{0}`$.

**Paso inductivo.** Supongamos que la SVD existe para matrices de
$`(n-1)\times(n-1)`$; en particular,
$`\tilde{A} = U'\Sigma'V'^T`$. Entonces

```math
U^TAV
= \begin{pmatrix} \sigma & \mathbf{0}^T \\ \mathbf{0} & U'\Sigma'V'^T \end{pmatrix}
= \begin{pmatrix} 1 & \mathbf{0}^T \\ \mathbf{0} & U' \end{pmatrix}
\begin{pmatrix} \sigma & \mathbf{0}^T \\ \mathbf{0} & \Sigma' \end{pmatrix}
\begin{pmatrix} 1 & \mathbf{0}^T \\ \mathbf{0} & V' \end{pmatrix}^T
```

(hagan la multiplicación para convencerse). Multiplicando por $U$ por
la izquierda y por $`V^T`$ por la derecha:

```math
A =
\left[ U \begin{pmatrix} 1 & \mathbf{0}^T \\ \mathbf{0} & U' \end{pmatrix} \right]
\begin{pmatrix} \sigma & \mathbf{0}^T \\ \mathbf{0} & \Sigma' \end{pmatrix}
\left[ V \begin{pmatrix} 1 & \mathbf{0}^T \\ \mathbf{0} & V' \end{pmatrix} \right]^T .
```

Las matrices entre corchetes son productos de matrices ortogonales, y
por lo tanto ortogonales; la de en medio es diagonal y no negativa.
Esto es una SVD de $A$. (Y quedan bien ordenados:
$`\sigma'_0 = \|\tilde{A}\|_2 \le \|U^TAV\|_2 = \sigma`$, porque
$`\tilde{A}`$ es un pedazo de $`U^TAV`$ y no puede estirar más.)

De paso, la prueba dice que $`\sigma_0 = \|A\|_2`$, y que
$`\mathbf{v}_0`$ es el vector que $A$ más estira.

### Cómo calcular la SVD

Una propiedad muy conocida: los eigenvalores de la matriz simétrica
$`A^TA`$ son los $`\sigma_i^2`$, y sus eigenvectores son los
$`\mathbf{v}_i`$. (Pruébenlo: $`A^TA = V\Sigma U^TU\Sigma V^T = V\Sigma^2V^T`$.)
Muchos libros calculan así la SVD, pero se puede mostrar
(problemas 4.41 y 4.42 de [1]) que es **inestable**: el resultado es mucho más
sensible a perturbaciones que el problema original. Lo que se gana en
intuición se pierde en la práctica.

Seguimos otro camino, igual de sencillo y más estable, que usa la misma
idea que las bibliotecas profesionales (que además, para ser
eficientes, llevan primero la matriz a forma *bidiagonal*). Partimos
de $`AV = U\Sigma`$. Transponiendo $`A = U\Sigma V^T`$ y multiplicando
por $U$ por la derecha sale la otra mitad:

```math
AV = U\Sigma,
\qquad
A^TU = V\Sigma
```

(usando $`\Sigma^T = \Sigma`$). Las dos juntas son **una sola**
ecuación con bloques de $2\times 2$:

```math
\begin{pmatrix} 0 & A^T \\ A & 0 \end{pmatrix}
\frac{1}{\sqrt{2}}
\begin{pmatrix} V & V \\ U & -U \end{pmatrix}
=
\frac{1}{\sqrt{2}}
\begin{pmatrix} V & V \\ U & -U \end{pmatrix}
\begin{pmatrix} \Sigma & 0 \\ 0 & -\Sigma \end{pmatrix} .
```

(Hagan las multiplicaciones por bloques: el bloque de arriba a la
izquierda da $`A^TU = V\Sigma`$, el de abajo a la izquierda
$`AV = U\Sigma`$, y los de la derecha lo mismo con signos.) Esto es
la diagonalización $`HW = W\Lambda`$ de la matriz de $2n\times 2n$

```math
H = \begin{pmatrix} 0 & A^T \\ A & 0 \end{pmatrix},
\qquad
W = \frac{1}{\sqrt{2}} \begin{pmatrix} V & V \\ U & -U \end{pmatrix},
\qquad
\Lambda = \begin{pmatrix} \Sigma & 0 \\ 0 & -\Sigma \end{pmatrix} .
```

En palabras: los eigenvalores de $H$ son los $`\sigma_i`$ y los
$`-\sigma_i`$. El eigenvector de $`+\sigma_i`$ (la columna $i$ de $W$)
tiene arriba $`\mathbf{v}_i/\sqrt{2}`$ y abajo $`\mathbf{u}_i/\sqrt{2}`$;
el $`1/\sqrt{2}`$ es para que tenga norma 1. Calcular la SVD de $A$ se
convirtió en el problema de eigenvalores de $H$, sin pasar por
$`A^TA`$. Y $H$ es **simétrica**, el mejor caso para eigenvalores.

### Implementación

`matriz_h(A)` arma $H$: ceros, $A$ en el bloque
de abajo a la izquierda y $`A^T`$ en el de arriba a la derecha. La SVD
la hace `svd(A)`:

```python
def svd(A, desplazamiento=1.0, eps=1e-12):
    n = A.rows
    H = matriz_h(A)
    H_desplazada = H - Identity(2 * n) * desplazamiento
    eigenvalores_h = [lam + desplazamiento for lam in metodo_qr(H_desplazada)]
    sigmas_qr = sorted([lam for lam in eigenvalores_h if lam >= 0], reverse=True)
    sigmas = []
    columnas_u = []
    columnas_v = []
    for sigma_qr in sigmas_qr:
        sigma, w, _ = potencia_inversa(H, desplazamiento=sigma_qr + eps)
        sigmas.append(sigma)
        columnas_v.append([raiz_cuadrada(2.0) * w_i for w_i in w[:n]])
        columnas_u.append([raiz_cuadrada(2.0) * w_i for w_i in w[n:]])
    return Matrix(columnas_u).transpose(), sigmas, Matrix(columnas_v)
```

Conceptualmente es eso: eigenvalores de $H$ y luego eigenvectores. Los
detalles:

- **El desplazamiento.** Los eigenvalores de $H$ vienen en pares
  $`\pm\sigma_i`$, con el **mismo valor absoluto**, y el método QR (como
  lo programamos) solo funciona con eigenvalores de valores absolutos
  distintos. El remedio más simple es el de la potencia inversa:
  aplicar QR a $`H - sI`$, cuyos eigenvalores $`\sigma_i - s`$ y
  $`-\sigma_i - s`$ ya tienen valores absolutos distintos, y al final
  volver a sumar $s$. Usamos $s = 1$.
- **Solo los no negativos.** De los $2n$ eigenvalores de $H$ nos
  quedamos con los $`\ge 0`$, que son los $`\sigma_i`$, y los ordenamos
  de mayor a menor.
- **Los vectores singulares.** Como en `eigenpares`, cada $`\sigma_i`$
  (más `eps`) es el desplazamiento de la potencia inversa, que da el
  eigenvector $\mathbf{w}$ de $H$ y refina $`\sigma_i`$. Sus primeros
  $n$ elementos, por $\sqrt{2}$, son $`\mathbf{v}_i`$; los últimos $n$,
  por $\sqrt{2}$, son $`\mathbf{u}_i`$. El signo de $\mathbf{w}$ no
  importa: si cambia, cambian $`\mathbf{u}_i`$ y $`\mathbf{v}_i`$ a la
  vez, y $`\sigma_i\mathbf{u}_i\mathbf{v}_i^T`$ sigue igual.
- **Qué regresa.** $U$ como matriz, los $`\sigma_i`$ como **lista** y
  $`V^T`$ (no $V$), como hacen las bibliotecas. Como `columnas_v` es
  la lista de las columnas de $V$, que son los renglones de $`V^T`$,
  ya es $`V^T`$ sin transponer.

### Resultados

[`svd.py`](svd.py) usa la matriz de prueba **menos la
identidad**, `crear_prueba(4, 21)` $- I$ (abajo vemos por qué).

**La parte 1** muestra el desplazamiento en acción. Con $H$ tal cual,
la diagonal de $`H^{(k)}`$ se queda en **ceros**: QR no puede separar
$\sigma$ de $-\sigma$. Con $s = 1$, los cuatro eigenvalores positivos
salen bien. Los negativos convergen mal (los tres más chicos quedan
en $`-\sigma_i - 1 \approx -2`$, casi empatados en valor absoluto),
pero esos se descartan. Incluso $`\sigma_0`$ sale de QR con un error de
$`\sim 10^{-7}`$, porque $`|\sigma_0 - 1| = 19.36`$ y
$`|-\sigma_0 - 1| = 21.36`$ están cerca (cociente $0.9$); la potencia
inversa lo corrige.

**La parte 2** calcula la SVD. Los valores singulares son
$20.3559952277$, $1.0652108300$, $1.0000326288$ y $0.9999998880$
(coinciden con `numpy.linalg.svd`). Se comprueba de varias formas:

- $`A\mathbf{v}_i = \sigma_i\mathbf{u}_i`$, con residuos $`\sim 10^{-14}`$.
- $`\sum_i \sigma_i^2 = \|A\|_F^2`$ y $`\prod_i \sigma_i = |\det A|`$, a
  10 cifras. (Las dos salen de que $U$ y $V$ son ortogonales; para la
  segunda, $`|\det U| = |\det V| = 1`$.)
- $`\|A - U\Sigma V^T\|`$, $`\|U^TU - I\|`$ y $`\|V^TV - I\|`$ son
  $`\sim 2\times 10^{-10}`$, no $`10^{-16}`$. No es un error del
  programa: $`\sigma_2`$ y $`\sigma_3`$ están a solo $`3\times 10^{-5}`$
  uno del otro, y como vimos en "Sensibilidad de los eigenvectores",
  un eigenvector con un vecino tan cercano es sensible: el redondeo lo
  puede girar hacia el otro unos
  $`\epsilon\|H\|/(3\times 10^{-5}) \approx 10^{-10}`$. Por eso cada
  par por separado cumple $`A\mathbf{v}_i = \sigma_i\mathbf{u}_i`$ muy
  bien, pero $`\mathbf{v}_2`$ y $`\mathbf{v}_3`$ ya no son del todo
  ortogonales entre sí.

**La parte 3** comprueba que $`\|A\|_2 = \sigma_0`$: con 10000 vectores
al azar, el mayor $`\|A\mathbf{x}\|/\|\mathbf{x}\|`$ es $20.3486$, un
poco menos que $`\sigma_0 = 20.3560`$, que se alcanza justo en
$`\mathbf{v}_0`$. De la misma forma, $`\|A^{-1}\|_2 = 1/\sigma_{n-1}`$
(los valores singulares de $`A^{-1} = V\Sigma^{-1}U^T`$ son los
$`1/\sigma_i`$), así que el número de condición en norma 2 es

```math
\kappa_2(A) = \|A\|_2\,\|A^{-1}\|_2 = \frac{\sigma_0}{\sigma_{n-1}} .
```

Para esta matriz, $`\kappa_2 \approx 20.4`$: restarle $I$ la volvió
bien condicionada.

**La parte 4** resuelve $`A\mathbf{x} = \mathbf{b}`$ con la SVD:
$`\mathbf{b}' = U^T\mathbf{b}`$, luego $`x'_i = b'_i/\sigma_i`$ y
$`\mathbf{x} = V\mathbf{x}'`$. Da lo mismo que LU. (Es más caro que
LU, pero si algún $`\sigma_i`$ es cero o casi, se ve de inmediato qué
componente de $`\mathbf{x}'`$ está mal determinada.)

**La parte 5** muestra el límite de nuestra `svd`. Sin restarle $I$,
los valores singulares de la matriz de prueba van de $21.35$ a
$`7.5\times 10^{-8}`$. Con $s = 1$, los dos más chicos dan
$`\sigma_i - 1 \approx -1`$ y sus compañeros negativos
$`-\sigma_i - 1 \approx -1`$: valores absolutos que difieren en
$`\sim 10^{-7}`$, y el método QR no los separa en 100 iteraciones. Se
nota en todo: $`\sigma_2`$ sale mal, $`\prod_i\sigma_i`$ no da el
determinante, $`\mathbf{v}_2\cdot\mathbf{v}_3 \approx -0.8`$ (la
potencia inversa encontró casi el mismo vector dos veces) y
$`\|A - U\Sigma V^T\| \sim 10^{-4}`$. Esto pasa en general cuando los
valores singulares abarcan muchos órdenes de magnitud, y **es culpa del
método para los eigenvalores**, no de la idea de usar $H$: con un método
para eigenvalores más robusto (por ejemplo `numpy.linalg.eigvals`),
la SVD sale correcta también para estas matrices [1].

**Para pensar.**

- Calculen la SVD por el camino de $`A^TA`$ (problema 4.41): los
  eigenvalores de $`A^TA`$ (con `eigenpares`) son los $`\sigma_i^2`$, y
  sus eigenvectores los $`\mathbf{v}_i`$; los $`\mathbf{u}_i`$ salen de
  $`A\mathbf{v}_i = \sigma_i\mathbf{u}_i`$. Comparen con `svd`. ¿Cuánto
  vale $`\kappa_2(A^TA)`$ comparado con $`\kappa_2(A)`$? (Problema 4.42.)
- ¿Hay algún desplazamiento $s$ con el que `svd` funcione para la
  matriz de la parte 5? Piensen en los valores absolutos de
  $`\sigma_i - s`$ y $`-\sigma_i - s`$ cuando $`\sigma_i`$ es muy chico.
- **Aproximación de rango bajo.** La SVD se puede escribir como una
  suma, $`A = \sum_i \sigma_i\mathbf{u}_i\mathbf{v}_i^T`$. Quédense con
  los primeros $k$ términos, $`A_k`$, y comparen $`\|A - A_k\|_F`$ con
  $`\sqrt{\sigma_k^2 + \cdots + \sigma_{n-1}^2}`$.
- Si $A$ es singular, su último valor singular es cero, y
  $`\mathbf{v}_{n-1}`$ resuelve el sistema homogéneo
  $`A\mathbf{x} = \mathbf{0}`$. ¿Por qué? Prueben con una matriz cuyo
  último renglón sea la suma de los otros. ¿Lo logra nuestra `svd`?

## Contenido
- [`fiscomp/algebra_lineal.py`](../../fiscomp/algebra_lineal.py): los
  métodos de la potencia (`producto_punto`, `normalizar`,
  `cociente_rayleigh`, `potencia`, `potencia_inversa`) y el método QR
  (`descomposicion_qr`, `metodo_qr`, `eigenpares`) y la SVD (`matriz_h`, `svd`). Es el mismo módulo de la unidad 09: la potencia
  inversa usa `descomposicion_lu` y `resolver_lu`, y todos los
  ejemplos usan su matriz de prueba (`crear_prueba`).
- [`condicion_eigenvalores.py`](condicion_eigenvalores.py): el número
  de condición $\kappa^{ev}$ de los eigenvalores de una matriz no
  simétrica contra una simétrica, comparado con el cambio real al
  perturbar la matriz; y cómo gira un eigenvector cuando dos
  eigenvalores están cerca. Otra vez $2\times 2$, con el polinomio
  característico.
- [`potencia.py`](potencia.py): el método de la potencia iteración por
  iteración, con el error del eigenvalor bajando por un factor
  $|\lambda_1/\lambda_0|$; los cuatro eigenvalores de la matriz de
  prueba con la potencia inversa y distintos desplazamientos,
  comprobados con la traza y el determinante; cuántas iteraciones
  hacen falta según qué tan cerca está el desplazamiento; y el cambio
  de signo con un eigenvalor negativo.
- [`qr.py`](qr.py): la descomposición QR de la matriz de prueba
  ($`\|A - QR\|`$ y la pérdida de ortogonalidad de $Q$); por qué la
  diagonal de $U$ no da los eigenvalores; la iteración simultánea sin
  y con ortonormalizar, comparada paso a paso con el método QR; el
  método QR volviendo triangular a $`A^{(k)}`$ (y casi diagonal a una
  matriz simétrica); y todos los eigenpares con `eigenpares`,
  comprobados con el residuo, la traza y el determinante.
- [`svd.py`](svd.py): los eigenvalores de $H$ con y sin desplazamiento;
  la SVD de la matriz de prueba menos la identidad, comprobada con la
  reconstrucción $`U\Sigma V^T`$, la norma de Frobenius y el
  determinante; la norma 2 y $`\kappa_2`$; resolver
  $`A\mathbf{x} = \mathbf{b}`$ en las bases de la SVD; y un caso en el
  que falla (valores singulares de muchos órdenes de magnitud).

Se corren desde la raíz del repositorio (con el `.venv` activado), por
ejemplo `python3 unidades/10_valores_y_vectores_propios/qr.py`.

## Referencias

[1] A. Gezerlis, *Numerical Methods in Physics with Python*, vol. 1.
Cambridge University Press, Cambridge, UK, 2023.
