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
  eigenpares combinando QR con la potencia inversa.
- Corresponde a la sección 4.4 (hasta la 4.4.4) del libro de Gezerlis,
  y a la parte de eigenvalores de la 4.2.
- Se apoya en la [unidad 09](../09_algebra_lineal/) (sistemas de
  ecuaciones lineales): normas, número de condición, descomposición LU
  y el criterio de paro de Jacobi.

## Notación

Los índices van de $0$ a $n-1$, como en el libro y en Python. Las
matrices se escriben con mayúsculas ($A$, $L$, $U$) y los vectores en
negritas minúsculas ($\mathbf{x}$, $\mathbf{b}$). Cuando hace falta,
los dos índices de un elemento se separan con una coma: $A_{1,n-1}$.
Usamos "eigenvalor" y "valor propio" como sinónimos (igual
"eigenvector" y "vector propio"); "eigenpar" es un eigenvalor junto
con su eigenvector.

En el código, una matriz es una **lista de listas** (un renglón por
sublista, `A[i][j]`) y un vector es una lista simple (`x[i]`). Es el
mismo formato que guarda la clase `Matrix` de la práctica 2 en su
atributo `.data`. En este tema **no usamos NumPy**: todos los métodos
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

En la [unidad 09](../09_algebra_lineal/notas.md) (sección "Análisis de
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

(La prueba es el problema 4.21 del libro: escriban $AV$ columna por
columna.) Por eso a resolver el problema de eigenvalores también se le
llama **diagonalizar** la matriz, y a $A = V\Lambda V^{-1}$ se le
llama la **eigendescomposición**.

No vamos a suponer que $A$ es simétrica: muchos problemas de física
dan matrices simétricas, pero no todos. Eso sí, una matriz no
simétrica puede tener eigenvalores complejos; aquí solo veremos
matrices no simétricas con eigenvalores reales (y distintos). La
matriz de prueba del libro, `crear_prueba(4, 21)`, es de ese tipo.

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
    q = normalizar([1.0] * len(A))
    for k in range(1, kmax):
        z = mat_vec(A, q)
        q = normalizar(z)
    return cociente_rayleigh(A, q), q
```

Igual que en el libro, el vector inicial es $`(1, 1, \ldots, 1)`$
normalizado, se hace un número fijo de iteraciones, y el cociente de
Rayleigh se calcula **una sola vez**, al final: calcularlo en cada
iteración sería desperdiciar un producto matriz-vector por vuelta.

La primera parte de [`potencia.py`](potencia.py) lo prueba con la
matriz de prueba del libro y reproduce su resultado,
$\lambda_0 = 21.3166626635$ y
$`\mathbf{v}_0 = (0.44439562, 0.48218122, 0.51720906, 0.55000781)`$.
El libro compara contra `numpy.linalg.eig`; como aquí no usamos NumPy,
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
eigenvectores. (El proyecto de espines ya lo hacía, en
`iteracion_inversa`.)

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
    n = len(A)
    A_desplazada = [[A[i][j] - (desplazamiento if i == j else 0.0) for j in range(n)] for i in range(n)]
    L, U = descomposicion_lu(A_desplazada)
    q = normalizar([1.0] * n)
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
  listas distintas (`normalizar` crea una lista nueva); el libro tiene
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
  como en `potencia_inversa` (problema 4.23 del libro). Prueben con
  `crear_prueba(20, 21)`: el libro dice que bastan 6 iteraciones.
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
(problema 4.27 del libro).

### Implementación

Como nuestras matrices son listas de **renglones** y Gram-Schmidt
trabaja con **columnas**, primero agregamos a
[`fiscomp/algebra_lineal.py`](../../fiscomp/algebra_lineal.py) la
`transpuesta` (las columnas de $A$ son los renglones de $A^T$) y el
producto de matrices `mat_mat` (que también necesitamos para el
método QR). Con ellas:

```python
def descomposicion_qr(A):
    n = len(A)
    columnas = transpuesta(A)
    q = []  # columnas de Q, ya ortonormales
    R = [[0.0] * n for _ in range(n)]
    for j in range(n):
        a_prima = columnas[j][:]
        for i in range(j):
            R[i][j] = producto_punto(q[i], columnas[j])
            a_prima = [a - R[i][j] * q_i for a, q_i in zip(a_prima, q[i])]
        R[j][j] = norma_euclidea(a_prima)
        q.append([a / R[j][j] for a in a_prima])
    return transpuesta(q), R
```

Es una traducción directa de las fórmulas: `a_prima` empieza siendo
una copia de $\mathbf{a}_j$ y se le van restando las componentes; el
ciclo interno no se ejecuta para $j = 0$. Al final, `q` es la lista de
columnas de $Q$, y la transponemos para regresar $Q$ como lista de
renglones. No hay `kmax` ni tolerancia: es un método directo.

¿Cómo sabemos si funciona? Con dos pruebas. La obvia,
$`\|A - QR\|`$, que debe ser cero. Y otra: en aritmética exacta $Q$ es
ortogonal por construcción, pero con redondeo hay que comprobarlo con
$`\|Q^TQ - I\|`$. La parte 1 de [`qr.py`](qr.py) lo hace con la matriz
de prueba (norma de Frobenius):

| $n$ | $`\|A - QR\|`$ | $`\|Q^TQ - I\|`$ |
|---|---|---|
| 4 | $`1.3\times 10^{-15}`$ | $`7.3\times 10^{-3}`$ |
| 6 | $`2.5\times 10^{-15}`$ | $2.4$ |
| 8 | $`6.2\times 10^{-15}`$ | $4.5$ |

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

Generalicemos el método de la potencia a $n$ vectores a la vez, para
sacar los $n$ eigenvectores. En vez de un vector inicial
$`\mathbf{z}^{(0)}`$ usamos $n$, como columnas de una matriz
$`Z^{(0)}`$; la elección natural es la identidad, $`Z^{(0)} = I`$
(columnas $\mathbf{e}_0, \ldots, \mathbf{e}_{n-1}$, linealmente
independientes). Y multiplicamos por $A$ una y otra vez:

```math
Z^{(k)} = A Z^{(k-1)} = A^k Z^{(0)}, \qquad k = 1, 2, \ldots
```

Cada columna se expande en eigenvectores, ahora con coeficientes de
dos índices (uno para la suma y otro para saber de qué columna se
trata):

```math
\mathbf{z}_j^{(0)} = \sum_{i=0}^{n-1} c_{ij}\mathbf{v}_i
\quad\Longrightarrow\quad
\mathbf{z}_j^{(k)} = \sum_{i=0}^{n-1} c_{ij}\lambda_i^k\mathbf{v}_i .
```

Uno esperaría que cada columna se quedara con "su" eigenvector. No
pasa: por el mismo argumento que en la potencia, en **todas** las
columnas domina el término de $\lambda_0$, y todas tienden a
$\mathbf{v}_0$. La parte 3a de [`qr.py`](qr.py) calcula $`A^{10}`$ y
normaliza sus columnas: las cuatro son
$`(0.4444, 0.4822, 0.5172, 0.5500)`$, el eigenvector dominante, repetido
cuatro veces. (Además está el problema del desbordamiento, como en la
potencia sin normalizar.)

### Iteración simultánea: ortonormalizando

Con un solo vector bastaba normalizar. Con varios hay que, además,
**quitarle a cada columna lo que tenga de las otras**: ortogonalizar.
Eso es justo lo que hace Gram-Schmidt, o sea, una descomposición QR.
La **iteración simultánea** (o *iteración ortogonal*) es: empezando
con $`\mathcal{Q}^{(0)} = I`$, para $k = 1, 2, \ldots$

```math
Z^{(k)} = A\,\mathcal{Q}^{(k-1)},
\qquad
Z^{(k)} = \mathcal{Q}^{(k)} R^{(k)} .
```

Compárenlo con la potencia normalizada,
$`\mathbf{z}^{(k)} = A\mathbf{q}^{(k-1)}`$,
$`\mathbf{q}^{(k)} = \mathbf{z}^{(k)}/\|\mathbf{z}^{(k)}\|`$: es lo
mismo con matrices en vez de vectores, y "ortonormalizar" en vez de
"normalizar". (La $\mathcal{Q}$ caligráfica es solo un nombre: la
matriz ortogonal que sale de descomponer $`Z^{(k)}`$. La razón de la
letra rara se ve en la siguiente sección.)

Definamos también el producto (en orden inverso) de todas las $R$:

```math
\mathcal{R}^{(k)} = R^{(k)} R^{(k-1)} \cdots R^{(2)} R^{(1)} ,
\qquad \mathcal{R}^{(0)} = I .
```

El producto de matrices triangulares superiores es triangular superior
(problema 4.2), así que $`\mathcal{R}^{(k)}`$ lo es. Entonces se
cumple algo notable:

```math
A^k = \mathcal{Q}^{(k)}\mathcal{R}^{(k)} .
```

(Cuidado: $`A^k`$ es una potencia; $`\mathcal{Q}^{(k)}`$ es el número
de iteración.) Es decir, la iteración simultánea va construyendo la
**descomposición QR de las potencias** de $A$. Se prueba por inducción.
Para $k = 0$ es $I = I\cdot I$. Si vale para $k - 1$:

```math
A^k = A\,A^{k-1}
= A\,\mathcal{Q}^{(k-1)}\mathcal{R}^{(k-1)}
= Z^{(k)}\mathcal{R}^{(k-1)}
= \mathcal{Q}^{(k)} R^{(k)}\mathcal{R}^{(k-1)}
= \mathcal{Q}^{(k)}\mathcal{R}^{(k)} .
```

(Se usó la hipótesis, luego las dos ecuaciones de la iteración, y al
final la definición de $`\mathcal{R}^{(k)}`$.)

¿Qué tiene que ver con los eigenvectores? Si $A$ es simétrica, sus
eigenvectores son ortonormales y $`\mathcal{Q}^{(k)}`$ tiende a la matriz
de eigenvectores. Si no es simétrica (nuestro caso), los
eigenvectores no son ortogonales, y $`\mathcal{Q}^{(k)}`$ tiende al
"factor ortogonal" $\tilde{Q}$ de la matriz de eigenvectores: la
descomposición QR $V = \tilde{Q}U$, con $U$ triangular superior
(problema 4.32).

**Los eigenvalores.** En la potencia, el eigenvalor salía del cociente
de Rayleigh $`[\mathbf{q}^{(k)}]^T A\mathbf{q}^{(k)}`$. La
generalización a $n$ vectores es

```math
A^{(k)} = [\mathcal{Q}^{(k)}]^T A\,\mathcal{Q}^{(k)} ,
```

una **transformación de similitud ortogonal**: $`A^{(k)}`$ tiene los
mismos eigenvalores que $A$. Y cuando $`\mathcal{Q}^{(k)} \to \tilde{Q}`$:

```math
A^{(k)} \to \tilde{Q}^T A\tilde{Q}
= \tilde{Q}^T A V U^{-1}
= \tilde{Q}^T V \Lambda U^{-1}
= \tilde{Q}^T \tilde{Q} U \Lambda U^{-1}
= U \Lambda U^{-1} .
```

(Se usó $`\tilde{Q} = VU^{-1}`$, luego $AV = V\Lambda$, luego otra vez
$V = \tilde{Q}U$, y por último $`\tilde{Q}^T\tilde{Q} = I`$.) El
resultado es un producto de tres matrices triangulares superiores (la
inversa de una triangular superior también lo es), así que es
**triangular superior**, con $\Lambda$ en la diagonal: los
eigenvalores se leen de la diagonal de $`A^{(k)}`$. Esta es la
descomposición de Schur que buscábamos. Si $A$ es simétrica, $`A^{(k)}`$
tiende a una matriz **diagonal**.

Resumiendo, la iteración simultánea da

```math
A^k = \mathcal{Q}^{(k)}\mathcal{R}^{(k)},
\qquad
A^{(k)} = [\mathcal{Q}^{(k)}]^T A\,\mathcal{Q}^{(k)} .
```

La primera tiene que ver con los eigenvectores y la segunda con los
eigenvalores. Como criterio de paro se puede usar un número fijo de
iteraciones o vigilar cuánto cambia $`A^{(k)}`$.

### El método QR: algoritmo

El método QR es sorprendentemente simple. Empezando con
$`A^{(0)} = A`$, para $k = 1, 2, \ldots$

```math
A^{(k-1)} = Q^{(k)} R^{(k)},
\qquad
A^{(k)} = R^{(k)} Q^{(k)} .
```

Se descompone, y se **multiplica en orden inverso**. Y otra vez, y
otra vez. (Ahora la $Q$ es normal, no caligráfica.)

Dicho así parece magia. Pero es **exactamente la iteración
simultánea**, solo que calcula primero $`A^{(k)}`$ en vez de
$`\mathcal{Q}^{(k)}`$. Para verlo, definimos las mismas
$`\mathcal{R}^{(k)}`$ de antes y

```math
\mathcal{Q}^{(k)} = Q^{(1)} Q^{(2)} \cdots Q^{(k-1)} Q^{(k)} ,
\qquad \mathcal{Q}^{(0)} = I ,
```

(en el orden **contrario** al de las $R$). El producto de matrices
ortogonales es ortogonal (problema 4.2), así que $`\mathcal{Q}^{(k)}`$
lo es. Hay que probar las mismas dos relaciones de antes.

**Segunda relación,** $`A^{(k)} = [\mathcal{Q}^{(k)}]^T A\,\mathcal{Q}^{(k)}`$.
Para $k = 0$ es $`A^{(0)} = A`$. Multiplicando la primera ecuación del
método por $`[Q^{(k)}]^T`$ (que es la inversa de $`Q^{(k)}`$):
$`R^{(k)} = [Q^{(k)}]^T A^{(k-1)}`$. Sustituyendo en la segunda:

```math
A^{(k)} = [Q^{(k)}]^T A^{(k-1)} Q^{(k)} .
```

Cada paso del método QR es una transformación de similitud ortogonal.
Usando la hipótesis de inducción para $`A^{(k-1)}`$:

```math
A^{(k)} = [Q^{(k)}]^T [\mathcal{Q}^{(k-1)}]^T A\,\mathcal{Q}^{(k-1)} Q^{(k)}
= [\mathcal{Q}^{(k)}]^T A\,\mathcal{Q}^{(k)} ,
```

donde se usó que la transpuesta de un producto es el producto de las
transpuestas en orden inverso. Por eso las $Q$ se multiplican en ese
orden.

**Primera relación,** $`A^k = \mathcal{Q}^{(k)}\mathcal{R}^{(k)}`$. Para
$k = 0$ es $I = I\cdot I$. Si vale para $k - 1$:

```math
A^k = A\,A^{k-1}
= A\,\mathcal{Q}^{(k-1)}\mathcal{R}^{(k-1)}
= \mathcal{Q}^{(k-1)} A^{(k-1)}\mathcal{R}^{(k-1)}
= \mathcal{Q}^{(k-1)} Q^{(k)} R^{(k)}\mathcal{R}^{(k-1)}
= \mathcal{Q}^{(k)}\mathcal{R}^{(k)} .
```

En el tercer paso se usó la segunda relación (ya probada) para $k-1$,
multiplicada por $`\mathcal{Q}^{(k-1)}`$ por la izquierda:
$`A\,\mathcal{Q}^{(k-1)} = \mathcal{Q}^{(k-1)}A^{(k-1)}`$. En el cuarto, la
primera ecuación del método QR.

Así que los dos métodos cumplen las mismas relaciones; con
$`\mathcal{Q}^{(0)} = I`$ dan, paso a paso, las mismas matrices. ¿Por qué
preferir el método QR? La iteración simultánea se concentra en los
eigenvectores ($`\mathcal{Q}^{(k)}`$ sale directo); para tener
$`A^{(k)}`$ hay que hacer dos productos de matrices extra. El método
QR calcula directamente $`A^{(k)}`$, de cuya diagonal salen los
eigenvalores; para tener $`\mathcal{Q}^{(k)}`$ habría que ir
acumulando el producto (problema 4.35). Y ya no hay magia: la segunda
ecuación, $`A^{(k)} = R^{(k)}Q^{(k)}`$, **es** una transformación de
similitud ortogonal.

Tal como está, el método QR es poco eficiente. En la práctica se
combina con desplazamientos (como en la potencia inversa) y con otros
trucos.

La parte 3b de [`qr.py`](qr.py) corre los dos métodos lado a lado y
compara las diagonales de $`A^{(k)}`$: coinciden hasta
$`\sim 10^{-13}`$. (Los elementos de arriba de la diagonal difieren
más, $`\sim 10^{-4}`$: la iteración simultánea usa
$`[\mathcal{Q}^{(k)}]^T`$ como si fuera la inversa, y con Gram-Schmidt
clásico $Q$ no es del todo ortogonal.)

### Implementación

Con `descomposicion_qr` y `mat_mat`, el método QR es casi trivial:

```python
def metodo_qr(A, kmax=100):
    A_k = [renglon[:] for renglon in A]
    for k in range(1, kmax):
        Q, R = descomposicion_qr(A_k)
        A_k = mat_mat(R, Q)
    return [A_k[i][i] for i in range(len(A))]
```

Como en el libro, se hace un número fijo de iteraciones, y por eso no
hace falta guardar los superíndices: `A_k`, `Q` y `R` se van
reemplazando. Se copia `A` al principio para no modificar la matriz
de quien llama a la función.

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
eigenvectores podríamos acumular $`\mathcal{Q}^{(k)}`$, pero para una
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
  `eigenvalores[i]`), que con listas de listas es lo más cómodo. El
  libro, con NumPy, los regresa como columnas de una matriz.
- `eigenpares` es tan corta porque usa todo lo anterior:
  `metodo_qr` usa `descomposicion_qr`, y `potencia_inversa` usa
  `descomposicion_lu` y las sustituciones.

La parte 5 de [`qr.py`](qr.py) la prueba. Los residuos
$`\|A\mathbf{v} - \lambda\mathbf{v}\|/\|\mathbf{v}\|`$ son todos
$`\sim 10^{-15}`$, la suma de los eigenvalores da la traza y su producto
el determinante, y el primer eigenvector es el de la potencia. Igual
que en el libro, la potencia inversa cambia un poco el eigenvalor más
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
- Modifiquen `metodo_qr` para que también acumule
  $`\mathcal{Q}^{(k)}`$ (problema 4.35). Con la matriz simétrica $S$,
  comparen sus columnas con los eigenvectores de `eigenpares`.
- Agreguen a `metodo_qr` un criterio de paro: por ejemplo, que la
  parte de abajo de la diagonal sea menor que una tolerancia.
- Prueben `eigenpares` con `eps = 0`.
## Contenido
- [`fiscomp/algebra_lineal.py`](../../fiscomp/algebra_lineal.py): los
  métodos de la potencia (`producto_punto`, `normalizar`,
  `cociente_rayleigh`, `potencia`, `potencia_inversa`) y el método QR
  (`transpuesta`, `mat_mat`, `descomposicion_qr`, `metodo_qr`,
  `eigenpares`). Es el mismo módulo de la unidad 09: la potencia
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

Se corren desde la raíz del repositorio (con el `.venv` activado), por
ejemplo `python3 unidades/10_valores_y_vectores_propios/qr.py`.
