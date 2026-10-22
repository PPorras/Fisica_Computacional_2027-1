# Unidad 09 — Matrices

## Temas
- **Física computacional:** por qué el álgebra lineal aparece en todos
  lados en física; los dos grandes problemas del tema, sistemas de
  ecuaciones lineales $A\mathbf{x} = \mathbf{b}$ y el problema de
  eigenvalores $A\mathbf{v} = \lambda\mathbf{v}$; análisis de error
  *a priori*: normas de vectores y de matrices, número de condición
  para sistemas lineales, número de condición para eigenvalores
  simples y sensibilidad de los eigenvectores; solución de sistemas
  triangulares (sustitución hacia adelante y hacia atrás) y conteo de
  operaciones; eliminación gaussiana y descomposición LU (Doolittle),
  su costo, y cómo reutilizar LU para muchos lados derechos (la
  inversa, el número de condición y el determinante); inestabilidad
  sin mal condicionamiento y pivoteo parcial; el método iterativo de
  Jacobi; el problema de eigenvalores numéricamente (diagonalización,
  por qué no usar el polinomio característico), el método de la
  potencia, el cociente de Rayleigh y el método de la potencia
  inversa con desplazamiento; la descomposición QR (Gram-Schmidt),
  las transformaciones de similitud, la iteración simultánea y el
  método QR para todos los eigenvalores; todos los eigenpares
  combinando QR con la potencia inversa.
- Corresponde a las secciones 4.1 a 4.4.4 del libro de Gezerlis.

## Notación

Los índices van de $0$ a $n-1$, como en el libro y en Python. Las
matrices se escriben con mayúsculas ($A$, $L$, $U$) y los vectores en
negritas minúsculas ($\mathbf{x}$, $\mathbf{b}$). Cuando hace falta,
los dos índices de un elemento se separan con una coma: $A_{1,n-1}$.

En el código, una matriz es una **lista de listas** (un renglón por
sublista, `A[i][j]`) y un vector es una lista simple (`x[i]`). Es el
mismo formato que guarda la clase `Matrix` de la práctica 2 en su
atributo `.data`. En este tema **no usamos NumPy**: todos los métodos
se programan desde cero.

## Motivación: ejemplos de física

El álgebra lineal aparece en casi todas las áreas de la física. Tres
ejemplos de licenciatura, que no requieren cómputo pesado pero sí
usan los mismos conceptos que vamos a necesitar:

**1. Rotaciones en dos dimensiones.** Rotar el punto
$\mathbf{r} = (x, y)^T$ un ángulo $\theta$ en sentido antihorario
alrededor del origen da el punto $\mathbf{r}' = (x', y')^T$, con

```math
\begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix}
\begin{pmatrix} x \\ y \end{pmatrix}
= \begin{pmatrix} x' \\ y' \end{pmatrix}.
```

Si conocemos $\mathbf{r}'$ y queremos $\mathbf{r}$, hay que resolver
un sistema de dos ecuaciones lineales. Nótese que la matriz de
rotación **no es simétrica** (sus dos elementos fuera de la diagonal
no son iguales). Las matrices simétricas aparecen por todos lados en
física, pero no siempre; por eso estudiaremos matrices generales.

**2. Potenciales electrostáticos.** Tenemos $n$ cargas $q_j$
*desconocidas* en posiciones $\mathbf{R}_j$ conocidas, y medimos el
potencial $\phi(\mathbf{r}_i)$ en $n$ puntos $\mathbf{r}_i$ conocidos.
Por el principio de superposición,

```math
\phi(\mathbf{r}_i) = \sum_{j=0}^{n-1} \frac{k}{|\mathbf{r}_i - \mathbf{R}_j|}\, q_j,
\qquad i = 0, 1, \ldots, n-1.
```

Es el problema **inverso** del proyecto
[`simulacion_potencial`](../../simulaciones/simulacion_potencial/): ahí
conocíamos las cargas y calculábamos el potencial; aquí conocemos el
potencial y queremos las cargas. Para cuatro cargas es el sistema
$4\times 4$

```math
\begin{pmatrix}
\frac{k}{|\mathbf{r}_0-\mathbf{R}_0|} & \frac{k}{|\mathbf{r}_0-\mathbf{R}_1|} & \frac{k}{|\mathbf{r}_0-\mathbf{R}_2|} & \frac{k}{|\mathbf{r}_0-\mathbf{R}_3|} \\
\frac{k}{|\mathbf{r}_1-\mathbf{R}_0|} & \frac{k}{|\mathbf{r}_1-\mathbf{R}_1|} & \frac{k}{|\mathbf{r}_1-\mathbf{R}_2|} & \frac{k}{|\mathbf{r}_1-\mathbf{R}_3|} \\
\frac{k}{|\mathbf{r}_2-\mathbf{R}_0|} & \frac{k}{|\mathbf{r}_2-\mathbf{R}_1|} & \frac{k}{|\mathbf{r}_2-\mathbf{R}_2|} & \frac{k}{|\mathbf{r}_2-\mathbf{R}_3|} \\
\frac{k}{|\mathbf{r}_3-\mathbf{R}_0|} & \frac{k}{|\mathbf{r}_3-\mathbf{R}_1|} & \frac{k}{|\mathbf{r}_3-\mathbf{R}_2|} & \frac{k}{|\mathbf{r}_3-\mathbf{R}_3|}
\end{pmatrix}
\begin{pmatrix} q_0 \\ q_1 \\ q_2 \\ q_3 \end{pmatrix}
= \begin{pmatrix} \phi(\mathbf{r}_0) \\ \phi(\mathbf{r}_1) \\ \phi(\mathbf{r}_2) \\ \phi(\mathbf{r}_3) \end{pmatrix},
```

cuya matriz tampoco es simétrica en general.

**3. Momentos principales de inercia.** Para un cuerpo rígido que rota
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

## Los dos problemas a resolver

### Sistemas de ecuaciones lineales

Tenemos $n$ incógnitas $x_i$, $n\times n$ coeficientes $A_{ij}$ y $n$
constantes $b_i$:

```math
\begin{pmatrix}
A_{00} & A_{01} & \cdots & A_{0,n-1} \\
A_{10} & A_{11} & \cdots & A_{1,n-1} \\
\vdots & \vdots & \ddots & \vdots \\
A_{n-1,0} & A_{n-1,1} & \cdots & A_{n-1,n-1}
\end{pmatrix}
\begin{pmatrix} x_0 \\ x_1 \\ \vdots \\ x_{n-1} \end{pmatrix}
= \begin{pmatrix} b_0 \\ b_1 \\ \vdots \\ b_{n-1} \end{pmatrix},
\qquad\text{o sea}\qquad
A\mathbf{x} = \mathbf{b}.
```

$A$ es la **matriz de coeficientes**. Suponemos que $|A| \neq 0$ (la
matriz es no singular, sus columnas son linealmente independientes),
de modo que la solución existe y es única. Parece un problema muy
sencillo, pero le vamos a dedicar bastante tiempo.

Casi siempre trabajaremos con la **matriz aumentada** $(A|\mathbf{b})$,
que pone juntos los elementos de $A$ y de $\mathbf{b}$ (sin escribir
explícitamente $\mathbf{x}$). Sobre ella se permiten las siguientes
**operaciones elementales de renglón**, que cambian $A$ y
$\mathbf{b}$ pero **no cambian la solución** $\mathbf{x}$:

- **Escalamiento:** multiplicar un renglón (una ecuación) por una
  constante. Multiplica $|A|$ por esa constante.
- **Pivoteo:** intercambiar dos renglones. Cambia el signo de $|A|$.
- **Eliminación:** sustituir un renglón por la suma de ese renglón con
  un múltiplo de otro. No cambia $|A|$.

Como se aplican sobre la matriz aumentada, al intercambiar dos
renglones de $A$ hay que intercambiar también los dos elementos
correspondientes de $\mathbf{b}$.

### El problema de eigenvalores

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

Como en la unidad 06, buscamos cotas de error *pesimistas* (el peor
caso). No vamos a analizar un método específico, sino **qué tan
sensible es el problema mismo** a cambios pequeños en los datos de
entrada: si el problema está **bien condicionado** o **mal
condicionado**. Esto basta, porque un resultado clásico (Wilkinson,
años sesenta, para eliminación gaussiana) muestra que los errores de
redondeo de un buen método equivalen a perturbar un poco los datos de
entrada (análisis de error *hacia atrás*, *backward error analysis*).

### De *a posteriori* a *a priori*: el ejemplo de Kahan

Hay dos maneras de analizar errores:

- ***A posteriori***: ya tenemos una solución y queremos saber qué tan
  buena es.
- ***A priori***: antes de resolver, queremos saber qué tan difícil es
  el problema.

Consideremos el sistema (debido a W. Kahan)

```math
(A|\mathbf{b}) = \left(\begin{array}{cc|c}
1.2969 & 0.8648 & 0.8642 \\
0.2161 & 0.1441 & 0.1440
\end{array}\right),
```

y supongamos que alguien nos da la solución aproximada
$\tilde{\mathbf{x}}^T = (0.9911, -0.4870)$. (Escribimos la transpuesta
para ahorrar espacio.)

**A posteriori.** Una forma de medir qué tan buena es: el **vector
residuo**

```math
\mathbf{r} = \mathbf{b} - A\tilde{\mathbf{x}}.
```

Si $\mathbf{x}$ es la solución exacta, $A\mathbf{x} = \mathbf{b}$ y el
residuo es cero; una buena aproximación debería dar un residuo chico.
Aquí $\mathbf{r}^T = (10^{-8}, -10^{-8})$: diminuto. Uno concluiría
que $\tilde{\mathbf{x}}$ es excelente. Pero la solución exacta es

```math
\mathbf{x}^T = (2, -2),
```

así que $\tilde{\mathbf{x}}$ **no tiene ni una sola cifra
significativa correcta**. Un residuo chico no garantiza una solución
buena.

**A priori.** ¿Podíamos haber sabido de antemano que este problema es
patológico? Una forma: perturbar un poco los datos de entrada, como si
no conociéramos los coeficientes de $A$ con toda precisión (en física,
$A$ y $\mathbf{b}$ muchas veces vienen de cálculos previos o de
mediciones, con su propio error). Si cambiamos **un solo** elemento,
$A_{00} = 1.2969 \to 1.2970$ (menos de 0.01%), la solución pasa de
$(2, -2)$ a aproximadamente $(0.0014, 0.9972)$. Un cambio minúsculo en
los datos produce un cambio enorme en la solución: el problema está
**mal condicionado**.

¿Qué hace que este problema se porte tan mal? Una idea que aparece a
veces en los libros: como una matriz singular tiene determinante cero,
quizá una matriz con determinante "cercano a cero" está "cerca" de ser
singular. Aquí, efectivamente, $\det(A) = 10^{-8}$. Pero, ¿qué quiere
decir "chico"? Intuitivamente debería compararse contra el tamaño de
los elementos de la matriz. Para eso necesitamos una forma de medir el
tamaño de una matriz: las normas.

### Normas de matrices y de vectores

Una **norma de matriz** mide la magnitud de $A$ con un solo número. Se
escribe con doble barra, $`\|A\|`$ (la barra sencilla $|A|$ es el
determinante, o el valor absoluto de un número). Usaremos dos:

**Norma de Frobenius:**

```math
\|A\|_F = \sqrt{\sum_{i=0}^{n-1}\sum_{j=0}^{n-1} |A_{ij}|^2}.
```

**Norma infinito** (máxima suma de renglón):

```math
\|A\|_\infty = \max_{0\le i\le n-1} \sum_{j=0}^{n-1} |A_{ij}|.
```

Cualquier norma de matrices cuadradas cumple:

```math
\begin{aligned}
&\|A\| \ge 0, \\
&\|A\| = 0 \text{ si y solo si todos los } A_{ij} = 0, \\
&\|kA\| = |k|\,\|A\|, \\
&\|A + B\| \le \|A\| + \|B\| \quad\text{(desigualdad del triángulo)}, \\
&\|AB\| \le \|A\|\,\|B\|.
\end{aligned}
```

Una norma es un **número**, no una matriz (igual que el determinante).

Para vectores usaremos:

- **Norma euclídea:** $`\|\mathbf{x}\|_E = \sqrt{\sum_{i=0}^{n-1} |x_i|^2}`$.
- **Norma infinito** (de máxima magnitud): $`\|\mathbf{x}\|_\infty = \max_{0\le i\le n-1} |x_i|`$.

Todas están en [`fiscomp/algebra_lineal.py`](../../fiscomp/algebra_lineal.py)
(`norma_frobenius`, `norma_infinito`, `norma_euclidea`,
`norma_infinito_vec`).

### ¿Determinante chico?

Con normas, el criterio "determinante chico" se podría escribir
$`|\det(A)| \ll \|A\|`$. En el ejemplo de Kahan,
$`|\det(A)| = 10^{-8} \ll \|A\|_\infty \approx 2.16`$, y el problema sí
está mal condicionado. **Pero el criterio es falso.** Tomemos
$`D = 0.1\, I`$ de $20\times 20$: $`\det(D) = 0.1^{20} = 10^{-20}`$ y
$`\|D\|_\infty = 0.1`$, así que $`|\det(D)| \ll \|D\|`$... y sin embargo
$D\mathbf{x} = \mathbf{b}$ se resuelve trivialmente ($`\mathbf{x} = 10\,\mathbf{b}`$)
y cambiar un poco $D$ o $\mathbf{b}$ cambia un poco
$\mathbf{x}$. El determinante es sensible a la escala y a la dimensión
de la matriz de una forma que no tiene nada que ver con el
condicionamiento. Necesitamos otra medida.

### Número de condición para sistemas lineales

Partimos del problema sin perturbar,

```math
A\mathbf{x} = \mathbf{b},
```

y cambiamos un poco $A$ (dejando $\mathbf{b}$ fija). La solución
cambia también:

```math
(A + \Delta A)(\mathbf{x} + \Delta\mathbf{x}) = \mathbf{b}.
```

(Como en la unidad 06, el error absoluto es "aproximado menos exacto":
esto es $\tilde{A}\tilde{\mathbf{x}} = \mathbf{b}$.) Desarrollando y
usando $A\mathbf{x} = \mathbf{b}$:

```math
A\,\Delta\mathbf{x} = -\Delta A\,(\mathbf{x} + \Delta\mathbf{x})
\quad\Longrightarrow\quad
\Delta\mathbf{x} = -A^{-1}\Delta A\,(\mathbf{x} + \Delta\mathbf{x}).
```

Tomando normas y usando $`\|kA\| = |k|\|A\|`$ (con $k=-1$) y
$`\|AB\| \le \|A\|\|B\|`$ (dos veces):

```math
\|\Delta\mathbf{x}\| \le \|A^{-1}\|\,\|\Delta A\|\,\|\mathbf{x} + \Delta\mathbf{x}\|
\le \|A^{-1}\|\,\|\Delta A\|\,\left(\|\mathbf{x}\| + \|\Delta\mathbf{x}\|\right)
\approx \|A^{-1}\|\,\|\Delta A\|\,\|\mathbf{x}\|,
```

donde en el segundo paso usamos la desigualdad del triángulo y en el
último despreciamos el término de segundo orden
$`\|\Delta A\|\,\|\Delta\mathbf{x}\|`$ (producto de dos cantidades
chicas). Dividiendo entre $`\|\mathbf{x}\|`$, y multiplicando y
dividiendo por $`\|A\|`$:

```math
\frac{\|\Delta\mathbf{x}\|}{\|\mathbf{x}\|} \le \|A\|\,\|A^{-1}\|\,\frac{\|\Delta A\|}{\|A\|}.
```

Es decir: una cota al cambio *relativo* en $A$ se traduce en una cota
al cambio *relativo* en $\mathbf{x}$, y el factor que los conecta es el
**número de condición**

```math
\kappa(A) = \|A\|\,\|A^{-1}\|.
```

- Si $\kappa(A)$ es de orden 1, una perturbación chica no se
  amplifica: el problema está **bien condicionado**. (Siempre
  $\kappa(A) \ge 1$.)
- Si $\kappa(A) \gg 1$, una perturbación chica se puede amplificar
  mucho: el problema está **mal condicionado**.

Para el ejemplo de Kahan, $\kappa_\infty(A) \approx 3.3\times 10^{8}$:
un error relativo de $10^{-4}$ en $A$ puede volverse un error relativo
de hasta $10^{4}$ en $\mathbf{x}$. Y el propio redondeo al guardar los
datos ($\sim 10^{-16}$) se puede amplificar hasta $\sim 10^{-8}$. Para
$`D = 0.1\,I`$, en cambio, $\kappa(D) = 0.1 \times 10 = 1$. **El número
de condición no involucra al determinante.**

Algunas observaciones:

- $\kappa(A)$ mide a la vez qué tan bien condicionado está resolver
  $A\mathbf{x} = \mathbf{b}$ y qué tan bien condicionado está invertir
  $A$ (conceptualmente, aunque no en la práctica, resolver el sistema
  es invertir la matriz).
- Su valor preciso depende de la norma que se use, pero su orden de
  magnitud no.
- Tiene algo de trampa: para calcular $\kappa(A)$ hace falta $A^{-1}$,
  que se calcula con los mismos métodos cuya confiabilidad queremos
  evaluar, y que cuesta más ($`O(n^3)`$) que resolver el sistema. En la
  práctica se usan métodos que *estiman* $\kappa(A)$ (dentro de un
  factor de 10 o así) con solo $O(n^2)$ operaciones.
- Se puede hacer lo mismo perturbando $\mathbf{b}$ en vez de $A$ (o
  ambas): el resultado vuelve a involucrar a $\kappa(A)$.

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

## Matrices triangulares

Con el análisis de error listo, empezamos a *resolver*
$A\mathbf{x} = \mathbf{b}$. Nos concentramos en **métodos directos**,
que transforman el problema original en otro más fácil de resolver.
(Los **métodos iterativos** parten de una propuesta de solución y la
refinan hasta converger; son útiles sobre todo para matrices
*ralas*, *sparse*, con muchos ceros.)

El caso más sencillo es el de las **matrices triangulares**, con ceros
arriba o abajo de la diagonal. No es un problema de juguete: los
métodos generales (eliminación gaussiana y descomposición LU, en las
secciones siguientes) transforman el problema hasta dejarlo en términos de
una o dos matrices triangulares. Aquí construimos esa base.

(En las bibliotecas profesionales no se guardan los ceros de una
matriz triangular, que no aportan información; incluso es común
guardar una triangular inferior y una superior juntas en una sola
matriz. Aquí, por claridad, guardamos la matriz completa, con todo y
ceros.)

### Sustitución hacia adelante

Sea $L$ triangular inferior (*lower*). Queremos resolver
$L\mathbf{x} = \mathbf{b}$. (Aquí $\mathbf{b}$ es simplemente "el
vector del lado derecho" y $\mathbf{x}$ "la incógnita"; más adelante
les daremos otros nombres.) Para $3\times 3$:

```math
\begin{pmatrix}
L_{00} & 0 & 0 \\
L_{10} & L_{11} & 0 \\
L_{20} & L_{21} & L_{22}
\end{pmatrix}
\begin{pmatrix} x_0 \\ x_1 \\ x_2 \end{pmatrix}
= \begin{pmatrix} b_0 \\ b_1 \\ b_2 \end{pmatrix},
\qquad\text{o sea}\qquad
\begin{aligned}
L_{00}x_0 &= b_0 \\
L_{10}x_0 + L_{11}x_1 &= b_1 \\
L_{20}x_0 + L_{21}x_1 + L_{22}x_2 &= b_2
\end{aligned}
```

La primera ecuación da $x_0$; con él, la segunda da $x_1$; con ambos,
la tercera da $x_2$:

```math
x_0 = \frac{b_0}{L_{00}}, \qquad
x_1 = \frac{b_1 - L_{10}x_0}{L_{11}}, \qquad
x_2 = \frac{b_2 - L_{20}x_0 - L_{21}x_1}{L_{22}}.
```

Se llama **sustitución hacia adelante** (*forward substitution*)
porque se empieza por la primera ecuación y se avanza. Para
$n\times n$:

```math
x_i = \frac{1}{L_{ii}}\left(b_i - \sum_{j=0}^{i-1} L_{ij}\,x_j\right),
\qquad i = 0, 1, \ldots, n-1,
```

entendiendo que la suma no tiene términos si $i = 0$, tiene uno si
$i = 1$, etc.

### Sustitución hacia atrás

Si en cambio $U$ es triangular superior (*upper*), $U\mathbf{x} = \mathbf{b}$
es, para $3\times 3$,

```math
\begin{aligned}
U_{00}x_0 + U_{01}x_1 + U_{02}x_2 &= b_0 \\
U_{11}x_1 + U_{12}x_2 &= b_1 \\
U_{22}x_2 &= b_2
\end{aligned}
```

Ahora se empieza por la **última** ecuación y se avanza hacia atrás:

```math
x_2 = \frac{b_2}{U_{22}}, \qquad
x_1 = \frac{b_1 - U_{12}x_2}{U_{11}}, \qquad
x_0 = \frac{b_0 - U_{01}x_1 - U_{02}x_2}{U_{00}}.
```

Es la **sustitución hacia atrás** (*back substitution*). Para
$n\times n$:

```math
x_i = \frac{1}{U_{ii}}\left(b_i - \sum_{j=i+1}^{n-1} U_{ij}\,x_j\right),
\qquad i = n-1, n-2, \ldots, 1, 0,
```

con la suma vacía si $i = n-1$, de un término si $i = n-2$, etc.

Nótese que ambas fórmulas dividen entre los elementos de la diagonal:
solo funcionan si ninguno es cero. (Para una matriz triangular, el
determinante es el producto de la diagonal, así que eso es lo mismo
que pedir que la matriz sea no singular.)

### Implementación

Las dos fórmulas se traducen casi literalmente a Python, en
[`fiscomp/algebra_lineal.py`](../../fiscomp/algebra_lineal.py):

```python
def sustitucion_adelante(L, b):
    n = len(b)
    x = [0.0] * n
    for i in range(n):
        suma = 0.0
        for j in range(i):  # para i = 0 este ciclo no hace nada
            suma += L[i][j] * x[j]
        x[i] = (b[i] - suma) / L[i][i]
    return x


def sustitucion_atras(U, b):
    n = len(b)
    x = [0.0] * n
    for i in reversed(range(n)):
        suma = 0.0
        for j in range(i + 1, n):  # para i = n-1 este ciclo no hace nada
            suma += U[i][j] * x[j]
        x[i] = (b[i] - suma) / U[i][i]
    return x
```

Algunos detalles:

- El tamaño del problema se toma del vector `b`.
- La solución `x` se crea primero llena de ceros y luego se llena
  elemento por elemento (en vez de empezar con una lista vacía y usar
  `append`). Así, al calcular `x[i]`, los `x[j]` que ya se conocen
  están en su lugar.
- `range(i)` no produce nada cuando `i = 0`, así que el caso de la
  primera ecuación (sin suma) sale solo; lo mismo con
  `range(i + 1, n)` cuando `i = n - 1`.
- `reversed(range(n))` recorre `n-1, n-2, ..., 0`.
- Solo se leen los elementos de la diagonal y de un lado de ella, así
  que las funciones dan el mismo resultado si se les pasa la matriz
  completa en vez de su parte triangular.

En el libro (código 4.1) estas funciones se escriben con NumPy, y la
suma interna se hace con un producto punto, `L[i,:i] @ xs[:i]`, sin el
ciclo sobre `j`. Nosotros escribimos ese ciclo explícitamente.

Para probarlas, `triangulares.py` arma la matriz de prueba del libro,
$A_{ij} = \sqrt{21 + 4i + j}$ ($4\times 4$, no simétrica), toma su
parte triangular inferior y superior (`triangular_inferior`,
`triangular_superior`, el equivalente de `numpy.tril` y `numpy.triu`)
y comprueba la solución de dos formas, sin cajas negras: calculando el
residuo, y resolviendo un sistema cuya solución exacta conocemos
(elegimos $\mathbf{x}$, calculamos $\mathbf{b} = L\mathbf{x}$ y vemos
si recuperamos $\mathbf{x}$). En ambos casos el error es del orden de
$\epsilon_{\text{mach}}$: con matrices tan chicas y bien portadas no
esperábamos otra cosa.

### Conteo de operaciones

En álgebra lineal es muy común **contar cuántas operaciones de punto
flotante** (*flops*: sumas, restas, multiplicaciones, divisiones)
requiere un cálculo, en función del tamaño $n$ del problema. Si el
costo crece muy rápido con $n$ (por ejemplo, exponencialmente), no
podremos resolver problemas mucho más grandes que los actuales; si
crece como una potencia baja, sí.

Se usa la notación $O$ (unidad 08): un método $`O(n^3)`$ es mejor que
uno $O(n^4)$ para $n$ grande, sin importar los prefactores. Cuando se
cuenta con más cuidado interesa también el prefactor ($2n^3$ es mejor
que $4n^3$), pero los términos de grado menor se suelen tirar:
$2n^3 - 7n^2 + 5n \sim 2n^3$.

**Producto matriz-vector.** $y_i = \sum_{j=0}^{n-1} A_{ij}x_j$. Cada
$y_i$ requiere $n$ multiplicaciones y $n-1$ sumas; hay $n$ de ellos,
así que en total son $n^2$ multiplicaciones y $n(n-1)$ sumas:

```math
2n^2 - n \sim 2n^2 = O(n^2)\ \text{operaciones}.
```

**Sustitución hacia adelante.** Cada $x_i$ requiere una división ($n$
en total). Además, para $x_i$ hay $i$ multiplicaciones
($L_{ij}x_j$ para $j = 0, \ldots, i-1$) e $i$ sumas/restas ($i-1$
sumas para juntar los $i$ productos, más la resta de $b_i$). Sumando
sobre $i$, con $\sum_{i=0}^{n-1} i = \frac{(n-1)n}{2}$:

```math
\text{sumas/restas: } \frac{n^2 - n}{2},
\qquad
\text{multiplicaciones/divisiones: } n + \frac{n^2 - n}{2} = \frac{n^2 + n}{2}.
```

En total, **exactamente $n^2$ operaciones**. Se podría escribir
$O(n^2)$, pero eso es menos informativo: el conteo explícito nos dice
que el prefactor es exactamente 1. (La sustitución hacia atrás cuesta
lo mismo.)

`triangulares.py` comprueba ambos conteos agregando un contador a cada
operación, y mide el tiempo de la sustitución hacia adelante al
duplicar $n$: como el costo es $\propto n^2$, el tiempo se cuadruplica.

Para comparar: el producto de dos vectores es $O(n)$ y el de dos
matrices $`O(n^3)`$ (calculen los prefactores exactos, junto con los
términos de grado menor). Y, como se ve en la sección "Eliminación
gaussiana", resolver un sistema general con ella es $`O(n^3)`$: mucho más caro que resolver uno
triangular.

## Eliminación gaussiana

Pasamos ahora al caso general: resolver $A\mathbf{x} = \mathbf{b}$
cuando $A$ **no** es triangular. El primer método es la **eliminación
gaussiana** (aunque se usaba en China dos mil años antes, y Newton la
conocía más de un siglo antes que Gauss). Usa solo la tercera
operación elemental de renglón de la sección "Sistemas de ecuaciones
lineales": sustituir un renglón por ese renglón más un múltiplo de
otro, que no cambia la solución.

El método tiene dos fases:

1. **Eliminación:** aplicar esa operación muchas veces, hasta que la
   matriz de coeficientes quede triangular superior.
2. **Sustitución hacia atrás:** resolver el sistema triangular que
   quedó, con `sustitucion_atras`, que ya tenemos.

### Ejemplo $3\times 3$

Antes del caso general, resolvamos un ejemplo a mano:

```math
\begin{aligned}
2x_0 + x_1 + x_2 &= 8 \\
x_0 + x_1 - 2x_2 &= -2 \\
5x_0 + 10x_1 + 5x_2 &= 10
\end{aligned}
\qquad\text{o sea}\qquad
\begin{pmatrix} 2 & 1 & 1 \\ 1 & 1 & -2 \\ 5 & 10 & 5 \end{pmatrix}
\begin{pmatrix} x_0 \\ x_1 \\ x_2 \end{pmatrix}
= \begin{pmatrix} 8 \\ -2 \\ 10 \end{pmatrix}.
```

Trabajamos con la matriz aumentada, que contiene todo lo que importa
(la $\mathbf{x}$ queda implícita):

```math
(A|\mathbf{b}) = \left(\begin{array}{ccc|c}
2 & 1 & 1 & 8 \\
1 & 1 & -2 & -2 \\
5 & 10 & 5 & 10
\end{array}\right).
```

Nótese que $A$ no es simétrica. La operación que vamos a repetir es

```math
\text{nuevo renglón } i = \text{renglón } i - \text{coeficiente} \times \text{renglón } j.
```

El renglón $j$ se llama **renglón pivote**; el renglón $i$ es el que
estamos transformando. El coeficiente se escoge para que el primer
elemento distinto de cero del renglón $i$ se vuelva $0$.

**Pivote $j = 0$.** Para $i = 1$ el coeficiente es $0.5$, porque
$1 - 0.5\times 2 = 0$. La operación se hace con **todo** el renglón,
incluyendo el elemento de $\mathbf{b}$:

```math
\left(\begin{array}{ccc|c}
2 & 1 & 1 & 8 \\
0 & 0.5 & -2.5 & -6 \\
5 & 10 & 5 & 10
\end{array}\right).
```

Para $i = 2$ el coeficiente es $2.5$, porque $5 - 2.5\times 2 = 0$:

```math
\left(\begin{array}{ccc|c}
2 & 1 & 1 & 8 \\
0 & 0.5 & -2.5 & -6 \\
0 & 7.5 & 2.5 & -10
\end{array}\right).
```

Ya terminamos con el pivote $j = 0$: la columna $0$ tiene ceros abajo
de la diagonal.

**Pivote $j = 1$.** Siempre se usa la versión **más reciente** de la
matriz, así que el renglón pivote es $`(0,\ 0.5,\ -2.5 \,|\, -6)`$. Los
renglones que se transforman están siempre abajo del pivote; aquí
solo queda $i = 2$. El coeficiente es $15$, porque
$7.5 - 15\times 0.5 = 0$:

```math
\left(\begin{array}{ccc|c}
2 & 1 & 1 & 8 \\
0 & 0.5 & -2.5 & -6 \\
0 & 0 & 40 & 80
\end{array}\right).
```

La matriz de coeficientes ya es triangular superior: terminó la
eliminación. Con la sustitución hacia atrás:

```math
x_2 = \frac{80}{40} = 2, \qquad
x_1 = \frac{-6 - (-2.5)\times 2}{0.5} = -2, \qquad
x_0 = \frac{8 - 1\times(-2) - 1\times 2}{2} = 4.
```

`eliminacion_gaussiana_lu.py` repite este ejemplo e imprime la matriz
aumentada después de cada paso.

### Caso general

Para una matriz $n\times n$, la eliminación modifica $A$ y
$\mathbf{b}$ hasta que $A$ queda triangular. A media eliminación, justo
cuando el renglón $j$ se convierte en pivote, la matriz aumentada se ve
así, con $n = 5$ y $j = 2$, donde $\ast$ es un número cualquiera:

```math
\left(\begin{array}{ccccc|c}
\ast & \ast & \ast & \ast & \ast & \ast \\
0 & \ast & \ast & \ast & \ast & \ast \\
0 & 0 & A_{22} & \ast & \ast & \ast \\
0 & 0 & \ast & \ast & \ast & \ast \\
0 & 0 & \ast & \ast & \ast & \ast
\end{array}\right).
```

Los renglones de arriba del pivote ya están listos y los de abajo
todavía tienen que transformarse. (Los valores son los actuales, ya
modificados por los pasos anteriores; el primer renglón nunca cambia.)

- El renglón pivote recorre $j = 0, 1, \ldots, n-2$: el último renglón
  que se transforma es el último, así que el último pivote es el
  penúltimo renglón.
- Para cada pivote, los renglones que se transforman son los de abajo:
  $i = j+1, j+2, \ldots, n-1$.
- El primer elemento distinto de cero del renglón $i$ es $A_{ij}$, y
  el del renglón $j$ es $A_{jj}$. Como
  $`A_{ij} - (A_{ij}/A_{jj})\,A_{jj} = 0`$, el coeficiente es

```math
\text{coeficiente} = \frac{A_{ij}}{A_{jj}}.
```

$A_{jj}$ se llama **elemento pivote**: es el que se divide para
eliminar los primeros elementos de los renglones de abajo. Con ese
coeficiente, el renglón $i$ se actualiza elemento por elemento:

```math
\begin{aligned}
A_{ik} &\leftarrow A_{ik} - \text{coeficiente}\times A_{jk}, \qquad k = j, j+1, \ldots, n-1, \\
b_i &\leftarrow b_i - \text{coeficiente}\times b_j.
\end{aligned}
```

(Las columnas $k < j$ ya son cero en ambos renglones, así que no hace
falta tocarlas.) Al final, $A$ es triangular superior y se aplica la
sustitución hacia atrás.

### Implementación

En [`fiscomp/algebra_lineal.py`](../../fiscomp/algebra_lineal.py):

```python
def eliminacion_gaussiana(A, b):
    A = [renglon[:] for renglon in A]
    b = b[:]
    n = len(b)
    for j in range(n - 1):
        for i in range(j + 1, n):
            coeficiente = A[i][j] / A[j][j]
            for k in range(j, n):
                A[i][k] -= coeficiente * A[j][k]
            b[i] -= coeficiente * b[j]
    return sustitucion_atras(A, b)
```

Algunos detalles:

- Las dos primeras líneas hacen **copias** de `A` y de `b`. Sin ellas,
  la función modificaría las listas de quien la llama (las listas se
  pasan por referencia), y después de resolver el sistema ya no
  tendríamos la matriz original. Es ineficiente, pero más seguro.
- Los dos ciclos externos son exactamente $j = 0, \ldots, n-2$ e
  $i = j+1, \ldots, n-1$. Aquí empieza a valer la pena haber numerado
  todo desde $0$.
- En el libro, con NumPy, el ciclo sobre `k` se escribe en una línea,
  `A[i,j:] -= coeff*A[j,j:]`: se actualiza todo un pedazo del renglón
  a la vez. Nosotros escribimos ese tercer ciclo explícitamente.
- `k` empieza en `j`, no en `j + 1`: calculamos $A_{ij}$ aunque ya
  sabemos que va a dar $0$. Es una operación de más por renglón, a
  cambio de un código más parecido a la fórmula.
- No revisamos si `A[i][j]` ya es cero para saltarnos el renglón:
  estaríamos haciendo esa comparación todo el tiempo para un caso que
  casi nunca ocurre (en matrices *ralas*, con muchos ceros, se usan
  otros métodos).
- Al final reutilizamos `sustitucion_atras`.

Con la matriz de prueba del libro, $A_{ij} = \sqrt{21 + 4i + j}$ de
$4\times 4$, el residuo es diminuto, del orden de $10^{-11}$. Pero si
resolvemos un sistema con solución conocida, eligiendo
$\mathbf{x} = (1, 2, 3, 4)$ y calculando $\mathbf{b} = A\mathbf{x}$,
recuperamos $\mathbf{x}$ solo con unos 8 dígitos correctos, no 16. El
libro lo deja como pregunta abierta; nosotros ya tenemos con qué
contestarla: es el número de condición. En la sección "Descomposición
LU" calculamos $A^{-1}$ y resulta $\kappa(A) \approx 4\times 10^{8}$.
Como vimos en "Análisis de error", el error relativo en $\mathbf{x}$
puede ser hasta $\kappa(A)$ veces el error relativo de los datos, que
es del orden de $\epsilon_{\text{mach}} \approx 10^{-16}$:
$`4\times 10^{8} \times 10^{-16} \approx 4\times 10^{-8}`$, justo lo que
observamos. No es culpa del método, sino de la matriz: sus renglones
son casi iguales.

### Conteo de operaciones

Separamos el costo en dos partes: (a) la eliminación, que convierte
$A$ en triangular y modifica $\mathbf{b}$, y (b) la sustitución hacia
atrás, que ya sabemos que cuesta $n^2$. Contemos (a).

**Pivote $j = 0$.** Hay que modificar los $n-1$ renglones de abajo.
Cada uno necesita:

- una división, para el coeficiente $A_{i0}/A_{00}$;
- $n$ multiplicaciones y $n$ restas, una por cada columna de $A$;
- una multiplicación y una resta para $b_i$.

En total, $(n-1) + 2(n+1)(n-1)$ operaciones.

**Pivote $j = 1$.** Ahora son $n-2$ renglones, cada uno con una
división, $n-1$ multiplicaciones y $n-1$ restas en $A$, y una de cada
una en $\mathbf{b}$: $(n-2) + 2n(n-2)$ operaciones.

**Pivote $j$ cualquiera.** El patrón es

```math
(n-1-j) + 2(n+1-j)(n-1-j).
```

Sumando sobre todos los pivotes, con el cambio de variable
$k = n-1-j$, que va de $n-1$ a $1$:

```math
N = \sum_{j=0}^{n-2}\left[(n-1-j) + 2(n+1-j)(n-1-j)\right]
= \sum_{k=1}^{n-1}\left[k + 2(k+2)k\right]
= \sum_{k=1}^{n-1}\left(2k^2 + 5k\right).
```

Con las sumas conocidas
$\sum_{k=0}^{n-1} k = \frac{(n-1)n}{2}$ y
$\sum_{k=0}^{n-1} k^2 = \frac{(n-1)n(2n-1)}{6}$:

```math
N = 2\,\frac{(n-1)n(2n-1)}{6} + 5\,\frac{(n-1)n}{2}
= \frac{2}{3}n^3 + \frac{3}{2}n^2 - \frac{13}{6}n
\sim \frac{2}{3}n^3.
```

Para $n = 3$ da $18 + 13.5 - 6.5 = 25$: compruébenlo contando las
operaciones del ejemplo de arriba. La sustitución hacia atrás agrega
$n^2$, pero para $n$ grande la eliminación domina por completo:

```math
\frac{2}{3}n^3 + n^2 \sim \frac{2}{3}n^3.
```

Es decir: resolver un sistema general cuesta $`O(n^3)`$, mientras que uno
triangular cuesta $n^2$. Para $n = 1000$, la eliminación es unas 670
veces más cara que la sustitución. `eliminacion_gaussiana_lu.py`
cuenta las operaciones con un contador, como `triangulares.py`, y
coinciden exactamente con la fórmula.

### Un pivote cero

La fórmula del coeficiente divide entre $A_{jj}$. Si algún pivote vale
cero, el método truena, **aunque la matriz no sea singular**. Por
ejemplo,

```math
\begin{pmatrix} 0 & -1 \\ 1 & 1 \end{pmatrix}
\begin{pmatrix} x_0 \\ x_1 \end{pmatrix}
= \begin{pmatrix} 1 \\ 2 \end{pmatrix}
```

tiene determinante $1$ y solución $\mathbf{x} = (3, -1)$, pero
`eliminacion_gaussiana` levanta `ZeroDivisionError` en el primer
paso. La solución es obvia: intercambiar los dos renglones (la segunda
operación elemental, el **pivoteo**). Lo vemos con cuidado en la
sección "Pivoteo", donde también aparece un problema menos visible: el
de los pivotes que no son cero pero sí muy chicos.

## Descomposición LU

La eliminación gaussiana funciona bien, pero tiene un defecto: si
queremos resolver $A\mathbf{x} = \mathbf{b}$ con la **misma** $A$ y
**otro** $\mathbf{b}$, hay que repetir toda la eliminación, que es la
parte cara, $`\sim 2n^3/3`$, aunque $A$ no haya cambiado.

¿Pasa eso en la práctica? Todo el tiempo:

- Para calcular la inversa $A^{-1}$ hay que resolver $n$ sistemas con
  la misma $A$ (abajo lo hacemos).
- Varios métodos para eigenvalores (que veremos en el tema de
  eigenvalores) resuelven un sistema con la misma matriz en cada
  iteración.
- Al resolver ecuaciones diferenciales con métodos implícitos (por
  ejemplo, la ecuación de calor), en cada paso de tiempo se resuelve un
  sistema con la misma matriz y un $\mathbf{b}$ nuevo.

Lo que queremos es **guardar** el resultado de la eliminación para
reutilizarlo. Eso es la descomposición LU.

### La descomposición de Doolittle

Supongamos que una matriz no singular $A$ se puede escribir como el
producto de una triangular inferior $L$ y una triangular superior $U$:

```math
A = LU.
```

Es la **descomposición LU** (o factorización LU) de $A$. (Más adelante,
con el pivoteo, veremos que la historia es un poco más complicada;
por ahora supongamos que se puede.) La descomposición no es única;
para fijarla pedimos que $L$ tenga **unos en la diagonal**, $L_{ii} = 1$.
Eso se llama descomposición de **Doolittle**.

Veamos cómo se construye en el caso general $3\times 3$:

```math
L = \begin{pmatrix} 1 & 0 & 0 \\ L_{10} & 1 & 0 \\ L_{20} & L_{21} & 1 \end{pmatrix},
\qquad
U = \begin{pmatrix} U_{00} & U_{01} & U_{02} \\ 0 & U_{11} & U_{12} \\ 0 & 0 & U_{22} \end{pmatrix}.
```

Multiplicándolas:

```math
A = LU = \begin{pmatrix}
U_{00} & U_{01} & U_{02} \\
L_{10}U_{00} & L_{10}U_{01} + U_{11} & L_{10}U_{02} + U_{12} \\
L_{20}U_{00} & L_{20}U_{01} + L_{21}U_{11} & L_{20}U_{02} + L_{21}U_{12} + U_{22}
\end{pmatrix}.
```

Apliquemos la eliminación gaussiana a **esta** $A$ (sin ningún
$\mathbf{b}$: solo nos interesa la matriz).

**Pivote $j = 0$, $i = 1$.** El coeficiente que anula el primer
elemento del renglón 1 es $L_{10}U_{00}/U_{00} = L_{10}$:
"renglón 1 $-\ L_{10}\times$ renglón 0" deja el renglón 1 como
$(0,\ U_{11},\ U_{12})$.

**Pivote $j = 0$, $i = 2$.** El coeficiente es $L_{20}$, y el renglón 2
queda como $(0,\ L_{21}U_{11},\ L_{21}U_{12} + U_{22})$.

**Pivote $j = 1$, $i = 2$.** El coeficiente es
$L_{21}U_{11}/U_{11} = L_{21}$, y el renglón 2 queda como
$(0,\ 0,\ U_{22})$. Al final:

```math
\begin{pmatrix} U_{00} & U_{01} & U_{02} \\ 0 & U_{11} & U_{12} \\ 0 & 0 & U_{22} \end{pmatrix} = U.
```

Leyendo esto al revés, llegamos a la conclusión importante:

- **Para descomponer $A = LU$ basta con hacer la eliminación
  gaussiana.** Lo que queda de $A$ al final es $U$.
- **Los elementos de $L$ abajo de la diagonal son los coeficientes que
  se usaron en la eliminación:** $L_{ij}$ es el coeficiente con el que
  el pivote $j$ eliminó al renglón $i$.

Con el ejemplo $3\times 3$ de la sección anterior, ya sin
$\mathbf{b}$: $U$ es la matriz triangular a la que llegamos, y $L$
junta los coeficientes $0.5$, $2.5$ y $15$:

```math
U = \begin{pmatrix} 2 & 1 & 1 \\ 0 & 0.5 & -2.5 \\ 0 & 0 & 40 \end{pmatrix},
\qquad
L = \begin{pmatrix} 1 & 0 & 0 \\ 0.5 & 1 & 0 \\ 2.5 & 15 & 1 \end{pmatrix}.
```

No hubo que calcular nada nuevo: solo guardar lo que ya aparecía.
(Multipliquen $L$ por $U$ para convencerse de que dan $A$;
`eliminacion_gaussiana_lu.py` lo hace con la clase `Matrix`.)

Guardar dos matrices $n\times n$ es un desperdicio: sabemos que la
diagonal de $L$ son puros unos y que la mitad de cada matriz son
ceros. En la práctica es común guardar $L$ y $U$ juntas en una sola
matriz, dejando implícitos los unos de la diagonal de $L$. Aquí, por
claridad, las guardamos por separado.

### Resolver un sistema con LU

Con $A = LU$, el sistema $A\mathbf{x} = \mathbf{b}$ se convierte en
$LU\mathbf{x} = \mathbf{b}$, que podemos escribir como
$L(U\mathbf{x}) = \mathbf{b}$. Llamando $\mathbf{y} = U\mathbf{x}$, el
problema se parte en dos sistemas **triangulares**:

```math
\begin{aligned}
L\mathbf{y} &= \mathbf{b} \qquad \text{(sustitución hacia adelante)}, \\
U\mathbf{x} &= \mathbf{y} \qquad \text{(sustitución hacia atrás)}.
\end{aligned}
```

Primero se resuelve el de $L$ para obtener $\mathbf{y}$, y con ese
$\mathbf{y}$ se resuelve el de $U$. Aquí es donde las dos
sustituciones de la sección "Matrices triangulares" se usan juntas.

### Implementación

En [`fiscomp/algebra_lineal.py`](../../fiscomp/algebra_lineal.py):

```python
def descomposicion_lu(A):
    n = len(A)
    U = [renglon[:] for renglon in A]
    L = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    for j in range(n - 1):
        for i in range(j + 1, n):
            coeficiente = U[i][j] / U[j][j]
            for k in range(j, n):
                U[i][k] -= coeficiente * U[j][k]
            L[i][j] = coeficiente
    return L, U


def resolver_lu(L, U, b):
    y = sustitucion_adelante(L, b)
    return sustitucion_atras(U, y)
```

`descomposicion_lu` es idéntica a `eliminacion_gaussiana`, salvo que:

- no hay $\mathbf{b}$: solo descomponemos la matriz;
- a la matriz que se va modificando la llamamos `U` en vez de `A`;
- cada coeficiente se guarda en `L[i][j]`, y `L` empieza como la
  identidad para que su diagonal ya tenga los unos.

**Una diferencia con el libro.** En el libro, `lusolve(A, bs)` recibe
la matriz $A$ y llama a `ludec` adentro, así que cada vez que se
resuelve un sistema se vuelve a descomponer $A$: justo lo que queríamos
evitar. Aquí lo separamos en dos funciones: `descomposicion_lu(A)` se
llama **una vez**, y `resolver_lu(L, U, b)` se llama con cada
$\mathbf{b}$ distinto:

```python
L, U = descomposicion_lu(A)        # ~2n^3/3, una sola vez
for b in lados_derechos:
    x = resolver_lu(L, U, b)       # 2n^2 cada vez
```

Con la matriz de prueba del libro, la solución por LU no es idéntica
bit a bit a la de eliminación gaussiana: coinciden en unos 13 dígitos.
Matemáticamente hacen las mismas operaciones con $\mathbf{b}$, pero en
distinto orden (la eliminación va restando de $b_i$ un término a la
vez; la sustitución hacia adelante junta primero la suma y luego la
resta), y en punto flotante el orden cambia el redondeo (unidad 06).

### Conteo de operaciones

La descomposición es la eliminación sin $\mathbf{b}$: para el pivote
$j$ se ahorran la multiplicación y la resta de $b_i$ en cada uno de los
$n-1-j$ renglones. Repitiendo la cuenta de la sección anterior con
ese cambio se obtiene

```math
N_{LU} = \frac{2}{3}n^3 + \frac{1}{2}n^2 - \frac{7}{6}n \sim \frac{2}{3}n^3.
```

(Háganlo: es el mismo procedimiento, con $2k^2 + 3k$ en lugar de
$2k^2 + 5k$. `eliminacion_gaussiana_lu.py` comprueba la fórmula
contando.) El término dominante es el mismo que en la eliminación
gaussiana. La diferencia está en lo que sigue:

| Costo aproximado | Un sistema | $m$ sistemas con la misma $A$ |
|---|---|---|
| Eliminación gaussiana | $\frac{2}{3}n^3 + n^2$ | $`m\,\frac{2}{3}n^3`$ |
| LU | $\frac{2}{3}n^3 + 2n^2$ | $\frac{2}{3}n^3 + 2mn^2$ |

Para un solo sistema da casi lo mismo (LU hace una sustitución más).
Para muchos, LU gana por mucho: con $n = 100$ y $m = 20$, el conteo
predice que LU es unas 12 veces más rápida, y el script mide alrededor
de 14.

Los flops no son lo único que importa: también la memoria. Tal como lo
implementamos, LU guarda dos matrices $n\times n$, y la eliminación
gaussiana solo una.

### La matriz inversa

En la práctica casi nunca hace falta invertir una matriz: para
resolver $A\mathbf{x} = \mathbf{b}$, calcular $A^{-1}$ y luego
$A^{-1}\mathbf{b}$ cuesta más y acumula más error que usar LU
directamente. Pero hay usos legítimos; el primero para nosotros es el
número de condición, $`\kappa(A) = \|A\|\,\|A^{-1}\|`$.

La inversa cumple $AA^{-1} = I$. El truco es ver la identidad como $n$
vectores columna $\mathbf{e}_i$, cada uno con un solo $1$ en el lugar
$i$, y a $A^{-1}$ también como $n$ columnas $\mathbf{x}_i$:

```math
A^{-1} = \begin{pmatrix} \mathbf{x}_0 & \mathbf{x}_1 & \cdots & \mathbf{x}_{n-1} \end{pmatrix},
\qquad
I = \begin{pmatrix} \mathbf{e}_0 & \mathbf{e}_1 & \cdots & \mathbf{e}_{n-1} \end{pmatrix}.
```

Así, en vez de atacar $AA^{-1} = I$ de un golpe, lo partimos en $n$
sistemas, uno por columna:

```math
A\mathbf{x}_i = \mathbf{e}_i, \qquad i = 0, 1, \ldots, n-1.
```

Es nuestro primer ejemplo de muchos sistemas con la misma $A$: se
descompone $A = LU$ **una vez**, y para cada columna se resuelve
$L(U\mathbf{x}_i) = \mathbf{e}_i$ con una sustitución hacia adelante y
una hacia atrás. Cada par cuesta $2n^2$ y hay $n$ columnas:

```math
\underbrace{\frac{2}{3}n^3}_{\text{LU}} + \underbrace{n\cdot 2n^2}_{\text{sustituciones}} = \frac{8}{3}n^3 \text{ operaciones}.
```

Con eliminación gaussiana tendríamos que repetir todo, $`O(n^3)`$, para
cada una de las $n$ columnas: $`O(n^4)`$, demasiado caro.

En [`fiscomp/algebra_lineal.py`](../../fiscomp/algebra_lineal.py):

```python
def inversa(A):
    n = len(A)
    L, U = descomposicion_lu(A)
    columnas = []
    for k in range(n):
        e_k = [1.0 if i == k else 0.0 for i in range(n)]
        columnas.append(resolver_lu(L, U, e_k))
    # columnas[k][i] es el elemento (i, k) de la inversa: transponemos.
    return [[columnas[k][i] for k in range(n)] for i in range(n)]
```

Ojo con el último paso: `columnas` guarda las columnas de $A^{-1}$ como
si fueran renglones, así que hay que transponer.

Con $A^{-1}$ podemos calcular el número de condición de cualquier
matriz, no solo de las $2\times 2$ de `analisis_de_error.py`. Para la
matriz de prueba del libro, `eliminacion_gaussiana_lu.py` obtiene
$`\kappa_\infty(A) \approx 3.9\times 10^{8}`$: la explicación de los
8 dígitos perdidos de la sección "Eliminación gaussiana".

### El determinante

Para una matriz triangular, el determinante es el producto de la
diagonal. Y el determinante de un producto es el producto de los
determinantes. Con $A = LU$:

```math
\det(A) = \det(L)\,\det(U)
= \left(\prod_{i=0}^{n-1} 1\right)\left(\prod_{i=0}^{n-1} U_{ii}\right)
= \prod_{i=0}^{n-1} U_{ii},
```

porque la diagonal de $L$ son puros unos. **El determinante es el
producto de la diagonal de** $U$, que ya teníamos al terminar la
eliminación. De paso, esto dice que $A$ es no singular si y solo si
ningún $U_{ii}$ vale cero.

Cuesta lo mismo que la descomposición, $`\sim 2n^3/3`$. La fórmula de
cofactores que se ve en los cursos de álgebra cuesta del orden de
$n!$: para $n = 20$, unas $10^{18}$ operaciones contra unas
$5\times 10^{3}$.

La función `determinante(A)` de `fiscomp/algebra_lineal.py` hace
exactamente esto. Para la matriz de prueba del libro da
$`\det(A) \approx -7\times 10^{-12}`$. Recordemos la sección "¿Determinante
chico?": eso no dice nada por sí solo. Lo que importa es
$`\kappa \approx 4\times 10^{8}`$. El proyecto de espines
([`simulaciones/espines/`](../../simulaciones/espines/)) usa esta
función para buscar eigenvalores como raíces de $\det(H - \lambda I)$.

## Pivoteo

Hasta aquí aplicamos la eliminación gaussiana (y su prima, la LU) a
casos sin complicaciones. Ahora veremos que es fácil que algo salga
mal, cómo se arregla lo más grave, y qué remedios más complicados
existen.

### Inestabilidad sin mal condicionamiento

En la sección "Análisis de error" vimos que un problema puede estar
**mal condicionado**: es una propiedad del *problema*, que se mide con
$\kappa(A)$. Ahora veremos lo contrario: problemas perfectamente bien
condicionados en los que el *método* falla. La culpa, en ese caso, es
del método. Veamos tres ejemplos.

**Ejemplo 1: un cero en el primer pivote.**

```math
\left(\begin{array}{cc|c}
0 & -1 & 1 \\
1 & 1 & 2
\end{array}\right)
```

Está bien condicionado, con $`\kappa_\infty(A) = 4`$, y su solución
es fácil de obtener a mano: $\mathbf{x} = (3, -1)$. Pero la
eliminación gaussiana empieza (y termina) con el pivote $j = 0$, y el
coeficiente $A_{10}/A_{00}$ divide entre $A_{00} = 0$. En el libro,
con NumPy, el resultado es `[nan nan]` (*not a number*); con nuestras
listas, Python levanta `ZeroDivisionError`.

**Ejemplo 2: un cero escondido.** Si en el ejemplo $3\times 3$ de la
sección "Eliminación gaussiana" cambiamos un par de elementos (sin
cambiar la solución, que sigue siendo $`\mathbf{x} = (4, -2, 2)`$):

```math
\left(\begin{array}{ccc|c}
2 & 1 & 1 & 8 \\
2 & 1 & -4 & -2 \\
5 & 10 & 5 & 10
\end{array}\right).
```

El primer pivote, $2$, no tiene nada de malo. Pero después de eliminar
la columna $0$:

```math
\left(\begin{array}{ccc|c}
2 & 1 & 1 & 8 \\
0 & 0 & -5 & -10 \\
0 & 7.5 & 2.5 & -10
\end{array}\right),
```

y el siguiente pivote, $A_{11}$, vale cero. El problema no se veía a
simple vista: apareció a media eliminación.

**Ejemplo 3: un pivote diminuto.** Cuando un método falla en un caso,
casi siempre falla también en los casos parecidos. Cambiemos el $0$ del
ejemplo 1 por un número muy chico:

```math
\left(\begin{array}{cc|c}
10^{-20} & -1 & 1 \\
1 & 1 & 2
\end{array}\right),
\qquad
\mathbf{x} = \left(\frac{3}{1 + 10^{-20}},\ \frac{-1 + 2\times 10^{-20}}{1 + 10^{-20}}\right).
```

Sigue bien condicionado, con $`\kappa_\infty(A) \approx 4`$, y como
$10^{-20}$ es mucho menor que $\epsilon_{\text{mach}}$, esperaríamos
obtener $(3, -1)$, igual que en el ejemplo 1. Ahora el pivote no es
cero, así que la eliminación procede, con coeficiente
$1/10^{-20} = 10^{20}$:

```math
\left(\begin{array}{cc|c}
10^{-20} & -1 & 1 \\
0 & 1 + 10^{20} & 2 - 10^{20}
\end{array}\right).
```

Es la situación de la unidad 06: en punto flotante,
$1 + 10^{20} = 10^{20}$ y $2 - 10^{20} = -10^{20}$; el $1$ y el $2$ se
pierden por completo. Nadie divide entre cero, así que la sustitución
hacia atrás procede normalmente, pero con la matriz

```math
\left(\begin{array}{cc|c}
10^{-20} & -1 & 1 \\
0 & 10^{20} & -10^{20}
\end{array}\right),
```

que da $\mathbf{x} = (0, -1)$. Completamente distinto de la respuesta
correcta, y sin ningún aviso.

### Pivoteo parcial

En el ejemplo 1, si hubiéramos escrito las ecuaciones en el otro
orden,

```math
\left(\begin{array}{cc|c}
1 & 1 & 2 \\
0 & -1 & 1
\end{array}\right),
```

no habría pasado nada: la matriz ya es triangular superior. Pero
reacomodar a mano los renglones que "se ven mal" no es un método.
Necesitamos una regla general.

En los tres ejemplos, el problema fue usar como pivote un número muy
chico (o cero). Al eliminar, el coeficiente es $A_{ij}/A_{jj}$: si
$A_{jj}$ es chico, el coeficiente es enorme. La regla del **pivoteo
parcial** es: justo antes de usar el renglón $j$ como pivote, buscar en
la **columna** $j$, del renglón $j$ para abajo, el elemento de mayor
valor absoluto. Si está en el renglón $k$,

```math
|A_{kj}| = \max_{j \le m \le n-1} |A_{mj}|,
```

se **intercambian los renglones** $j$ y $k$ (de $A$ y de $\mathbf{b}$)
y se elimina como siempre. (Si hay empate, se toma el primero, el $k$
más chico.) Como ahora el pivote es el más grande de su columna, los
coeficientes nunca pasan de 1 en valor absoluto.

Ojo con los índices: el segundo índice, la columna $j$, está fijo; lo
que recorremos son los renglones $m$, de $j$ para abajo.

Intercambiar renglones es la segunda operación elemental de la
sección "Sistemas de ecuaciones lineales". No cambia la solución,
pero sí el signo del determinante. Se llama pivoteo *parcial* porque
solo intercambiamos renglones; también se podrían intercambiar
columnas. (Una advertencia de vocabulario: "pivote" es el renglón o el
elemento con el que eliminamos, y "pivotear" es intercambiar
renglones. Están relacionados, pero no son lo mismo.)

### Implementación

En [`fiscomp/algebra_lineal.py`](../../fiscomp/algebra_lineal.py),
`eliminacion_gaussiana_pivoteo` es `eliminacion_gaussiana` con un
bloque nuevo al inicio de cada pivote:

```python
def eliminacion_gaussiana_pivoteo(A, b):
    A = [renglon[:] for renglon in A]
    b = b[:]
    n = len(b)
    for j in range(n - 1):
        k = j
        for m in range(j + 1, n):
            if abs(A[m][j]) > abs(A[k][j]):
                k = m
        if k != j:
            A[j], A[k] = A[k], A[j]
            b[j], b[k] = b[k], b[j]
        for i in range(j + 1, n):
            coeficiente = A[i][j] / A[j][j]
            for k_col in range(j, n):
                A[i][k_col] -= coeficiente * A[j][k_col]
            b[i] -= coeficiente * b[j]
    return sustitucion_atras(A, b)
```

Algunos detalles:

- Para encontrar $k$, el libro usa `numpy.argmax`. Nosotros recorremos
  la columna con un ciclo. Como la comparación es `>` estricta, en un
  empate se queda el primero.
- `if k != j`: intercambiar un renglón consigo mismo no hace nada.
- El intercambio usa la asignación múltiple de Python,
  `A[j], A[k] = A[k], A[j]`, sin variable temporal. Con listas de
  listas esto intercambia los renglones completos de un golpe: lo que
  se intercambia son las referencias a las dos sublistas. (El libro
  tiene que escribir `A[j,:].copy()`, porque con NumPy un pedazo de un
  arreglo es una *vista* del original, y sin la copia se sobreescribe
  en vez de intercambiar. Con listas no hay ese problema.)
- El ciclo interno usa `k_col` en vez de `k` para no confundirlo con el
  renglón `k` del pivoteo.

`pivoteo.py` resuelve los tres ejemplos: sin pivoteo, los dos primeros
levantan `ZeroDivisionError` y el tercero da $(0, -1)$; con pivoteo,
los tres dan la respuesta correcta. Con la matriz de prueba del libro,
en cambio, los dos métodos pierden unos 8 dígitos. Ahí el problema es
$\kappa(A)$, y eso no lo arregla ningún método. En el libro, el
pivoteo acerca el resultado al de `numpy.linalg.solve`, que también
pivotea; con nuestra matriz de $4\times 4$ las diferencias son solo de
redondeo.

**Para pensar.**

1. Agreguen pivoteo parcial a `descomposicion_lu` (problema 4.16 del
   libro). La LU no modifica $\mathbf{b}$, así que hay que guardar qué
   renglones se intercambiaron, para luego intercambiar igual los
   elementos de $\mathbf{b}$ antes de las sustituciones.
2. Con pivoteo, `determinante` tiene que cambiar el signo en cada
   intercambio. ¿Cómo lo modificarían?
3. En vez de intercambiar físicamente los renglones, se puede dejar la
   matriz como está y llevar la cuenta del orden en que se usaron los
   renglones como pivote. Al final queda una triangular "revuelta".
   Es más eficiente, pero requiere más contabilidad.

### Más allá del pivoteo parcial

El pivoteo parcial no siempre basta. Tomen el ejemplo 3 y multipliquen
el segundo renglón por $10^{-20}$:

```math
\left(\begin{array}{cc|c}
10^{-20} & -1 & 1 \\
10^{-20} & 10^{-20} & 2\times 10^{-20}
\end{array}\right).
```

Son exactamente las mismas ecuaciones, así que la respuesta debería
ser la misma. Pero ahora los dos candidatos a pivote tienen la misma
magnitud, no hay intercambio, y volvemos a obtener $(0, -1)$. El libro
reporta que `numpy.linalg.solve` también falla aquí.

El remedio es el **pivoteo parcial escalado**: elegir el pivote con el
mayor tamaño *relativo* a su propio renglón. Para cada renglón $i$ se
define un factor de escala, su elemento más grande,

```math
s_i = \max_{0 \le j \le n-1} |A_{ij}|,
```

y se busca el $k$ más chico tal que

```math
\frac{|A_{kj}|}{s_k} = \max_{j \le m \le n-1} \frac{|A_{mj}|}{s_m}.
```

En el ejemplo, el primer renglón tiene $s_0 = 1$ y razón $10^{-20}$, y
el segundo $s_1 = 10^{-20}$ y razón $1$: sí se intercambian.
Implementarlo es el problema 4.15 del libro.

Ni siquiera el pivoteo escalado garantiza estabilidad. Lo que sí la
garantiza es el **pivoteo completo**, que intercambia renglones y
columnas, pero es más caro, y en la práctica el pivoteo parcial
(escalado o no) casi nunca falla. Por eso la mayoría de las
bibliotecas no lo implementan.

**Cuándo no hace falta pivotear.** La eliminación gaussiana funciona
sin pivoteo para matrices simétricas positivas definidas y para
matrices **diagonalmente dominantes**: aquellas en las que, en cada
renglón, el elemento de la diagonal es al menos tan grande, en valor
absoluto, como la suma de los demás,

```math
|A_{ii}| \ge \sum_{j \ne i} |A_{ij}|, \qquad i = 0, 1, \ldots, n-1.
```

Por ejemplo, el ejemplo $3\times 3$ de la sección "Eliminación
gaussiana" con el segundo y tercer renglón intercambiados:

```math
\left(\begin{array}{ccc|c}
2 & 1 & 1 & 8 \\
5 & 10 & 5 & 10 \\
1 & 1 & -2 & -2
\end{array}\right).
```

Con ella, `eliminacion_gaussiana_pivoteo` no intercambia nada y da
exactamente el mismo resultado que `eliminacion_gaussiana`.

¿Por qué querríamos evitar el pivoteo? Porque buscar el máximo e
intercambiar renglones cuesta, y porque el intercambio destruye
estructuras útiles: una matriz simétrica o de banda (con ceros lejos
de la diagonal) deja de serlo, y los métodos especializados para esas
estructuras ya no se pueden usar.

## Método iterativo de Jacobi

Todos los métodos que hemos visto son **directos**: hacen un número
fijo de operaciones, conocido de antemano. Los métodos **iterativos**
(o *de relajación*) parten de una propuesta de solución y la mejoran
una y otra vez hasta que deja de cambiar. Cuántas iteraciones hacen
falta depende de la matriz, del método, de la propuesta inicial y del
criterio de convergencia; y a veces no convergen nunca.

Por eso, en general, son más lentos que los directos. Su ventaja
aparece con matrices **ralas**, con casi todos sus elementos iguales a
cero, que son justo las que aparecen al resolver ecuaciones
diferenciales parciales (capítulo 8). Además, se *autocorrigen*: cada
iteración no arrastra el error de redondeo de la anterior. Aquí vemos
el más sencillo, el **método de Jacobi**; es también un primer
contacto con la idea de iterar hasta converger, que reaparecerá en los
eigenvalores y en la búsqueda de raíces.

### El algoritmo

Escribamos $A\mathbf{x} = \mathbf{b}$ con índices y despejemos la
componente $x_i$, la que multiplica al elemento de la diagonal $A_{ii}$:

```math
x_i = \frac{1}{A_{ii}}\left(b_i - \sum_{j=0}^{i-1} A_{ij}x_j - \sum_{j=i+1}^{n-1} A_{ij}x_j\right),
\qquad i = 0, 1, \ldots, n-1.
```

Se parece a las fórmulas de las sustituciones hacia adelante y hacia
atrás, pero hay una diferencia esencial: allá, al calcular $x_i$, ya
conocíamos los $x_j$ que aparecen a la derecha. Aquí no conocemos
ninguno.

Jacobi lo resuelve así: se propone un vector inicial
$\mathbf{x}^{(0)}$ (el superíndice entre paréntesis es el número de
iteración), se mete del lado derecho, y sale una propuesta mejorada
$\mathbf{x}^{(1)}$. Esa se vuelve a meter, y así sucesivamente. En la
iteración $k$:

```math
x_i^{(k)} = \frac{1}{A_{ii}}\left(b_i - \sum_{j=0}^{i-1} A_{ij}x_j^{(k-1)} - \sum_{j=i+1}^{n-1} A_{ij}x_j^{(k-1)}\right),
\qquad i = 0, 1, \ldots, n-1.
```

Todas las componentes nuevas se calculan con el vector **anterior**
completo.

**¿Cuándo parar?** No conocemos la solución exacta, así que no podemos
calcular el error. Lo que sí podemos medir es cuánto cambia cada
componente de una iteración a la siguiente. Como en capítulos
anteriores, usamos un cambio *relativo*, sumado sobre las componentes:
paramos cuando

```math
\sum_{i=0}^{n-1}\left|\frac{x_i^{(k)} - x_i^{(k-1)}}{x_i^{(k)}}\right| \le \epsilon,
```

con $\epsilon$ una tolerancia que escogemos. Dividimos entre
$x_i^{(k)}$ y no entre $x_i^{(k-1)}$ porque es nuestra mejor
estimación hasta el momento, y porque es común empezar con
$\mathbf{x}^{(0)} = \mathbf{0}$, que haría dividir entre cero en la
primera iteración. (Hay otras opciones razonables, como el cambio
máximo o la suma en cuadratura. Un criterio más sistemático usa el
residuo de la sección "Análisis de error":
$`\|\mathbf{b} - A\mathbf{x}^{(k)}\| \le \epsilon\left(\|A\|\,\|\mathbf{x}^{(k)}\| + \|\mathbf{b}\|\right)`$,
que mide directamente qué tan bien se cumple $A\mathbf{x} = \mathbf{b}$.)

**¿Converge?** No siempre, aunque la propuesta inicial sea buena. Una
condición **suficiente** es que $A$ sea diagonalmente dominante: así
converge para cualquier $\mathbf{x}^{(0)}$ (problema 4.20 del libro).
Un sistema que no lo es a veces se vuelve dominante reacomodando las
ecuaciones; y algunos sistemas que no lo son convergen de todos modos
con ciertas propuestas iniciales.

### Implementación

En [`fiscomp/algebra_lineal.py`](../../fiscomp/algebra_lineal.py), una
iteración, el criterio de paro y el método completo van en tres
funciones:

```python
def paso_jacobi(A, b, x):
    n = len(b)
    x_nuevo = []
    for i in range(n):
        suma = 0.0
        for j in range(n):
            if j != i:
                suma += A[i][j] * x[j]
        x_nuevo.append((b[i] - suma) / A[i][i])
    return x_nuevo


def cambio_relativo(x_viejo, x_nuevo):
    return sum(abs((xn - xv) / xn) for xv, xn in zip(x_viejo, x_nuevo))


def jacobi(A, b, kmax=50, tol=1e-6):
    x = [0.0] * len(b)
    for k in range(1, kmax):
        x_nuevo = paso_jacobi(A, b, x)
        error = cambio_relativo(x, x_nuevo)
        x = x_nuevo
        if error < tol:
            break
    else:
        return None, kmax
    return x, k
```

Algunos detalles:

- `paso_jacobi` junta las dos sumas de la fórmula en una sola,
  saltándose $j = i$. El libro las escribe por separado con NumPy.
- `paso_jacobi` construye una lista **nueva** y no toca `x`. Si fuéramos
  modificando `x` en su lugar, al calcular $x_i$ ya estaríamos usando
  algunos valores nuevos (eso es otro método, Gauss-Seidel). Y si
  escribiéramos `x_viejo = x` creyendo que es una copia, los dos nombres
  apuntarían a la misma lista, el cambio relativo daría siempre cero,
  y el método "convergería" en la primera iteración. El libro insiste
  en este punto: hay que usar `np.copy`.
- `jacobi` usa una construcción de Python que no habíamos visto: un
  `for` con `else`. El bloque `else` de un ciclo se ejecuta solo si
  el ciclo termina **sin** pasar por un `break`. Aquí eso quiere decir
  que se acabaron las iteraciones sin converger, y en ese caso la
  función regresa `(None, kmax)`. Si converge, el `break` se salta el
  `else` y la función regresa la solución y el número de iteraciones.

`jacobi.py` lo prueba con la matriz de prueba del libro más $21I$ (que
la vuelve diagonalmente dominante): converge en 33 iteraciones, "unas
tres docenas" como dice el libro, y coincide con la eliminación
gaussiana en unos 7 dígitos, lo que permite la tolerancia de
$10^{-6}$. Con la matriz original, que no es dominante, las
componentes crecen en cada iteración y `jacobi` regresa `None`.

**Costo.** Cada iteración cuesta unas $2n^2$ operaciones, contra
$`\sim 2n^3/3`$ de la eliminación gaussiana. Si el número de
iteraciones no crece mucho con $n$, Jacobi gana para $n$ grande.
`jacobi.py` lo mide con una matriz tridiagonal (diagonal $4$, vecinas
$-1$): siempre unas 40 iteraciones, y para $n = 400$ Jacobi ya es
unas 3 veces más rápido. Si además no guardáramos los ceros, cada
iteración costaría solo unas $5n$ operaciones.

**Para pensar.** El método de **Gauss-Seidel** cambia un solo detalle:
al calcular $x_i^{(k)}$, usa los valores *nuevos*
$x_0^{(k)}, \ldots, x_{i-1}^{(k)}$, que ya se calcularon en esta misma
iteración, en vez de los viejos. Impleméntenlo (problema 4.19 del
libro) y comparen cuántas iteraciones necesita con la matriz de prueba
más $21I$.

## El problema de eigenvalores, numéricamente

Pasamos a la "segunda mitad" del álgebra lineal:

```math
A\mathbf{v}_i = \lambda_i\mathbf{v}_i .
```

El truco principal de las secciones anteriores ya no sirve: restarle
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
Jacobi):

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
ideal para LU (sección "Resolver un sistema con LU"). Se factoriza
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
  singular, como ya nos decía $`\kappa \approx 4\times 10^{8}`$.)
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
- [`fiscomp/algebra_lineal.py`](../../fiscomp/algebra_lineal.py):
  producto matriz-vector (`mat_vec`), residuo, normas de vectores y de
  matrices, partes triangulares de una matriz, las sustituciones
  hacia adelante y hacia atrás, eliminación gaussiana sin y con
  pivoteo (`eliminacion_gaussiana`, `eliminacion_gaussiana_pivoteo`),
  descomposición LU (`descomposicion_lu`, `resolver_lu`), la inversa y
  el determinante (`inversa`, `determinante`), el método de Jacobi
  (`paso_jacobi`, `cambio_relativo`, `jacobi`), los métodos de la
  potencia para eigenvalores (`producto_punto`, `normalizar`,
  `cociente_rayleigh`, `potencia`, `potencia_inversa`), el método QR
  (`transpuesta`, `mat_mat`, `descomposicion_qr`, `metodo_qr`,
  `eigenpares`) y la matriz de prueba del libro (`crear_prueba`).
  Todo con listas de listas, sin NumPy. Aquí se irá acumulando el
  resto del tema.
- [`analisis_de_error.py`](analisis_de_error.py): el ejemplo de Kahan
  completo (sección "Análisis de error"): residuo diminuto con una
  solución totalmente equivocada; el cambio drástico de la solución al
  perturbar un solo elemento; por qué el determinante no sirve como
  criterio (con $`D = 0.1\,I`$); el número de condición $\kappa(A)$ con
  las normas infinito y de Frobenius, y la comprobación de la cota
  $`\|\Delta\mathbf{x}\|/\|\mathbf{x}\| \le \kappa(A)\,\|\Delta A\|/\|A\|`$.
  Se escribió antes de tener un método general, así que todo es
  $2\times 2$ y se resuelve con la regla de Cramer; con LU, el número
  de condición ya se calcula para cualquier $n$ (ver
  `eliminacion_gaussiana_lu.py`).
- [`condicion_eigenvalores.py`](condicion_eigenvalores.py): el número
  de condición $\kappa^{ev}$ de los eigenvalores de una matriz no
  simétrica contra una simétrica, comparado con el cambio real al
  perturbar la matriz; y cómo gira un eigenvector cuando dos
  eigenvalores están cerca. Otra vez $2\times 2$, con el polinomio
  característico.
- [`triangulares.py`](triangulares.py): sustitución hacia adelante y
  hacia atrás con la matriz de prueba del libro, comprobación de las
  soluciones, conteo de operaciones y medición de tiempos.
- [`eliminacion_gaussiana_lu.py`](eliminacion_gaussiana_lu.py): el
  ejemplo $3\times 3$ paso a paso (imprime la matriz aumentada después
  de cada operación), su $L$ y su $U$ con la comprobación $LU = A$, el
  sistema de prueba del libro resuelto con los dos métodos, la
  inversa, $\kappa(A)$ y el determinante con LU, el conteo de
  operaciones contra las fórmulas, el tiempo de resolver muchos
  sistemas con la misma $A$, y qué pasa con un pivote cero.
- [`pivoteo.py`](pivoteo.py): los tres ejemplos del libro de
  inestabilidad sin mal condicionamiento (un cero en el primer pivote,
  uno escondido y un pivote diminuto), resueltos sin y con pivoteo
  parcial; la matriz de prueba con los dos métodos; el caso en el que
  el pivoteo parcial no basta; y una matriz diagonalmente dominante,
  con la que el pivoteo no intercambia nada.
- [`jacobi.py`](jacobi.py): el método de Jacobi iteración por
  iteración con la matriz de prueba más $21I$, comparado con
  eliminación gaussiana; lo que pasa sin dominancia diagonal; y el
  tiempo contra la eliminación gaussiana con una matriz tridiagonal.
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
ejemplo `python3 unidades/09_algebra_lineal/triangulares.py`.
