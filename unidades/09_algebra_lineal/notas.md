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
  inversa y el número de condición).
- Corresponde a las secciones 4.1 a 4.3.3 del libro de Gezerlis.

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

$$
\begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix}
\begin{pmatrix} x \\ y \end{pmatrix}
= \begin{pmatrix} x' \\ y' \end{pmatrix}.
$$

Si conocemos $\mathbf{r}'$ y queremos $\mathbf{r}$, hay que resolver
un sistema de dos ecuaciones lineales. Nótese que la matriz de
rotación **no es simétrica** (sus dos elementos fuera de la diagonal
no son iguales). Las matrices simétricas aparecen por todos lados en
física, pero no siempre; por eso estudiaremos matrices generales.

**2. Potenciales electrostáticos.** Tenemos $n$ cargas $q_j$
*desconocidas* en posiciones $\mathbf{R}_j$ conocidas, y medimos el
potencial $\phi(\mathbf{r}_i)$ en $n$ puntos $\mathbf{r}_i$ conocidos.
Por el principio de superposición,

$$
\phi(\mathbf{r}_i) = \sum_{j=0}^{n-1} \frac{k}{|\mathbf{r}_i - \mathbf{R}_j|}\, q_j,
\qquad i = 0, 1, \ldots, n-1.
$$

Es el problema **inverso** del proyecto
[`simulacion_potencial`](../../simulaciones/simulacion_potencial/): ahí
conocíamos las cargas y calculábamos el potencial; aquí conocemos el
potencial y queremos las cargas. Para cuatro cargas es el sistema
$4\times 4$

$$
\begin{pmatrix}
\frac{k}{|\mathbf{r}_0-\mathbf{R}_0|} & \frac{k}{|\mathbf{r}_0-\mathbf{R}_1|} & \frac{k}{|\mathbf{r}_0-\mathbf{R}_2|} & \frac{k}{|\mathbf{r}_0-\mathbf{R}_3|} \\
\frac{k}{|\mathbf{r}_1-\mathbf{R}_0|} & \frac{k}{|\mathbf{r}_1-\mathbf{R}_1|} & \frac{k}{|\mathbf{r}_1-\mathbf{R}_2|} & \frac{k}{|\mathbf{r}_1-\mathbf{R}_3|} \\
\frac{k}{|\mathbf{r}_2-\mathbf{R}_0|} & \frac{k}{|\mathbf{r}_2-\mathbf{R}_1|} & \frac{k}{|\mathbf{r}_2-\mathbf{R}_2|} & \frac{k}{|\mathbf{r}_2-\mathbf{R}_3|} \\
\frac{k}{|\mathbf{r}_3-\mathbf{R}_0|} & \frac{k}{|\mathbf{r}_3-\mathbf{R}_1|} & \frac{k}{|\mathbf{r}_3-\mathbf{R}_2|} & \frac{k}{|\mathbf{r}_3-\mathbf{R}_3|}
\end{pmatrix}
\begin{pmatrix} q_0 \\ q_1 \\ q_2 \\ q_3 \end{pmatrix}
= \begin{pmatrix} \phi(\mathbf{r}_0) \\ \phi(\mathbf{r}_1) \\ \phi(\mathbf{r}_2) \\ \phi(\mathbf{r}_3) \end{pmatrix},
$$

cuya matriz tampoco es simétrica en general.

**3. Momentos principales de inercia.** Para un cuerpo rígido que rota
alrededor de un eje arbitrario, el tensor de inercia es

$$
I_{\alpha\beta} = \int \rho(\mathbf{r})\left(\delta_{\alpha\beta}\, r^2 - r_\alpha r_\beta\right) d^3r,
\qquad
I = \begin{pmatrix}
I_{xx} & I_{xy} & I_{xz} \\
I_{yx} & I_{yy} & I_{yz} \\
I_{zx} & I_{zy} & I_{zz}
\end{pmatrix},
$$

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

$$
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
$$

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

$$
A\mathbf{v} = \lambda\mathbf{v}.
$$

La diferencia crucial con $A\mathbf{x} = \mathbf{b}$ es que aquí
**tanto** el número $\lambda$ (el **eigenvalor**) **como** el vector
$\mathbf{v}$ (el **eigenvector**) son incógnitas. Pasando todo al lado
izquierdo,

$$
(A - \lambda I)\mathbf{v} = \mathbf{0},
$$

con $I$ la matriz identidad $n\times n$. Es un sistema lineal con
coeficientes $A - \lambda I$ y lado derecho cero, pero con $n+1$
incógnitas ($\lambda, v_0, \ldots, v_{n-1}$), así que no podemos
esperar una solución única. La solución trivial $\mathbf{v} = \mathbf{0}$
siempre existe; para que haya una no trivial, la matriz
$A - \lambda I$ tiene que ser **singular**:

$$
|A - \lambda I| = 0.
$$

Desarrollando el determinante se obtiene la **ecuación
característica**, un polinomio de grado $n$ en $\lambda$:

$$
(-1)^n\lambda^n + c_{n-1}\lambda^{n-1} + \cdots + c_1\lambda + c_0 = 0.
$$

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
$\|\mathbf{v}_i\| = 1$).

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

$$
(A|\mathbf{b}) = \left(\begin{array}{cc|c}
1.2969 & 0.8648 & 0.8642 \\
0.2161 & 0.1441 & 0.1440
\end{array}\right),
$$

y supongamos que alguien nos da la solución aproximada
$\tilde{\mathbf{x}}^T = (0.9911, -0.4870)$. (Escribimos la transpuesta
para ahorrar espacio.)

**A posteriori.** Una forma de medir qué tan buena es: el **vector
residuo**

$$
\mathbf{r} = \mathbf{b} - A\tilde{\mathbf{x}}.
$$

Si $\mathbf{x}$ es la solución exacta, $A\mathbf{x} = \mathbf{b}$ y el
residuo es cero; una buena aproximación debería dar un residuo chico.
Aquí $\mathbf{r}^T = (10^{-8}, -10^{-8})$: diminuto. Uno concluiría
que $\tilde{\mathbf{x}}$ es excelente. Pero la solución exacta es

$$
\mathbf{x}^T = (2, -2),
$$

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
escribe con doble barra, $\|A\|$ (la barra sencilla $|A|$ es el
determinante, o el valor absoluto de un número). Usaremos dos:

**Norma de Frobenius:**

$$
\|A\|_F = \sqrt{\sum_{i=0}^{n-1}\sum_{j=0}^{n-1} |A_{ij}|^2}.
$$

**Norma infinito** (máxima suma de renglón):

$$
\|A\|_\infty = \max_{0\le i\le n-1} \sum_{j=0}^{n-1} |A_{ij}|.
$$

Cualquier norma de matrices cuadradas cumple:

$$
\begin{aligned}
&\|A\| \ge 0, \\
&\|A\| = 0 \text{ si y solo si todos los } A_{ij} = 0, \\
&\|kA\| = |k|\,\|A\|, \\
&\|A + B\| \le \|A\| + \|B\| \quad\text{(desigualdad del triángulo)}, \\
&\|AB\| \le \|A\|\,\|B\|.
\end{aligned}
$$

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
$D\mathbf{x} = \mathbf{b}$ se resuelve trivialmente ($\mathbf{x} = 10\,\mathbf{b}$)
y cambiar un poco $D$ o $\mathbf{b}$ cambia un poco
$\mathbf{x}$. El determinante es sensible a la escala y a la dimensión
de la matriz de una forma que no tiene nada que ver con el
condicionamiento. Necesitamos otra medida.

### Número de condición para sistemas lineales

Partimos del problema sin perturbar,

$$
A\mathbf{x} = \mathbf{b},
$$

y cambiamos un poco $A$ (dejando $\mathbf{b}$ fija). La solución
cambia también:

$$
(A + \Delta A)(\mathbf{x} + \Delta\mathbf{x}) = \mathbf{b}.
$$

(Como en la unidad 06, el error absoluto es "aproximado menos exacto":
esto es $\tilde{A}\tilde{\mathbf{x}} = \mathbf{b}$.) Desarrollando y
usando $A\mathbf{x} = \mathbf{b}$:

$$
A\,\Delta\mathbf{x} = -\Delta A\,(\mathbf{x} + \Delta\mathbf{x})
\quad\Longrightarrow\quad
\Delta\mathbf{x} = -A^{-1}\Delta A\,(\mathbf{x} + \Delta\mathbf{x}).
$$

Tomando normas y usando $\|kA\| = |k|\|A\|$ (con $k=-1$) y
$\|AB\| \le \|A\|\|B\|$ (dos veces):

$$
\|\Delta\mathbf{x}\| \le \|A^{-1}\|\,\|\Delta A\|\,\|\mathbf{x} + \Delta\mathbf{x}\|
\le \|A^{-1}\|\,\|\Delta A\|\,\left(\|\mathbf{x}\| + \|\Delta\mathbf{x}\|\right)
\approx \|A^{-1}\|\,\|\Delta A\|\,\|\mathbf{x}\|,
$$

donde en el segundo paso usamos la desigualdad del triángulo y en el
último despreciamos el término de segundo orden
$\|\Delta A\|\,\|\Delta\mathbf{x}\|$ (producto de dos cantidades
chicas). Dividiendo entre $\|\mathbf{x}\|$, y multiplicando y
dividiendo por $\|A\|$:

$$
\frac{\|\Delta\mathbf{x}\|}{\|\mathbf{x}\|} \le \|A\|\,\|A^{-1}\|\,\frac{\|\Delta A\|}{\|A\|}.
$$

Es decir: una cota al cambio *relativo* en $A$ se traduce en una cota
al cambio *relativo* en $\mathbf{x}$, y el factor que los conecta es el
**número de condición**

$$
\kappa(A) = \|A\|\,\|A^{-1}\|.
$$

- Si $\kappa(A)$ es de orden 1, una perturbación chica no se
  amplifica: el problema está **bien condicionado**. (Siempre
  $\kappa(A) \ge 1$.)
- Si $\kappa(A) \gg 1$, una perturbación chica se puede amplificar
  mucho: el problema está **mal condicionado**.

Para el ejemplo de Kahan, $\kappa_\infty(A) \approx 3.3\times 10^{8}$:
un error relativo de $10^{-4}$ en $A$ puede volverse un error relativo
de hasta $10^{4}$ en $\mathbf{x}$. Y el propio redondeo al guardar los
datos ($\sim 10^{-16}$) se puede amplificar hasta $\sim 10^{-8}$. Para
$D = 0.1\,I$, en cambio, $\kappa(D) = 0.1 \times 10 = 1$. **El número
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

$$
A\mathbf{v}_i = \lambda_i\mathbf{v}_i.
$$

No vamos a suponer que $A$ es simétrica (en física hay casos
importantes que no lo son). Entonces conviene distinguir entre los
**eigenvectores derechos** $\mathbf{v}_i$ de arriba y los
**eigenvectores izquierdos** $\mathbf{u}_i$, definidos por

$$
\mathbf{u}_i^T A = \lambda_i\,\mathbf{u}_i^T.
$$

($\mathbf{u}_i^T$ es un vector renglón $1\times n$.) Transponiendo,
$A^T\mathbf{u}_i = \lambda_i\mathbf{u}_i$: **los eigenvectores
izquierdos de $A$ son los eigenvectores derechos de $A^T$**, con los
mismos eigenvalores. Si $A$ no es simétrica, en general
$\mathbf{u}_i \ne \mathbf{v}_i$.

Un dato que usaremos abajo: para eigenvalores distintos,
**$\mathbf{u}_k^T\mathbf{v}_i = 0$ si $k \ne i$**. Demostración:
$\mathbf{u}_k^T A \mathbf{v}_i$ se puede evaluar de dos formas,
$\lambda_i\,\mathbf{u}_k^T\mathbf{v}_i$ (usando $A\mathbf{v}_i$) o
$\lambda_k\,\mathbf{u}_k^T\mathbf{v}_i$ (usando $\mathbf{u}_k^T A$);
restando, $(\lambda_k - \lambda_i)\,\mathbf{u}_k^T\mathbf{v}_i = 0$.

Perturbamos $A$, y con ella cambian eigenvalores y eigenvectores:

$$
(A + \Delta A)(\mathbf{v}_i + \Delta\mathbf{v}_i) = (\lambda_i + \Delta\lambda_i)(\mathbf{v}_i + \Delta\mathbf{v}_i).
$$

Desarrollando, cancelando $A\mathbf{v}_i = \lambda_i\mathbf{v}_i$ y
despreciando los términos de segundo orden ($\Delta\times\Delta$):

$$
A\,\Delta\mathbf{v}_i + \Delta A\,\mathbf{v}_i = \lambda_i\,\Delta\mathbf{v}_i + \Delta\lambda_i\,\mathbf{v}_i.
$$

Multiplicando por la izquierda por $\mathbf{u}_i^T$, el primer término
de cada lado se cancela (porque
$\mathbf{u}_i^T A = \lambda_i\mathbf{u}_i^T$), y queda

$$
\mathbf{u}_i^T\,\Delta A\,\mathbf{v}_i = \Delta\lambda_i\,\mathbf{u}_i^T\mathbf{v}_i
\quad\Longrightarrow\quad
|\Delta\lambda_i| = \frac{|\mathbf{u}_i^T\,\Delta A\,\mathbf{v}_i|}{|\mathbf{u}_i^T\mathbf{v}_i|}.
$$

Por la desigualdad de Cauchy–Schwarz,
$|\mathbf{u}_i^T\,\Delta A\,\mathbf{v}_i| \le \|\mathbf{u}_i\|\,\|\Delta A\|\,\|\mathbf{v}_i\|$,
y si normalizamos los eigenvectores, $\|\mathbf{u}_i\| = \|\mathbf{v}_i\| = 1$
(como hacen las bibliotecas estándar):

$$
|\Delta\lambda_i| \le \frac{1}{|\mathbf{u}_i^T\mathbf{v}_i|}\,\|\Delta A\|.
$$

(Aquí $\|\mathbf{u}\|$ es la norma euclídea; $\|\Delta A\|$ es la norma
de matrices asociada a ella, que siempre es menor o igual que la de
Frobenius, así que la cota vale también con $\|\Delta A\|_F$.)

El factor que amplifica la perturbación es el **número de condición
para el eigenvalor simple** $\lambda_i$:

$$
\kappa^{ev}_{\lambda_i}(A) = \frac{1}{|\mathbf{u}_i^T\mathbf{v}_i|}.
$$

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

$$
A\,\Delta\mathbf{v}_i + \Delta A\,\mathbf{v}_i = \lambda_i\,\Delta\mathbf{v}_i + \Delta\lambda_i\,\mathbf{v}_i,
$$

y desarrollamos la perturbación del eigenvector en términos de los
*demás* eigenvectores:

$$
\Delta\mathbf{v}_i = \sum_{j\ne i} t_{ji}\,\mathbf{v}_j,
$$

con coeficientes $t_{ji}$ por determinar. (No hace falta un término
$j = i$: una componente de $\Delta\mathbf{v}_i$ a lo largo de
$\mathbf{v}_i$ solo cambia su normalización.) Sustituyendo, y usando
$A\mathbf{v}_j = \lambda_j\mathbf{v}_j$:

$$
\sum_{j\ne i}(\lambda_j - \lambda_i)\,t_{ji}\,\mathbf{v}_j + \Delta A\,\mathbf{v}_i = \Delta\lambda_i\,\mathbf{v}_i.
$$

Multiplicamos por la izquierda por $\mathbf{u}_k^T$ con $k \ne i$.
Como $\mathbf{u}_k^T\mathbf{v}_j = 0$ para $j \ne k$, de la suma solo
sobrevive el término $j = k$, y el lado derecho se anula:

$$
(\lambda_k - \lambda_i)\,t_{ki}\,\mathbf{u}_k^T\mathbf{v}_k + \mathbf{u}_k^T\,\Delta A\,\mathbf{v}_i = 0.
$$

Despejando $t_{ki}$ y sustituyendo en el desarrollo de
$\Delta\mathbf{v}_i$:

$$
\Delta\mathbf{v}_i = \sum_{j\ne i}\frac{\mathbf{u}_j^T\,\Delta A\,\mathbf{v}_i}{(\lambda_i - \lambda_j)\,\mathbf{u}_j^T\mathbf{v}_j}\,\mathbf{v}_j.
$$

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
métodos generales (eliminación gaussiana, descomposición LU, que
veremos después) transforman el problema hasta dejarlo en términos de
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

$$
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
$$

La primera ecuación da $x_0$; con él, la segunda da $x_1$; con ambos,
la tercera da $x_2$:

$$
x_0 = \frac{b_0}{L_{00}}, \qquad
x_1 = \frac{b_1 - L_{10}x_0}{L_{11}}, \qquad
x_2 = \frac{b_2 - L_{20}x_0 - L_{21}x_1}{L_{22}}.
$$

Se llama **sustitución hacia adelante** (*forward substitution*)
porque se empieza por la primera ecuación y se avanza. Para
$n\times n$:

$$
x_i = \frac{1}{L_{ii}}\left(b_i - \sum_{j=0}^{i-1} L_{ij}\,x_j\right),
\qquad i = 0, 1, \ldots, n-1,
$$

entendiendo que la suma no tiene términos si $i = 0$, tiene uno si
$i = 1$, etc.

### Sustitución hacia atrás

Si en cambio $U$ es triangular superior (*upper*), $U\mathbf{x} = \mathbf{b}$
es, para $3\times 3$,

$$
\begin{aligned}
U_{00}x_0 + U_{01}x_1 + U_{02}x_2 &= b_0 \\
U_{11}x_1 + U_{12}x_2 &= b_1 \\
U_{22}x_2 &= b_2
\end{aligned}
$$

Ahora se empieza por la **última** ecuación y se avanza hacia atrás:

$$
x_2 = \frac{b_2}{U_{22}}, \qquad
x_1 = \frac{b_1 - U_{12}x_2}{U_{11}}, \qquad
x_0 = \frac{b_0 - U_{01}x_1 - U_{02}x_2}{U_{00}}.
$$

Es la **sustitución hacia atrás** (*back substitution*). Para
$n\times n$:

$$
x_i = \frac{1}{U_{ii}}\left(b_i - \sum_{j=i+1}^{n-1} U_{ij}\,x_j\right),
\qquad i = n-1, n-2, \ldots, 1, 0,
$$

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

$$
2n^2 - n \sim 2n^2 = O(n^2)\ \text{operaciones}.
$$

**Sustitución hacia adelante.** Cada $x_i$ requiere una división ($n$
en total). Además, para $x_i$ hay $i$ multiplicaciones
($L_{ij}x_j$ para $j = 0, \ldots, i-1$) e $i$ sumas/restas ($i-1$
sumas para juntar los $i$ productos, más la resta de $b_i$). Sumando
sobre $i$, con $\sum_{i=0}^{n-1} i = \frac{(n-1)n}{2}$:

$$
\text{sumas/restas: } \frac{n^2 - n}{2},
\qquad
\text{multiplicaciones/divisiones: } n + \frac{n^2 - n}{2} = \frac{n^2 + n}{2}.
$$

En total, **exactamente $n^2$ operaciones**. Se podría escribir
$O(n^2)$, pero eso es menos informativo: el conteo explícito nos dice
que el prefactor es exactamente 1. (La sustitución hacia atrás cuesta
lo mismo.)

`triangulares.py` comprueba ambos conteos agregando un contador a cada
operación, y mide el tiempo de la sustitución hacia adelante al
duplicar $n$: como el costo es $\propto n^2$, el tiempo se cuadruplica.

Para comparar: el producto de dos vectores es $O(n)$ y el de dos
matrices $`O(n^3)`$ (calculen los prefactores exactos, junto con los
términos de grado menor). Y, como veremos, resolver un sistema general
con eliminación gaussiana es $`O(n^3)`$: mucho más caro que resolver uno
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

$$
\begin{aligned}
2x_0 + x_1 + x_2 &= 8 \\
x_0 + x_1 - 2x_2 &= -2 \\
5x_0 + 10x_1 + 5x_2 &= 10
\end{aligned}
\qquad\text{o sea}\qquad
\begin{pmatrix} 2 & 1 & 1 \\ 1 & 1 & -2 \\ 5 & 10 & 5 \end{pmatrix}
\begin{pmatrix} x_0 \\ x_1 \\ x_2 \end{pmatrix}
= \begin{pmatrix} 8 \\ -2 \\ 10 \end{pmatrix}.
$$

Trabajamos con la matriz aumentada, que contiene todo lo que importa
(la $\mathbf{x}$ queda implícita):

$$
(A|\mathbf{b}) = \left(\begin{array}{ccc|c}
2 & 1 & 1 & 8 \\
1 & 1 & -2 & -2 \\
5 & 10 & 5 & 10
\end{array}\right).
$$

Nótese que $A$ no es simétrica. La operación que vamos a repetir es

$$
\text{nuevo renglón } i = \text{renglón } i - \text{coeficiente} \times \text{renglón } j.
$$

El renglón $j$ se llama **renglón pivote**; el renglón $i$ es el que
estamos transformando. El coeficiente se escoge para que el primer
elemento distinto de cero del renglón $i$ se vuelva $0$.

**Pivote $j = 0$.** Para $i = 1$ el coeficiente es $0.5$, porque
$1 - 0.5\times 2 = 0$. La operación se hace con **todo** el renglón,
incluyendo el elemento de $\mathbf{b}$:

$$
\left(\begin{array}{ccc|c}
2 & 1 & 1 & 8 \\
0 & 0.5 & -2.5 & -6 \\
5 & 10 & 5 & 10
\end{array}\right).
$$

Para $i = 2$ el coeficiente es $2.5$, porque $5 - 2.5\times 2 = 0$:

$$
\left(\begin{array}{ccc|c}
2 & 1 & 1 & 8 \\
0 & 0.5 & -2.5 & -6 \\
0 & 7.5 & 2.5 & -10
\end{array}\right).
$$

Ya terminamos con el pivote $j = 0$: la columna $0$ tiene ceros abajo
de la diagonal.

**Pivote $j = 1$.** Siempre se usa la versión **más reciente** de la
matriz, así que el renglón pivote es $(0,\ 0.5,\ -2.5 \,|\, -6)$. Los
renglones que se transforman están siempre abajo del pivote; aquí
solo queda $i = 2$. El coeficiente es $15$, porque
$7.5 - 15\times 0.5 = 0$:

$$
\left(\begin{array}{ccc|c}
2 & 1 & 1 & 8 \\
0 & 0.5 & -2.5 & -6 \\
0 & 0 & 40 & 80
\end{array}\right).
$$

La matriz de coeficientes ya es triangular superior: terminó la
eliminación. Con la sustitución hacia atrás:

$$
x_2 = \frac{80}{40} = 2, \qquad
x_1 = \frac{-6 - (-2.5)\times 2}{0.5} = -2, \qquad
x_0 = \frac{8 - 1\times(-2) - 1\times 2}{2} = 4.
$$

`eliminacion_gaussiana_lu.py` repite este ejemplo e imprime la matriz
aumentada después de cada paso.

### Caso general

Para una matriz $n\times n$, la eliminación modifica $A$ y
$\mathbf{b}$ hasta que $A$ queda triangular. A media eliminación, justo
cuando el renglón $j$ se convierte en pivote, la matriz aumentada se ve
así (con $n = 5$ y $j = 2$; $*$ es un número cualquiera):

$$
\left(\begin{array}{ccccc|c}
* & * & * & * & * & * \\
0 & * & * & * & * & * \\
0 & 0 & A_{22} & * & * & * \\
0 & 0 & * & * & * & * \\
0 & 0 & * & * & * & *
\end{array}\right).
$$

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
  $A_{ij} - (A_{ij}/A_{jj})\,A_{jj} = 0$, el coeficiente es

$$
\text{coeficiente} = \frac{A_{ij}}{A_{jj}}.
$$

$A_{jj}$ se llama **elemento pivote**: es el que se divide para
eliminar los primeros elementos de los renglones de abajo. Con ese
coeficiente, el renglón $i$ se actualiza elemento por elemento:

$$
\begin{aligned}
A_{ik} &\leftarrow A_{ik} - \text{coeficiente}\times A_{jk}, \qquad k = j, j+1, \ldots, n-1, \\
b_i &\leftarrow b_i - \text{coeficiente}\times b_j.
\end{aligned}
$$

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

$$
(n-1-j) + 2(n+1-j)(n-1-j).
$$

Sumando sobre todos los pivotes, con el cambio de variable
$k = n-1-j$, que va de $n-1$ a $1$:

$$
N = \sum_{j=0}^{n-2}\left[(n-1-j) + 2(n+1-j)(n-1-j)\right]
= \sum_{k=1}^{n-1}\left[k + 2(k+2)k\right]
= \sum_{k=1}^{n-1}\left(2k^2 + 5k\right).
$$

Con las sumas conocidas
$\sum_{k=0}^{n-1} k = \frac{(n-1)n}{2}$ y
$\sum_{k=0}^{n-1} k^2 = \frac{(n-1)n(2n-1)}{6}$:

$$
N = 2\,\frac{(n-1)n(2n-1)}{6} + 5\,\frac{(n-1)n}{2}
= \frac{2}{3}n^3 + \frac{3}{2}n^2 - \frac{13}{6}n
\sim \frac{2}{3}n^3.
$$

Para $n = 3$ da $18 + 13.5 - 6.5 = 25$: compruébenlo contando las
operaciones del ejemplo de arriba. La sustitución hacia atrás agrega
$n^2$, pero para $n$ grande la eliminación domina por completo:

$$
\frac{2}{3}n^3 + n^2 \sim \frac{2}{3}n^3.
$$

Es decir: resolver un sistema general cuesta $`O(n^3)`$, mientras que uno
triangular cuesta $n^2$. Para $n = 1000$, la eliminación es unas 670
veces más cara que la sustitución. `eliminacion_gaussiana_lu.py`
cuenta las operaciones con un contador, como `triangulares.py`, y
coinciden exactamente con la fórmula.

### Un pivote cero

La fórmula del coeficiente divide entre $A_{jj}$. Si algún pivote vale
cero, el método truena, **aunque la matriz no sea singular**. Por
ejemplo,

$$
\begin{pmatrix} 0 & 1 \\ 1 & 1 \end{pmatrix}
\begin{pmatrix} x_0 \\ x_1 \end{pmatrix}
= \begin{pmatrix} 1 \\ 2 \end{pmatrix}
$$

tiene determinante $-1$ y solución $\mathbf{x} = (1, 1)$, pero
`eliminacion_gaussiana` levanta `ZeroDivisionError` en el primer
paso. La solución es obvia: intercambiar los dos renglones (la segunda
operación elemental, el **pivoteo**). Lo veremos con cuidado más
adelante; también resuelve un problema menos visible, el de los
pivotes que no son cero pero sí muy chicos.

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

$$
A = LU.
$$

Es la **descomposición LU** (o factorización LU) de $A$. (Más adelante,
con el pivoteo, veremos que la historia es un poco más complicada;
por ahora supongamos que se puede.) La descomposición no es única;
para fijarla pedimos que $L$ tenga **unos en la diagonal**, $L_{ii} = 1$.
Eso se llama descomposición de **Doolittle**.

Veamos cómo se construye en el caso general $3\times 3$:

$$
L = \begin{pmatrix} 1 & 0 & 0 \\ L_{10} & 1 & 0 \\ L_{20} & L_{21} & 1 \end{pmatrix},
\qquad
U = \begin{pmatrix} U_{00} & U_{01} & U_{02} \\ 0 & U_{11} & U_{12} \\ 0 & 0 & U_{22} \end{pmatrix}.
$$

Multiplicándolas:

$$
A = LU = \begin{pmatrix}
U_{00} & U_{01} & U_{02} \\
L_{10}U_{00} & L_{10}U_{01} + U_{11} & L_{10}U_{02} + U_{12} \\
L_{20}U_{00} & L_{20}U_{01} + L_{21}U_{11} & L_{20}U_{02} + L_{21}U_{12} + U_{22}
\end{pmatrix}.
$$

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

$$
\begin{pmatrix} U_{00} & U_{01} & U_{02} \\ 0 & U_{11} & U_{12} \\ 0 & 0 & U_{22} \end{pmatrix} = U.
$$

Leyendo esto al revés, llegamos a la conclusión importante:

- **Para descomponer $A = LU$ basta con hacer la eliminación
  gaussiana.** Lo que queda de $A$ al final es $U$.
- **Los elementos de $L$ abajo de la diagonal son los coeficientes que
  se usaron en la eliminación:** $L_{ij}$ es el coeficiente con el que
  el pivote $j$ eliminó al renglón $i$.

Con el ejemplo $3\times 3$ de la sección anterior, ya sin
$\mathbf{b}$: $U$ es la matriz triangular a la que llegamos, y $L$
junta los coeficientes $0.5$, $2.5$ y $15$:

$$
U = \begin{pmatrix} 2 & 1 & 1 \\ 0 & 0.5 & -2.5 \\ 0 & 0 & 40 \end{pmatrix},
\qquad
L = \begin{pmatrix} 1 & 0 & 0 \\ 0.5 & 1 & 0 \\ 2.5 & 15 & 1 \end{pmatrix}.
$$

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

$$
\begin{aligned}
L\mathbf{y} &= \mathbf{b} \qquad \text{(sustitución hacia adelante)}, \\
U\mathbf{x} &= \mathbf{y} \qquad \text{(sustitución hacia atrás)}.
\end{aligned}
$$

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

### La inversa y el número de condición

La columna $k$ de $A^{-1}$ es la solución de $A\mathbf{x} = \mathbf{e}_k$,
donde $\mathbf{e}_k$ es la columna $k$ de la identidad, porque
$AA^{-1} = I$. Calcular $A^{-1}$ es resolver $n$ sistemas con la
misma $A$: el caso ideal para LU. Una descomposición y $n$ pares de
sustituciones:

$$
\frac{2}{3}n^3 + n\cdot 2n^2 = \frac{8}{3}n^3 \text{ operaciones},
$$

contra $n\cdot\frac{2}{3}n^3$ si hiciéramos una eliminación gaussiana
por columna.

Con $A^{-1}$ por fin podemos calcular el número de condición
$\kappa(A) = \|A\|\,\|A^{-1}\|$ de cualquier matriz, no solo de las
$2\times 2$ de `analisis_de_error.py`. Para la matriz de prueba del
libro, `eliminacion_gaussiana_lu.py` obtiene
$`\kappa_\infty(A) \approx 3.9\times 10^{8}`$: la explicación de los
8 dígitos perdidos de la sección anterior. (Calcular $A^{-1}$ para
resolver un sistema, con $\mathbf{x} = A^{-1}\mathbf{b}$, es mala idea:
cuesta más y acumula más error que LU. Aquí la queremos solo para
$\kappa$.)

### Conteo de operaciones

La descomposición es la eliminación sin $\mathbf{b}$: para el pivote
$j$ se ahorran la multiplicación y la resta de $b_i$ en cada uno de los
$n-1-j$ renglones. Repitiendo la cuenta de la sección anterior con
ese cambio se obtiene

$$
N_{LU} = \frac{2}{3}n^3 + \frac{1}{2}n^2 - \frac{7}{6}n \sim \frac{2}{3}n^3.
$$

(Háganlo: es el mismo procedimiento, con $2k^2 + 3k$ en lugar de
$2k^2 + 5k$. `eliminacion_gaussiana_lu.py` comprueba la fórmula
contando.) El término dominante es el mismo que en la eliminación
gaussiana. La diferencia está en lo que sigue:

| Costo aproximado | Un sistema | $m$ sistemas con la misma $A$ |
|---|---|---|
| Eliminación gaussiana | $\frac{2}{3}n^3 + n^2$ | $m\,\frac{2}{3}n^3$ |
| LU | $\frac{2}{3}n^3 + 2n^2$ | $\frac{2}{3}n^3 + 2mn^2$ |

Para un solo sistema da casi lo mismo (LU hace una sustitución más).
Para muchos, LU gana por mucho: con $n = 100$ y $m = 20$, el conteo
predice que LU es unas 12 veces más rápida, y el script mide alrededor
de 14.

Los flops no son lo único que importa: también la memoria. Tal como lo
implementamos, LU guarda dos matrices $n\times n$, y la eliminación
gaussiana solo una.

## Contenido
- [`fiscomp/algebra_lineal.py`](../../fiscomp/algebra_lineal.py):
  producto matriz-vector (`mat_vec`), residuo, normas de vectores y de
  matrices, partes triangulares de una matriz, las sustituciones
  hacia adelante y hacia atrás, eliminación gaussiana
  (`eliminacion_gaussiana`), descomposición LU (`descomposicion_lu`,
  `resolver_lu`) y la matriz de prueba del libro (`crear_prueba`).
  Todo con listas de listas, sin NumPy. Aquí se irá acumulando el
  resto del tema.
- [`analisis_de_error.py`](analisis_de_error.py): el ejemplo de Kahan
  completo (sección "Análisis de error"): residuo diminuto con una
  solución totalmente equivocada; el cambio drástico de la solución al
  perturbar un solo elemento; por qué el determinante no sirve como
  criterio (con $D = 0.1\,I$); el número de condición $\kappa(A)$ con
  las normas infinito y de Frobenius, y la comprobación de la cota
  $\|\Delta\mathbf{x}\|/\|\mathbf{x}\| \le \kappa(A)\,\|\Delta A\|/\|A\|$.
  Como aún no tenemos un método general, todo es $2\times 2$ y se
  resuelve con la regla de Cramer.
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
  sistema de prueba del libro resuelto con los dos métodos, la inversa
  y $\kappa(A)$ con LU, el conteo de operaciones contra las fórmulas,
  el tiempo de resolver muchos sistemas con la misma $A$, y qué pasa
  con un pivote cero.

Se corren desde la raíz del repositorio (con el `.venv` activado), por
ejemplo `python3 unidades/09_algebra_lineal/triangulares.py`.
