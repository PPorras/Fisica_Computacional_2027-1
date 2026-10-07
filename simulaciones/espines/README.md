# Espines 1/2: la ecuación de Schrödinger como problema de eigenvalores

Proyecto basado en la sección 4.6 del libro de Gezerlis (códigos 4.12
a 4.14). Usa la clase `Matrix` y la descomposición LU de la
[unidad 09](../../unidades/09_algebra_lineal/notas.md), sin NumPy.

El libro calcula los eigenvalores con su función `qrmet()` (el método
QR), que todavía no hemos visto. Aquí hacemos todo lo demás con lo que
ya tenemos: construir los hamiltonianos, y sacar sus eigenvalores con
el determinante calculado con LU. Al final se explica qué cosas sí
necesitan QR.

## La física

La ecuación de Schrödinger independiente del tiempo,
$\hat{H}|\psi\rangle = E|\psi\rangle$, es un problema de eigenvalores.
Para partículas de espín 1/2, sin grados de libertad espaciales (sin
energía cinética), se vuelve directamente un problema de matrices:

```math
H\boldsymbol{\psi} = E\boldsymbol{\psi},
```

con $H$ una matriz y $\boldsymbol{\psi}$ un vector columna. (Que esto
sea así es un resultado, no algo obvio; el apéndice D.2 del libro, en
www.numphyspy.org, lo explica con detalle.)

**Una partícula.** $H$ es $2\times 2$ y proporcional a la matriz de
Pauli $\sigma_z$. Las tres matrices de Pauli son

```math
\sigma_x = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}, \qquad
\sigma_y = \begin{pmatrix} 0 & -i \\ i & 0 \end{pmatrix}, \qquad
\sigma_z = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix},
```

y el operador de espín es $\mathbf{S} = \frac{\hbar}{2}\boldsymbol{\sigma}$.

**Varias partículas: el producto de Kronecker.** Para dos o más
partículas, las matrices de una partícula se combinan con el
**producto de Kronecker** $\otimes$. Si $U$ es $n\times n$ y $V$ es
$p\times p$, $W = U\otimes V$ es $np\times np$, con

```math
W_{ab} = U_{ik}\,V_{jl}, \qquad a = pi + j, \quad b = pk + l,
```

donde $i, k$ van de $0$ a $n-1$ y $j, l$ de $0$ a $p-1$. En palabras:
$W$ está formada por bloques de $p\times p$, y el bloque $(i, k)$ es
$U_{ik}V$. El producto **no conmuta**: $U\otimes V \neq V\otimes U$.

Con él se construyen los operadores de cada partícula, poniendo la
matriz de esa partícula en su lugar e identidades $I$ de $2\times 2$ en
los demás. Para dos partículas, que el libro llama I y II:

```math
S_{Iz} = S_z\otimes I, \qquad S_{IIx} = I\otimes S_x,
```

y para tres, por ejemplo, $S_{IIy} = I\otimes S_y\otimes I$.

**Dos espines**, con una matriz de $4\times 4$:

```math
H = -\omega_I S_{Iz} - \omega_{II} S_{IIz} +
\gamma\left(S_{Ix}S_{IIx} + S_{Iy}S_{IIy} + S_{Iz}S_{IIz}\right).
```

Cada espín interactúa con un campo magnético externo (los términos con
$\omega$, proporcionales al campo) y los dos interactúan entre sí (el
término con $\gamma$; se abrevia como $`\gamma\,\mathbf{S}_I\cdot\mathbf{S}_{II}`$).

**Tres espines**, con una matriz de $8\times 8$:

```math
H = -\omega_I S_{Iz} - \omega_{II} S_{IIz} - \omega_{III} S_{IIIz} +
\gamma\left(\mathbf{S}_I\cdot\mathbf{S}_{II} + \mathbf{S}_I\cdot\mathbf{S}_{III} +
\mathbf{S}_{II}\cdot\mathbf{S}_{III}\right).
```

En todo el proyecto usamos $\hbar = 1$, y en las gráficas también
$\gamma = 1$: las energías quedan en unidades de $\gamma$.

## La programación

### Matrices complejas con `Matrix`

$\sigma_y$ tiene números complejos. Python los trae integrados (se
escriben `1j`), y `Matrix` los acepta en `data` sin cambiar nada:
suma, resta y producto de matrices funcionan igual. Al final, los
productos $S_{Iy}S_{IIy}$ tienen $i\cdot i = -1$ y todo $H$ queda real;
`parte_real` se queda con la parte real (en el libro, `H.real`).

### Operadores de espín

El libro escribe los operadores de cada partícula uno por uno, por
ejemplo `kron(kron(iden, pa), iden)` para la partícula II de tres.
`operadores_de_espin(posicion, numero_de_particulas)` hace lo mismo
para cualquier número de partículas: arma la lista de factores (la
matriz de Pauli en `posicion`, identidades en los demás) y los une con
`kron`. Así, el caso de cuatro espines (problema 4.50) sale sin
escribir nada nuevo.

### Un tropiezo: `sum()` con matrices

El libro suma los tres productos de $`\mathbf{S}_I\cdot\mathbf{S}_{II}`$
con `sum([...])`. Con `Matrix` eso falla:

```
TypeError: unsupported operand type(s) for +: 'int' and 'Matrix'
```

porque `sum()` empieza en `0` y hace `0 + matriz`, y `Matrix` no sabe
sumarse a un número (le falta `__radd__`). `producto_punto_espin`
suma los tres términos a mano, como sugiere el propio libro.

## Los eigenvalores, sin QR

### Dos espines: el bloque $2\times 2$

En la base $(\uparrow\uparrow, \uparrow\downarrow, \downarrow\uparrow, \downarrow\downarrow)$,
el hamiltoniano de dos espines es

```math
H = \begin{pmatrix}
-\frac{\omega_I + \omega_{II}}{2} + \frac{\gamma}{4} & 0 & 0 & 0 \\
0 & -\frac{\omega_I - \omega_{II}}{2} - \frac{\gamma}{4} & \frac{\gamma}{2} & 0 \\
0 & \frac{\gamma}{2} & \frac{\omega_I - \omega_{II}}{2} - \frac{\gamma}{4} & 0 \\
0 & 0 & 0 & \frac{\omega_I + \omega_{II}}{2} + \frac{\gamma}{4}
\end{pmatrix}.
```

(Compruébenlo imprimiendo `dos_espines(1.0, 2.0, 0.5)`.) Los estados
con los dos espines alineados no se mezclan con nada: sus energías son
los elementos de la diagonal. Los otros dos forman un bloque
$2\times 2$, cuyos eigenvalores salen del polinomio característico,
como en `condicion_eigenvalores.py` de la unidad 09:

```math
E_\pm = -\frac{\gamma}{4} \pm \sqrt{\left(\frac{\omega_I - \omega_{II}}{2}\right)^2 + \left(\frac{\gamma}{2}\right)^2}.
```

Con esto se obtiene la gráfica de dos espines completa
(`graficar_dos_espines.gp`, panel izquierdo de la Fig. 4.1):

- En $\omega = 0$ (sin campo), tres niveles coinciden en $\gamma/4$: es
  el **triplete** de espín total 1. El cuarto, $-3\gamma/4$, es el
  singulete.
- Con $\omega_{II} = 2\omega_I$, el nivel $\uparrow\uparrow$ cruza al
  nivel menor del bloque en $\omega_I = 3\gamma/4$. Los dos niveles se
  pueden cruzar porque no se mezclan: uno está fuera del bloque.
  (Es una situación parecida al efecto Zeeman en el hidrógeno.)

### Tres espines: el determinante con LU

Con tres espines ya no hay atajo: hay bloques de $3\times 3$, que
llevarían a resolver cúbicas. Pero los eigenvalores son las raíces de

```math
p(\lambda) = \det(H - \lambda I),
```

y el determinante se calcula fácil con LU: como $L$ tiene unos en la
diagonal, $\det(A) = \det(L)\det(U)$ es el producto de la diagonal de
$U$. Eso cuesta $`\sim 2n^3/3`$ operaciones, contra $n!$ de la fórmula
con cofactores.

El método (en [`eigenvalores.py`](eigenvalores.py)):

1. **Dónde buscar.** $H$ es real y simétrica, así que sus eigenvalores
   son reales, y todos cumplen $`|\lambda| \le \|H\|_\infty`$ (la norma
   de la unidad 09). Basta barrer ese intervalo.
2. **Barrido.** Se evalúa $p(\lambda)$ en una malla de puntos y se
   buscan los cambios de signo: entre dos puntos donde $p$ cambia de
   signo hay un eigenvalor. Con 2000 puntos, unos 3 decimales.
3. **Refinamiento.** Cada intervalo con cambio de signo se divide en
   10, nos quedamos con el pedazo donde sigue cambiando el signo, y se
   repite: un dígito más por ronda.
4. **Iteración inversa** (un adelanto del tema de eigenvalores). Con
   un valor aproximado $\sigma$, se factoriza $H - \sigma I = LU$ una
   sola vez y se repite
   $`(H - \sigma I)\,\mathbf{x}_{\text{nuevo}} = \mathbf{x}`$, resolviendo
   con sustitución hacia adelante y hacia atrás (`resolver_lu`), y
   normalizando. $\mathbf{x}$ converge al eigenvector cuyo eigenvalor
   está más cerca de $\sigma$, y el eigenvalor sale con precisión de
   máquina del cociente de Rayleigh, $\mathbf{x}\cdot H\mathbf{x} / \mathbf{x}\cdot\mathbf{x}$.
   Es el ejemplo perfecto de para qué sirve LU: una factorización
   cara, $`O(n^3)`$, y muchas sustituciones baratas, $`O(n^2)`$.

`niveles_de_energia.py` comprueba el resultado sin conocer la
respuesta: el residuo $`\|H\mathbf{v} - \lambda\mathbf{v}\|`$ es del
orden de $10^{-16}$; la suma de los eigenvalores es la traza de $H$
(cero); y los estados con los tres espines alineados, que no se
mezclan con nada, dan exactamente
$\mp(\omega_I + \omega_{II} + \omega_{III})/2 + 3\gamma/4$.

### Lo que todavía no se puede (hace falta QR)

- **Niveles degenerados.** Sin campo, con $\omega = 0$, los ocho niveles
  de tres espines forman dos grupos de cuatro: $p(\lambda)$ tiene dos
  raíces de multiplicidad 4. Una raíz de multiplicidad **par** no
  cambia el signo de $p$, así que el barrido no la ve: en $\omega = 0$
  no encuentra ningún eigenvalor. (Lo mismo pasa con el $0$ doble del
  caso sin interacción, $\gamma = 0$.) Contar multiplicidades es el
  problema 4.48.
- **Cruces de niveles.** Si dos eigenvalores caen en la misma celda de
  la malla, sus dos cambios de signo se cancelan y se pierden ambos.
  En la gráfica de tres espines eso deja huecos (por ejemplo, cerca de
  $`\omega_I = 0.45`$).
- **Pivoteo.** `descomposicion_lu` no pivotea: si en algún $\lambda$
  aparece un pivote cero, ese punto de la malla se salta.
- **Muchas partículas.** Con $N$ espines, $H$ es de $2^N\times 2^N$ y
  cada determinante cuesta $`\sim \frac{2}{3}8^N`$ operaciones. Para
  $N = 10$ (problema 4.51), con $1024\times 1024$, una sola LU en
  Python puro tarda minutos, y el barrido necesita miles.

## Archivos

- [`espines.py`](espines.py): matrices de Pauli, `kron`,
  `operadores_de_espin`, `producto_punto_espin`, y los hamiltonianos
  `dos_espines` y `tres_espines`. Corriéndolo solo, hace la prueba del
  libro: $\sigma_x$ por Kronecker con una matriz de unos de
  $3\times 3$, y al revés.
- [`eigenvalores.py`](eigenvalores.py): el polinomio característico
  (con `determinante` de `fiscomp/algebra_lineal.py`), barrido,
  refinamiento, iteración inversa, y el atajo del bloque $2\times 2$
  para dos espines.
- [`niveles_de_energia.py`](niveles_de_energia.py): el programa
  principal; imprime los resultados y las comprobaciones, y guarda los
  datos de las gráficas en `datos/`.
- [`graficar_dos_espines.gp`](graficar_dos_espines.gp) y
  [`graficar_tres_espines.gp`](graficar_tres_espines.gp): los dos
  paneles de la Fig. 4.1 del libro.

## Cómo correrlo

Con el entorno virtual activado (ver
[`recursos/notas_entorno_virtual.md`](../../recursos/notas_entorno_virtual.md)),
desde esta carpeta:

```bash
python3 espines.py
python3 niveles_de_energia.py
gnuplot graficar_dos_espines.gp
gnuplot graficar_tres_espines.gp
```

## Para explorar

1. **Las matrices.** Impriman $S_{Iz}$, $S_{IIx}$ y $S_{IIIx}$ (con
   `operadores_de_espin` y `parte_real`) y compárenlas con lo que
   calculen a mano. ¿Qué cambia en $\sigma_x\otimes(\text{unos})$ si
   ponen la matriz de Pauli en el segundo lugar? ¿Por qué salen bloques
   de $3\times 3$ en un caso y de $2\times 2$ en el otro?
2. **El cruce.** Deduzcan a mano que, con $\omega_{II} = 2\omega_I$, el
   nivel $\uparrow\uparrow$ cruza al nivel menor del bloque en
   $\omega_I = 3\gamma/4$, y compárenlo con la gráfica.
3. **Los huecos.** Aumenten `puntos` en el barrido de tres espines.
   ¿Desaparecen los huecos cerca de $\omega_I = 0.45$? ¿Y el de
   $\omega_I = 0$? ¿Por qué uno sí y el otro no?
4. **Iteración inversa cerca de un cruce.** Cuenten cuántas
   iteraciones necesita `iteracion_inversa` para que el eigenvalor deje
   de cambiar, lejos de un cruce y cerca de uno. ¿Por qué converge más
   lento cuando hay otro eigenvalor cerca de $\sigma$?
5. **Un término nuevo.** Agreguen a `tres_espines` un término
   $-\beta S_{Ix}$ (un campo en la dirección $x$ sobre la partícula I).
   Con eso, el atajo de los bloques ya no sirve ni siquiera para dos
   espines, pero el barrido funciona igual. (Es justo el punto del
   libro: con el polinomio característico a mano, cualquier término
   nuevo complica todo; con la computadora, no.)
6. **Cuatro espines** (problema 4.50). Escriban `cuatro_espines` con
   `operadores_de_espin(posicion, 4)`, que da matrices de
   $16\times 16$, y repitan la
   gráfica, con $\omega_{IV} = 4\omega_I$. Midan cuánto tarda contra el
   caso de tres espines y compárenlo con lo que predice el conteo de
   operaciones.
