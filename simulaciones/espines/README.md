# Espines 1/2: la ecuación de Schrödinger como problema de eigenvalores

Este proyecto se basa en la sección 4.6 de [1] (códigos 4.12 a 4.14 y
problemas 4.48 a 4.51). Construye los hamiltonianos con la clase
`Matrix` y calcula sus eigenvalores con el método QR de la
[unidad 10](../../unidades/10_valores_y_vectores_propios/), sin NumPy.

## La física

La ecuación de Schrödinger independiente del tiempo,
$\hat{H}|\psi\rangle = E|\psi\rangle$, es un problema de eigenvalores.
Para partículas de espín 1/2, sin grados de libertad espaciales (sin
energía cinética), se vuelve directamente un problema de matrices:

```math
H\boldsymbol{\psi} = E\boldsymbol{\psi},
```

con $H$ una matriz y $\boldsymbol{\psi}$ un vector columna. (Que esto
sea así es un resultado, no algo obvio; el apéndice D.2 de [1], en
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
los demás. A las partículas les llamamos I, II, III (con números
romanos, para no confundirlas con los índices de las componentes
cartesianas ni con los de las matrices):

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

**$N$ espines**, con una matriz de $`2^N\times 2^N`$: un término de campo
por partícula y un término de interacción por cada **par**,

```math
H = -\sum_{k} \omega_k S_{kz} + \gamma\sum_{a < b} \mathbf{S}_a\cdot\mathbf{S}_b .
```

En todo el proyecto usamos $\hbar = 1$, y en las gráficas también
$\gamma = 1$: las energías quedan en unidades de $\gamma$.

## La programación

### Matrices complejas con `Matrix`

$\sigma_y$ tiene números complejos. Python los trae integrados (se
escriben `1j`), y `Matrix` los acepta en `data` sin cambiar nada:
suma, resta y producto de matrices funcionan igual. Al final, los
productos $S_{Iy}S_{IIy}$ tienen $i\cdot i = -1$ y todo $H$ queda real;
`parte_real` se queda con la parte real (con NumPy sería `H.real`).
Nos conviene: el método QR que programamos solo trabaja con números
reales.

### Operadores de espín

Para tres partículas, los operadores de la partícula II se pueden
escribir uno por uno, `kron(kron(iden, pa), iden)`. En cambio,
`operadores_de_espin(posicion, numero_de_particulas)` sirve para
cualquier número de partículas: arma la lista de factores (la matriz
de Pauli en `posicion`, identidades en los demás) y los une con
`kron`. Con ella, `hamiltoniano(omegas, gamma)` arma el $H$ de $N$
espines de arriba, para cualquier $N$; `dos_espines` y `tres_espines`
lo escriben explícitamente, término por término, para que se vea la
correspondencia con las ecuaciones.

### Un tropiezo: `sum()` con matrices

Con NumPy, los tres productos de $`\mathbf{S}_I\cdot\mathbf{S}_{II}`$
se podrían sumar con `sum([...])`. Con `Matrix` eso falla:

```
TypeError: unsupported operand type(s) for +: 'int' and 'Matrix'
```

porque `sum()` empieza en `0` y hace `0 + matriz`, y `Matrix` no sabe
sumarse a un número (le falta `__radd__`). `producto_punto_espin`
suma los tres términos a mano.

## Los eigenvalores, con el método QR

### El desplazamiento

El método QR (`metodo_qr`, unidad 10) da todos los eigenvalores, pero
tal como lo programamos tiene dos limitaciones, y los hamiltonianos de
espín se topan con las dos (problema 4.48 de [1]):

- **Sin interacción** ($\gamma = 0$), $H$ es diagonal. Con
  $`\omega_I + \omega_{II} = \omega_{III}`$, como en $(1, 2, 3)$, los
  estados $\uparrow\uparrow\downarrow$ y $\downarrow\downarrow\uparrow$
  tienen energía $0$: $H$ es **singular**, y Gram-Schmidt divide entre
  $`R_{jj} = 0`$ (`ZeroDivisionError`).
- **Sin campo** ($\omega = 0$), los eigenvalores de tres espines son
  $`\pm 3\gamma/4`$: el **mismo valor absoluto** con signo contrario. El
  método QR no los separa, y la diagonal de $`A^{(k)}`$ da números que
  no son eigenvalores.

El remedio es el mismo de la SVD: aplicar QR a $H - sI$ y sumar $s$
al final. Aquí escogemos $`s = -(\|H\|_\infty + 0.1)`$. Como $H$ es
simétrica, sus eigenvalores cumplen $`|\lambda_i| \le \|H\|_\infty`$,
así que los de $H - sI$ quedan entre $0.1$ y $`2\|H\|_\infty + 0.1`$:
**todos positivos**. Entonces dos eigenvalores distintos tienen valores
absolutos distintos, y ninguno es cero. Los eigenvalores repetidos de
verdad (los degenerados) no le estorban a QR. Eso hace
`eigenvalores_qr(H)` en [`eigenvalores.py`](eigenvalores.py).

Hay un precio: los cocientes $`|\lambda_{i+1} - s|/|\lambda_i - s|`$
quedan más cerca de 1, y QR converge más lento. Por eso
`eigenvalores_qr` hace 500 iteraciones, en vez de las 100 de
`metodo_qr`. (Cerca de $\omega = 0$, con cuatro espines, hacen falta
casi todas: muchos niveles están muy juntos.)

### Dos espines

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

(Compruébenlo con la parte 1 de `niveles_de_energia.py`.) Los estados
con los dos espines alineados no se mezclan con nada: sus energías son
los elementos de la diagonal. Los otros dos forman un bloque
$2\times 2$, cuyos eigenvalores salen del polinomio característico:

```math
E_\pm = -\frac{\gamma}{4} \pm \sqrt{\left(\frac{\omega_I - \omega_{II}}{2}\right)^2 + \left(\frac{\gamma}{2}\right)^2}.
```

Para dos espines tenemos entonces la respuesta exacta, y sirve para
comprobar QR: coinciden en todos los dígitos. En la gráfica
(`graficar_dos_espines.gp`, panel izquierdo de la Fig. 4.1 de [1]),
los puntos son de QR y las líneas son las exactas:

- En $\omega = 0$ (sin campo), tres niveles coinciden en $\gamma/4$: es
  el **triplete** de espín total 1. El cuarto, $-3\gamma/4$, es el
  singulete (espín total 0).
- Con $\omega_{II} = 2\omega_I$, el nivel $\uparrow\uparrow$ cruza al
  nivel menor del bloque en $\omega_I = 3\gamma/4$. Los dos niveles se
  pueden cruzar porque no se mezclan: uno está fuera del bloque. (Es
  una situación parecida al efecto Zeeman en el hidrógeno.)

### Tres espines

Con tres espines ya no hay atajo fácil: hay bloques de $3\times 3$,
que llevarían a resolver cúbicas. QR da los ocho eigenvalores sin más
trabajo, y `niveles_de_energia.py` los comprueba sin conocer la
respuesta: la suma de los eigenvalores es la traza de $H$ (cero), y
los estados con los tres espines alineados, que no se mezclan con
nada, dan exactamente $\mp(\omega_I + \omega_{II} + \omega_{III})/2 + 3\gamma/4$.

En la gráfica (`graficar_tres_espines.gp`, panel derecho de la
Fig. 4.1 de [1]), con $\omega_{II} = 2\omega_I$ y
$\omega_{III} = 3\omega_I$, hay varios cruces de niveles. En
$\omega = 0$ parece que hay dos grupos de cuatro niveles, en
$\pm 3\gamma/4$. La historia es un poco más complicada (problema 4.48 de
[1]). El **espín total** $`\mathbf{S} = \mathbf{S}_I + \mathbf{S}_{II} + \mathbf{S}_{III}`$
cumple

```math
\mathbf{S}^2 = \mathbf{S}_I^2 + \mathbf{S}_{II}^2 + \mathbf{S}_{III}^2
+ 2\left(\mathbf{S}_I\cdot\mathbf{S}_{II} + \mathbf{S}_I\cdot\mathbf{S}_{III} + \mathbf{S}_{II}\cdot\mathbf{S}_{III}\right),
```

y cada $`\mathbf{S}_k^2 = \frac{1}{2}\left(\frac{1}{2} + 1\right) = \frac{3}{4}`$.
Así que sin campo

```math
H = \frac{\gamma}{2}\left(\mathbf{S}^2 - \frac{9}{4}\right),
```

y como los eigenvalores de $`\mathbf{S}^2`$ son $S(S+1)$, la energía
solo depende del espín total $S$: $E = 3\gamma/4$ para $S = 3/2$ y
$E = -3\gamma/4$ para $S = 1/2$. `espin_total_cuadrado(3)` arma
$`\mathbf{S}^2`$ y QR da sus eigenvalores: $15/4$ cuatro veces y $3/4$
cuatro veces. El grupo de arriba es un **cuadruplete** ($S = 3/2$,
con $M = -3/2, \ldots, 3/2$), pero el de abajo **no**: son **dos
dobletes** ($S = 1/2$, con $M = \pm 1/2$ cada uno). Con dos espines,
lo mismo da el triplete ($`\mathbf{S}^2 = 2`$) y el singulete
($`\mathbf{S}^2 = 0`$).

### Cuatro espines

Con `hamiltoniano` no hay que escribir nada nuevo para cuatro espines
(problema 4.50 de [1]): $H$ es de $16\times 16$, y
`graficar_cuatro_espines.gp` grafica los 16 niveles con
$`\omega_k = k\,\omega_I`$. Sin campo, $H = \frac{\gamma}{2}(\mathbf{S}^2 - 3)$,
y los niveles se agrupan según el espín total: $S = 2$ (5 estados) en
$3\gamma/2$, $S = 1$ (tres tripletes, 9 estados) en $-\gamma/2$ y
$S = 0$ (dos singuletes) en $-3\gamma/2$.

### ¿Y los eigenvectores?

La función `eigenpares` de la unidad 10 (QR más potencia inversa)
también daría los eigenvectores, pero con estos hamiltonianos falla
(problema 4.49 de [1]):

- $H$ no mezcla estados con distinto número de espines hacia arriba:
  con tres espines, en la base
  $(\uparrow\uparrow\uparrow, \uparrow\uparrow\downarrow, \ldots, \downarrow\downarrow\downarrow)$,
  queda formada por bloques de $1$, $3$, $3$ y $1$ estados (se ve en la
  parte 2 de `niveles_de_energia.py`). Por eso cada eigenvector tiene
  **ceros exactos** fuera de su bloque, y el criterio de paro de
  `potencia_inversa`, `cambio_relativo`, divide entre cada componente:
  `ZeroDivisionError`.
- Sin campo, con eigenvalores repetidos, la potencia inversa con el
  mismo desplazamiento y el mismo vector inicial daría el mismo vector
  varias veces, no una base del subespacio.

Para la física de este proyecto (los niveles de energía) bastan los
eigenvalores.

### Antes de QR: el polinomio característico

La primera versión de este proyecto se hizo antes de ver QR, con lo
que había en la unidad 09: las raíces de

```math
p(\lambda) = \det(H - \lambda I),
```

con el determinante calculado con LU, buscando cambios de signo de $p$
en una malla de valores de $\lambda$ y refinando cada intervalo. Las
funciones siguen en `eigenvalores.py` (`eigenvalores_por_barrido`)
para comparar, y en los casos normales coinciden con QR en unos 8
dígitos. Pero tienen problemas que QR no tiene:

- Una raíz de multiplicidad **par** no cambia el signo de $p$, así que
  el barrido no la ve: sin campo, con tres espines, no encuentra
  **ninguno** de los ocho eigenvalores.
- Si dos eigenvalores caen en la misma celda de la malla (cerca de un
  cruce), sus dos cambios de signo se cancelan y se pierden ambos.
- `descomposicion_lu` no pivotea: si en algún $\lambda$ aparece un
  pivote cero, ese punto de la malla se salta.

Y en general, encontrar raíces de polinomios es un problema mal
condicionado: por eso no se calculan así los eigenvalores.

### Diez espines

El problema 4.51 de [1] pide repetir todo para diez espines. Con
`hamiltoniano` el código es el mismo, pero $H$ es de
$`2^{10}\times 2^{10} = 1024\times 1024`$. Cada iteración de QR cuesta
$`\sim 4n^3`$ operaciones, así que pasar de $16\times 16$ a
$1024\times 1024$ multiplica el tiempo por $`64^3 \approx 2.6\times 10^5`$:
la parte 6 de `niveles_de_energia.py` lo estima en **un día por cada
valor de $\omega$**, en Python puro (sin contar construir $H$, que
también tiene productos de matrices de $1024\times 1024$). Con NumPy,
que hace las operaciones en C, son segundos. Y quien estudia sistemas
de muchos espines usa además que $H$ es casi toda ceros (se dice que es
**rala**) y sus bloques.

## Archivos

- [`espines.py`](espines.py): matrices de Pauli, `kron`,
  `operadores_de_espin`, `producto_punto_espin`, los hamiltonianos
  `dos_espines`, `tres_espines` y `hamiltoniano` ($N$ espines), y
  `espin_total_cuadrado`. Corriéndolo solo, hace la prueba del código
  4.12 de [1]: $\sigma_x$ por Kronecker con una matriz de unos de
  $3\times 3$, y al revés.
- [`eigenvalores.py`](eigenvalores.py): `eigenvalores_qr` (el método QR
  desplazado), `multiplicidades` para comprobar, el
  barrido del polinomio característico (la versión anterior a QR) y
  los eigenvalores exactos de dos espines (bloque $2\times 2$).
- [`niveles_de_energia.py`](niveles_de_energia.py): el programa
  principal; imprime los resultados y las comprobaciones (incluyendo
  los casos límite, el espín total y por qué falla `eigenpares`),
  guarda los datos de las gráficas en `datos/` y estima el caso de
  diez espines. Tarda unos 40 segundos, casi todo en los cuatro
  espines.
- [`graficar_dos_espines.gp`](graficar_dos_espines.gp) y
  [`graficar_tres_espines.gp`](graficar_tres_espines.gp): los dos
  paneles de la Fig. 4.1 de [1].
- [`graficar_cuatro_espines.gp`](graficar_cuatro_espines.gp): lo mismo
  para cuatro espines.

## Cómo correrlo

Con el entorno virtual activado (ver
[`recursos/notas_entorno_virtual.md`](../../recursos/notas_entorno_virtual.md)),
desde esta carpeta:

```bash
python3 espines.py
python3 niveles_de_energia.py
gnuplot graficar_dos_espines.gp
gnuplot graficar_tres_espines.gp
gnuplot graficar_cuatro_espines.gp
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
3. **El desplazamiento.** Cambien en `eigenvalores_qr` el desplazamiento
   por $s = 0$ (es decir, `metodo_qr` tal cual) y por $s = 1$, y vuelvan
   a hacer la gráfica de tres espines. ¿En qué valores de $\omega_I$
   falla cada uno? Para cada falla, busquen la causa: dos eigenvalores
   distintos con $`|\lambda_i - s|`$ y $`|\lambda_j - s|`$ iguales (o muy
   parecidos), o un $`\lambda_i = s`$.
4. **Cuántas iteraciones.** Con cuatro espines, comparen
   `eigenvalores_qr(H, kmax)` con `kmax` = 100, 200 y 500 en
   $\omega_I = 0.1$ y en $\omega_I = 1$. ¿Por qué cerca de $\omega = 0$
   hacen falta más?
5. **Un término nuevo.** Agreguen a `tres_espines` un término
   $-\beta S_{Ix}$ (un campo en la dirección $x$ sobre la partícula I).
   Ahora $H$ sí mezcla estados con distinto número de espines hacia
   arriba, el atajo de los bloques ya no sirve, pero QR funciona
   igual. (Con el polinomio característico a mano, cualquier término
   nuevo complica todo; con la computadora, no.) ¿Funciona ahora
   `eigenpares`?
6. **El espín total con cuatro espines.** Calculen los eigenvalores de
   `espin_total_cuadrado(4)` con sus multiplicidades y comprueben que
   hay 5 estados con $S = 2$, 9 con $S = 1$ y 2 con $S = 0$.
7. **¿Hasta dónde llegamos?** Midan cuánto tarda `eigenvalores_qr` con
   5 y 6 espines ($32\times 32$ y $64\times 64$). ¿Crece como $n^3$?

## Referencias

[1] A. Gezerlis, *Numerical Methods in Physics with Python*, vol. 1.
Cambridge University Press, Cambridge, UK, 2023.
