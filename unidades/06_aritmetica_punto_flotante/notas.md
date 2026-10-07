# Unidad 06 — Aritmética de punto flotante

## Temas
- **Física computacional:** representación de números en punto
  flotante (IEEE 754), error de redondeo, por qué `0.1 + 0.2 != 0.3`,
  comparación de flotantes con una tolerancia (`abs(a - b) < eps`).

## Definiciones

Un sistema de punto flotante queda determinado por:

- $\beta$ — base o radix (*base or radix*)
- $t$ — precisión (*precision*)
- $[L, U]$ — rango del exponente (*exponent range*)

Por definición, cualquier número $x$ en el sistema de punto flotante se
representa como

```math
x = \pm\left(d_0 + \frac{d_1}{\beta} + \frac{d_2}{\beta^2} + \cdots + \frac{d_{t-1}}{\beta^{t-1}}\right)\beta^e,
```

donde

```math
0 \le d_i \le \beta - 1, \qquad i = 0, \dots, t-1,
```

y

```math
L \le e \le U.
```

### Normalización

Un sistema de punto flotante está **normalizado** si el dígito líder
$d_0$ siempre es distinto de cero (a menos que el número representado
sea el cero). Así, en un sistema normalizado, la mantisa $m$ de un
número de punto flotante distinto de cero siempre satisface

```math
1 \le m < \beta.
```

(Una convención alternativa es que $d_0$ siempre sea cero, en cuyo
caso un número de punto flotante está normalizado si $d_1 \neq 0$, y
en su lugar $\beta^{-1} \le m < 1$.)

Los sistemas de punto flotante suelen normalizarse porque:

- La representación de cada número es entonces única.
- No se desperdician dígitos en ceros a la izquierda, con lo que se
  maximiza la precisión.
- En un sistema binario ($\beta = 2$), el bit líder siempre es 1 y por
  lo tanto no necesita almacenarse, ganando así un bit extra de
  precisión para un ancho de campo dado.

### Propiedades

Un sistema de números de punto flotante es **finito y discreto**.

**Número de elementos.** El número de números de punto flotante
normalizados es

```math
2(\beta - 1)\,\beta^{t-1}\,(U - L + 1) + 1,
```

porque hay dos posibles signos, $\beta - 1$ opciones para el dígito
líder de la mantisa, $\beta$ opciones para cada uno de los $t - 1$
dígitos restantes de la mantisa, y $U - L + 1$ valores posibles para
el exponente. El $+1$ final se agrega porque el número también puede
ser cero.

Existe un menor número positivo normalizado,

```math
\text{Underflow level} = \text{UFL} = \beta^{L},
```

que tiene un 1 como dígito líder y 0 en el resto de los dígitos de la
mantisa, junto con el menor valor posible del exponente. Existe un
mayor número de punto flotante,

```math
\text{Overflow level} = \text{OFL} = \beta^{U+1}(1 - \beta^{-t}),
```

que tiene $\beta - 1$ como valor de cada dígito de la mantisa y el
mayor valor posible del exponente. Cualquier número mayor que OFL no
se puede representar en el sistema de punto flotante dado, ni tampoco
ningún número positivo menor que UFL.

**Cálculo de OFL.** El mayor número de punto flotante se obtiene con
la mayor mantisa posible ($d_i = \beta - 1$ para todo $i$) y el mayor
exponente posible ($e = U$). La mantisa máxima es

```math
m_{\max} = (\beta - 1) + \frac{\beta - 1}{\beta} + \frac{\beta - 1}{\beta^2} + \cdots + \frac{\beta - 1}{\beta^{t-1}}
= (\beta - 1) \sum_{i=0}^{t-1} \beta^{-i}.
```

La suma geométrica vale

```math
\sum_{i=0}^{t-1} \beta^{-i} = \frac{1 - \beta^{-t}}{1 - \beta^{-1}} = \frac{\beta}{\beta - 1}\left(1 - \beta^{-t}\right),
```

así que

```math
m_{\max} = (\beta - 1) \cdot \frac{\beta}{\beta - 1}\left(1 - \beta^{-t}\right) = \beta\left(1 - \beta^{-t}\right).
```

Multiplicando por $\beta^{U}$ (el mayor factor de escala) se obtiene

```math
\text{OFL} = m_{\max} \cdot \beta^{U} = \beta^{U+1}\left(1 - \beta^{-t}\right).
```

*Ejemplo (doble precisión IEEE 754, $`\beta = 2`$, $`t = 53`$, $`U = 1023`$):*

```math
\text{OFL} = 2^{1024}\left(1 - 2^{-53}\right) \approx 1.7977 \times 10^{308},
```

que coincide con el valor máximo representable de un `float` de Python
(`sys.float_info.max`).

Los números de punto flotante **no están distribuidos uniformemente**
en su rango: están igualmente espaciados solo entre potencias
sucesivas de $\beta$. No todos los números reales son representables
exactamente en un sistema de punto flotante; a los números reales que
sí lo son en un sistema dado se les llama a veces **números de
máquina** (*machine numbers*).

### Redondeo (Rounding)

Si un número real $x$ dado no es representable exactamente como
número de punto flotante, entonces debe aproximarse por algún número
de punto flotante "cercano". A la aproximación de punto flotante de un
número real $x$ la denotamos $fl(x)$. El proceso de elegir un número
de punto flotante cercano $fl(x)$ para aproximar un número real $x$ se
llama **redondeo** (*rounding*), y el error introducido por dicha
aproximación se llama **error de redondeo** (*rounding error* o
*roundoff error*).

Dos de las reglas de redondeo más usadas son:

- **Truncamiento (Chop):** la expansión en base $\beta$ de $x$ se
  trunca después del dígito $(t-1)$-ésimo. Como $fl(x)$ es entonces el
  siguiente número de punto flotante hacia el cero desde $x$, a esta
  regla también se le llama a veces **redondeo hacia cero** (*round
  toward zero*).
- **Redondeo al más cercano (Round to nearest):** $fl(x)$ es el número
  de punto flotante más cercano a $x$; en caso de empate, se usa el
  número de punto flotante cuyo último dígito almacenado es par. Por
  esta última propiedad, a esta regla también se le llama **redondeo
  al par** (*round to even*).

El redondeo al más cercano es el más preciso, pero es un poco más
costoso de implementar correctamente. Algunos sistemas en el pasado
usaron reglas de redondeo más baratas de implementar, como el
truncamiento, pero el redondeo al más cercano es la regla de redondeo
por default en los sistemas del estándar IEEE.

**Ejemplo 1.6 — Reglas de redondeo.** Redondeando los siguientes
números decimales a dos dígitos con cada una de las reglas de
redondeo se obtienen los siguientes resultados:

| Número | Chop | Round to nearest |
|:------:|:----:|:-----------------:|
| 1.649  | 1.6  | 1.6                |
| 1.650  | 1.6  | 1.6                |
| 1.651  | 1.6  | 1.7                |
| 1.699  | 1.6  | 1.7                |
| 1.749  | 1.7  | 1.7                |
| 1.750  | 1.7  | 1.8                |
| 1.751  | 1.7  | 1.8                |
| 1.799  | 1.7  | 1.8                |

Nótese el caso de empate en 1.650 y 1.750: en 1.650, el número de
punto flotante más cercano con dos dígitos y último dígito par es
1.6 (no 1.7); en 1.750, es 1.8 (no 1.7). Esto ilustra la regla del
"par más cercano" cuando $x$ queda exactamente a la mitad entre dos
números de punto flotante representables.

### Precisión de máquina (Machine Precision)

La exactitud de un sistema de punto flotante puede caracterizarse con
una cantidad conocida como el **unit roundoff**, **machine
precision** o **machine epsilon**, que denotamos $\epsilon_{\text{mach}}$.
Su valor depende de la regla de redondeo usada. Con redondeo por
truncamiento (chop),

```math
\epsilon_{\text{mach}} = \beta^{1-t},
```

mientras que con redondeo al más cercano (round to nearest),

```math
\epsilon_{\text{mach}} = \tfrac{1}{2}\beta^{1-t}.
```

El unit roundoff es importante porque determina el máximo error
relativo posible al representar un número real $x$ distinto de cero en
un sistema de punto flotante:

```math
\left|\frac{fl(x) - x}{x}\right| \le \epsilon_{\text{mach}}.
```

Una caracterización alternativa del unit roundoff, que a veces se ve,
es que es el número $\epsilon$ más pequeño tal que

```math
fl(1 + \epsilon) > 1,
```

aunque esto no es del todo equivalente a la definición anterior si se
usa la regla de redondeo al par. Otra definición que a veces se usa es
que $\epsilon_{\text{mach}}$ es la distancia de 1 al siguiente número de
punto flotante más grande, pero esto puede diferir de las otras dos
definiciones. Aunque difieren en detalle, las tres definiciones de
$\epsilon_{\text{mach}}$ tienen la misma intención básica: medir la
granularidad de un sistema de punto flotante.

Para un sistema de punto flotante "de juguete" ilustrativo (como el
del Ejemplo 1.5 del texto, con $\beta = 2$, $t = 3$),
$\epsilon_{\text{mach}} = 0.25$ con redondeo por truncamiento, y
$\epsilon_{\text{mach}} = 0.125$ con redondeo al más cercano. Para los
sistemas de punto flotante binarios IEEE,
$\epsilon_{\text{mach}} = 2^{-24} \approx 10^{-7}$ en precisión simple
y $\epsilon_{\text{mach}} = 2^{-53} \approx 10^{-16}$ en precisión
doble. Por esto decimos que los sistemas de punto flotante IEEE de
precisión simple y doble tienen aproximadamente 7 y 16 dígitos
decimales de precisión, respectivamente.

Aunque ambas son "pequeñas", el unit roundoff **no debe confundirse**
con el underflow level. El unit roundoff $\epsilon_{\text{mach}}$ está
determinado por el número de dígitos del campo de la mantisa de un
sistema de punto flotante, mientras que el underflow level UFL está
determinado por el número de dígitos del campo del exponente. En todo
sistema de punto flotante práctico,

```math
0 < \text{UFL} < \epsilon_{\text{mach}} < \text{OFL}.
```

En Python, `sys.float_info.epsilon` da el valor de
$\epsilon_{\text{mach}}$ con redondeo al más cercano en doble
precisión ($2^{-52}$, la distancia de 1 al siguiente flotante más
grande — la tercera definición de arriba); el paquete del curso lo
expone en [`fiscomp.precision_numerica.EPS`](../../fiscomp/precision_numerica.py).

### Deducción de $\epsilon_{\text{mach}}$

Tomemos $x$ normalizado, positivo y con exponente $e$ (es decir,
$1 \le m < \beta$, así que $\beta^e \le x < \beta^{e+1}$):

```math
x = \left(d_0 + \frac{d_1}{\beta} + \cdots + \frac{d_{t-1}}{\beta^{t-1}} + \frac{d_t}{\beta^t} + \cdots\right)\beta^e.
```

**Truncamiento (chop).** $fl(x)$ conserva solo los primeros $t$
dígitos ($d_0, \dots, d_{t-1}$) y descarta el resto. El error absoluto
es lo que queda descartado:

```math
|x - fl(x)| = \left(\frac{d_t}{\beta^t} + \frac{d_{t+1}}{\beta^{t+1}} + \cdots\right)\beta^e
\le (\beta - 1)\sum_{i=t}^{\infty}\beta^{-i}\,\beta^e = \beta^{1-t}\,\beta^e.
```

(la suma geométrica $\sum_{i=t}^{\infty}\beta^{-i} = \beta^{1-t}/(\beta-1)$
es la misma cuenta que usamos para OFL, pero sumando desde $i=t$ en
vez de hasta $t-1$). Como $x \ge \beta^e$ por ser normalizado, el
error relativo queda acotado por

```math
\frac{|x - fl(x)|}{x} \le \frac{\beta^{1-t}\,\beta^e}{\beta^e} = \beta^{1-t} = \epsilon_{\text{mach}}.
```

**Redondeo al más cercano.** $\beta^{1-t}\beta^e$ es exactamente el
tamaño del "hueco" entre dos números de punto flotante consecutivos
con exponente $e$ (la distancia de un dígito en la última posición,
$d_{t-1}$). Redondear a la nearest solo puede acercarse, en el peor
caso, hasta la **mitad** de ese hueco (si $x$ cayera más lejos, ya
habría un flotante todavía más cercano). Entonces

```math
|x - fl(x)| \le \tfrac{1}{2}\,\beta^{1-t}\,\beta^e
\quad\Longrightarrow\quad
\frac{|x - fl(x)|}{x} \le \tfrac{1}{2}\beta^{1-t} = \epsilon_{\text{mach}}.
```

El factor $\tfrac{1}{2}$ es exactamente lo que diferencia las dos
fórmulas: redondear al más cercano nunca comete más de la mitad del
error que truncar, porque en vez de siempre "quedarse corto" hacia
cero, se elige el flotante más próximo en cualquier dirección.

**Por qué el peor caso ocurre en $m = 1$.** La cota
$|x - fl(x)|/x \le \epsilon_{\text{mach}}$ es más ajustada (se acerca
más a igualdad) cuando $x$ es lo más pequeño posible para ese
exponente, es decir $x = \beta^e$ (mantisa $m=1$): ahí el error
absoluto del último dígito pesa lo más posible en relación con $x$.
Para mantisas más grandes ($m \to \beta$), el mismo error absoluto es
una fracción relativa menor. Por eso $\epsilon_{\text{mach}}$ no
depende de $e$: al ser un cociente, el factor $\beta^e$ se cancela, y
solo sobrevive la dependencia en $\beta$ y $t$ (cuántos dígitos hay en
la mantisa).

**Verificación con el sistema de juguete** ($\beta = 2$, $t = 3$):

```math
\epsilon_{\text{mach}}^{\text{chop}} = 2^{1-3} = 2^{-2} = 0.25,
\qquad
\epsilon_{\text{mach}}^{\text{nearest}} = \tfrac{1}{2}\cdot 2^{-2} = 2^{-3} = 0.125,
```

que son justo los valores mencionados arriba.

## Error de redondeo acumulado en sumas largas

Cada suma de punto flotante $s_n = fl(s_{n-1} + x_n)$ introduce un
nuevo error de redondeo, acotado (relativo a $s_n$) por
$\epsilon_{\text{mach}}$. Al sumar $n$ términos uno tras otro (suma
naive, "corriendo" un acumulador), esos errores se van arrastrando:
el resultado calculado $\hat{s}_n$ satisface, en el peor caso,

```math
|\hat{s}_n - s_n| \lesssim (n - 1)\,\epsilon_{\text{mach}} \sum_{i=1}^{n} |x_i|.
```

Es decir, el error **crece con $n$**: no es que cada suma individual
sea muy imprecisa (cada una comete a lo más $\epsilon_{\text{mach}}$
de error relativo), sino que sumar muchos términos da muchas
oportunidades de acumular ese error. Esto es distinto del error de un
solo redondeo (como al leer un dato o hacer una sola operación): aquí
el error crece con el número de operaciones.

### El orden de la suma importa

Cuando el acumulador $s_{n-1}$ ya es mucho más grande que el
siguiente término $x_n$, sumar $x_n$ pierde precisión: $x_n$ se tiene
que representar a la escala de $s_{n-1}$, y sus dígitos de orden
más chico simplemente no caben en la mantisa (efectivamente lo mismo
que la cancelación catastrófica, pero al revés: en vez de restar dos
números casi iguales, se suma un número muy chico a uno muy grande y
la mayor parte de su información se descarta).

Por eso, sumar una lista de números positivos **de menor a mayor**
suele dar un resultado más preciso que sumarla de mayor a menor: así
el acumulador se mantiene comparable en tamaño a los términos que le
siguen durante más tiempo, en vez de volverse gigante desde el
principio y "tragarse" los términos chicos que vienen después.

### Suma compensada (algoritmo de Kahan)

El **algoritmo de suma de Kahan** (*Kahan summation* o *compensated
summation*) reduce el error acumulado sin cambiar el orden de los
términos ni la cantidad de operaciones aritméticas de fondo: en cada
paso guarda, en una variable aparte (`compensacion`), una estimación
de los dígitos de bajo orden que se perdieron en la suma anterior, y
la resta del siguiente término antes de sumarlo — así el error que se
iba a perder se "reintroduce" en la siguiente suma en vez de
descartarse. Con esto, el error acumulado deja de crecer con $n$ y
se queda del orden de $\epsilon_{\text{mach}}$, prácticamente
independiente del número de términos.

## Contenido
- [`aritmetica_punto_flotante.py`](aritmetica_punto_flotante.py): script
  principal con los ejemplos de esta unidad —
  underflow (UFL) y la zona de subnormales, overflow (OFL) y por qué
  `*` regresa `inf` pero `**` y `math.exp` lanzan `OverflowError`,
  dos ejemplos de cancelación catastrófica (raíces de una cuadrática
  con $b^2 \gg 4ac$, y $1 - \cos(x)$ para $x$ chica), error de
  redondeo acumulado en sumas largas (suma naive vs. `math.fsum`, el
  orden de la suma en la serie armónica, y la suma compensada de
  Kahan), y un ejercicio con la derivada numérica que muestra el
  error de redondeo compitiendo contra el error de truncamiento.
