# Unidad 08 — Diferencias finitas

## Temas
- **Física computacional:** diferenciación numérica — diferencia hacia
  adelante (*forward difference*), diferencia hacia atrás (*backward
  difference*) y diferencia central (*central difference*); derivación
  de cada una a partir de la serie de Taylor; el error de truncamiento
  (*approximation error*) compitiendo contra el error de redondeo
  (*roundoff error*, ver unidad 06) al elegir el tamaño de paso
  $h$; fórmulas más precisas (segunda diferencia hacia adelante,
  $O(h^2)$, y segunda diferencia central, $O(h^4)$) a cambio de más
  evaluaciones de $f$; la segunda derivada con diferencia central.

## Diferenciación analítica

La derivada de $f$ en un punto $\tilde{x}$ se define, como siempre, por
el límite

$$
\left.\frac{df(x)}{dx}\right|_{\tilde{x}} = \lim_{h \to 0} \frac{f(\tilde{x}+h) - f(\tilde{x})}{h}.
$$

En la práctica casi nunca evaluamos este límite directamente: usamos
las reglas de diferenciación (potencias, cociente, producto, regla de
la cadena, derivadas de funciones especiales) que ya conocemos, por
ejemplo

$$
\frac{d}{dx}e^{\sin(2x)} = 2\cos(2x)\,e^{\sin(2x)}.
$$

Estas reglas dan la derivada exacta, como una expresión simbólica.
Cuando la evaluación numérica de esa expresión oscurece más de lo que
ayuda, o cuando el cálculo a mano es tedioso, se recurre a un paquete
de álgebra computacional (*computer algebra system*) como SymPy o Sage,
que manipula expresiones simbólicamente en vez de solo con números.

Pero en cómputo científico casi siempre queremos, en cambio, evaluar
derivadas *numéricamente*, a partir de la función evaluada en unos
cuantos puntos (por ejemplo, cuando la función viene de datos
experimentales, o cuando ni siquiera tenemos una expresión cerrada
para ella). Para eso sirven las **diferencias finitas**.

## Diferencias finitas

La idea ingenua es tomar la definición de la derivada y, en vez de
tomar el límite $h \to 0$, simplemente usar una $h$ "chica":

$$
f'(\tilde{x}) \approx \frac{f(\tilde{x}+h) - f(\tilde{x})}{h}.
$$

Esta receta ad hoc deja varias preguntas sin responder: ¿qué tan
precisa es? ¿qué tan "chica" debe ser $h$? Nótese además que, al hacer
$h$ más chica, el numerador $f(\tilde{x}+h) - f(\tilde{x})$ *también*
se hace más chico (estamos evaluando $f$ en dos puntos cada vez más
próximos entre sí): cualquier error al evaluar ese numerador se
amplifica al dividirlo entre una $h$ cada vez más chica. Para responder
estas preguntas de forma sistemática, en vez de partir de la
definición del límite, derivamos cada fórmula de diferencias finitas a
partir de la **serie de Taylor** de $f$, lo que nos da control expreso
sobre el error cometido.

### Diferencia hacia adelante (Forward difference)

Partimos de la expansión de Taylor de $f(x+h)$ alrededor de $x$:

$$
f(x+h) = f(x) + h f'(x) + \frac{h^2}{2}f''(x) + \frac{h^3}{6}f'''(x) + \cdots
$$

Despejando $f'(x)$:

$$
f'(x) = \frac{f(x+h) - f(x)}{h} - \frac{h}{2}f''(x) - \cdots
$$

y quedándonos solo con el primer término del lado derecho obtenemos la
**aproximación de diferencia hacia adelante**:

$$
f'(x) = \frac{f(x+h) - f(x)}{h} + O(h).
$$

Se llama "hacia adelante" porque arranca en $x$ y se mueve en la
dirección positiva hasta $x+h$: geométricamente, es la pendiente de la
recta que une $f(x)$ con $f(x+h)$. El error que se comete al truncar la
serie de Taylor es $O(h)$: es decir, el error decrece linealmente con
$h$ (a grandes rasgos, si $h$ no es demasiado chica, reducir $h$ a la
mitad reduce el error aproximadamente a la mitad).

### Diferencia hacia atrás (Backward difference)

De forma análoga, partiendo de la expansión de Taylor de $f(x-h)$:

$$
f(x-h) = f(x) - h f'(x) + \frac{h^2}{2}f''(x) - \frac{h^3}{6}f'''(x) + \cdots
$$

se despeja

$$
f'(x) = \frac{f(x) - f(x-h)}{h} + O(h),
$$

la **aproximación de diferencia hacia atrás**, que se mueve desde $x$
en dirección negativa hasta $x-h$. Forward y backward son ambas
**diferencias no centradas** (*non-central differences*): usan dos
evaluaciones de $f$ que no están distribuidas simétricamente alrededor
de $x$.

### Análisis de error de la diferencia hacia adelante

Queremos $h$ chica, pero no demasiado chica, porque hay dos fuentes de
error compitiendo:

- **Error de truncamiento** $E_{\text{app}}$, por cortar la serie de
  Taylor:
  $$
  E_{\text{app}} = \frac{h}{2}|f''(x)|.
  $$
  Este error *decrece* con $h$.
- **Error de redondeo** $E_{\text{ro}}$: el numerador
  $f(x+h) - f(x)$ resta dos números muy cercanos entre sí (ver
  cancelación catastrófica, unidad 06); el error absoluto de esa
  resta es aproximadamente $2|f(x)|\epsilon_{\text{mach}}$, y al
  dividir entre $h$:
  $$
  E_{\text{ro}} = \frac{2|f(x)|\epsilon_{\text{mach}}}{h}.
  $$
  Este error *crece* al achicar $h$.

El error total es $E = E_{\text{app}} + E_{\text{ro}}$. Minimizando
respecto a $h$ (derivando e igualando a cero) se obtiene el paso
óptimo y el error mínimo correspondiente:

$$
h_{\text{opt}} = \sqrt{\frac{4\epsilon_{\text{mach}}\,f(x)}{f''(x)}},
\qquad
E_{\text{opt}} = \sqrt{4\epsilon_{\text{mach}}\,f(x)\,f''(x)}.
$$

Si $f(x)$ y $f''(x)$ son de orden 1, con $\epsilon_{\text{mach}}
\approx 2^{-52} \approx 2\times10^{-16}$ (doble precisión) esto da
$h_{\text{opt}} \approx \sqrt[]{4\epsilon_{\text{mach}}} \approx
3\times10^{-8}$ y un error mínimo $E_{\text{opt}} \approx
3\times10^{-8}$ también. Ese error no es nada impresionante: la
diferenciación analítica (sección anterior) no tiene error de
truncamiento en absoluto, solo el error de evaluar $f$.

Ojo: afirmaciones como "el error de truncamiento es $O(h)$" solo valen
para funciones razonablemente bien portadas; si la derivada de orden
correspondiente no existe o diverge, esta escala no aplica.

### Diferencia central (Central difference)

Para mejorar sobre forward/backward, partimos de la expansión de
Taylor centrada en $x$, pero con pasos de tamaño $h/2$ hacia cada lado:

$$
f\!\left(x+\tfrac{h}{2}\right) = f(x) + \frac{h}{2}f'(x) + \frac{h^2}{8}f''(x) + \frac{h^3}{48}f'''(x) + \cdots
$$

$$
f\!\left(x-\tfrac{h}{2}\right) = f(x) - \frac{h}{2}f'(x) + \frac{h^2}{8}f''(x) - \frac{h^3}{48}f'''(x) + \cdots
$$

**Restando** la segunda de la primera, se cancelan $f(x)$ y todas las
derivadas de orden par; despejando $f'(x)$:

$$
f'(x) = \frac{f\!\left(x+\tfrac{h}{2}\right) - f\!\left(x-\tfrac{h}{2}\right)}{h} - \frac{h^2}{24}f'''(x) - \cdots
$$

lo que da la **aproximación de diferencia central**:

$$
f'(x) = \frac{f\!\left(x+\tfrac{h}{2}\right) - f\!\left(x-\tfrac{h}{2}\right)}{h} + O(h^2).
$$

Se llama "central" porque las dos evaluaciones, en $x - h/2$ y
$x + h/2$, están centradas alrededor de $x$ (y siguen separadas entre
sí por una distancia $h$, igual que forward/backward). También usa
solo dos evaluaciones de $f$, igual que forward, pero el error de
truncamiento es $O(h^2)$ en vez de $O(h)$: como $h$ es chica, $h^2$ es
todavía más chica, así que para la misma $h$ la diferencia central es
más precisa (a grandes rasgos, reducir $h$ a la mitad *cuadruplica* la
calidad de la aproximación, en vez de solo duplicarla).

Cuando solo se cuenta con una tabla de datos discretos
$(x_i, f(x_i))$, $i=0,\dots,n-1$, la diferencia central no se puede
usar en los extremos $x_0$ ni $x_{n-1}$ (no hay punto "del otro lado");
ahí es necesario usar forward o backward. En los puntos intermedios sí
se puede usar central, evaluando en $x_i \pm h/2$ (o, de forma
equivalente y más común en la práctica con datos en rejilla, usando
los vecinos $x_{i-1}$ y $x_{i+1}$).

### Análisis de error de la diferencia central

Con el mismo tipo de cuentas que en la sección anterior (ver problema
en el material del curso), el error total resulta

$$
E = E_{\text{app}} + E_{\text{ro}} = \frac{h^2}{24}|f'''(x)| + \frac{2|f(x)|\epsilon_{\text{mach}}}{h},
$$

y minimizando respecto a $h$:

$$
h_{\text{opt}} = \left(\frac{24\,\epsilon_{\text{mach}}\,f(x)}{f'''(x)}\right)^{1/3},
\qquad
E_{\text{opt}} = \left(\frac{9}{8}\epsilon_{\text{mach}}^2\,[f(x)]^2\,|f'''(x)|\right)^{1/3}.
$$

Con $f(x)$ y $f'''(x)$ de orden 1, esto da $h_{\text{opt}} \approx
(24\epsilon_{\text{mach}})^{1/3} \approx 2\times10^{-5}$ y
$E_{\text{opt}} \approx (9\epsilon_{\text{mach}}/8)^{1/3} \approx
4\times10^{-11}$: mucho mejor que el $10^{-8}$ de la diferencia hacia
adelante, aunque todavía peor que la diferenciación analítica.

Un resultado llamativo: el paso óptimo de la diferencia central
($10^{-5}$) es *mucho más grande* que el de la diferencia hacia
adelante ($10^{-8}$), a pesar de que el error final es mucho menor. El
mejor algoritmo permite usar una $h$ más grande y aun así obtener un
mejor resultado.

## Diferencias finitas más precisas

Hasta aquí, tanto la diferencia hacia adelante ($O(h)$) como la central
($O(h^2)$) usan **dos** evaluaciones de $f$. Si estamos dispuestos a
evaluar $f$ en más puntos, podemos combinar más series de Taylor para
cancelar más términos del error, y obtener fórmulas de orden más alto.

### Segunda diferencia hacia adelante

Usamos tres puntos, todos a la derecha de $x$ (o en $x$ mismo): $x$,
$x+h/2$ y $x+h$. Las series de Taylor de los dos últimos son

$$
f\!\left(x+\tfrac{h}{2}\right) = f(x) + \frac{h}{2}f'(x) + \frac{h^2}{8}f''(x) + \frac{h^3}{48}f'''(x) + \cdots
$$

$$
f(x+h) = f(x) + h f'(x) + \frac{h^2}{2}f''(x) + \frac{h^3}{6}f'''(x) + \cdots
$$

Queremos una combinación en la que se cancele el término de $f''$.
Multiplicando la primera por 4 y restándole la segunda, los términos
$\tfrac{4h^2}{8}f''$ y $\tfrac{h^2}{2}f''$ se cancelan; restando además
$3f(x)$ desaparece también $f(x)$:

$$
4f\!\left(x+\tfrac{h}{2}\right) - f(x+h) - 3f(x) = h f'(x) - \frac{h^3}{12}f'''(x) + \cdots
$$

Despejando $f'(x)$ se obtiene la **segunda diferencia hacia adelante**
(*second forward difference*):

$$
f'(x) = \frac{4f\!\left(x+\tfrac{h}{2}\right) - f(x+h) - 3f(x)}{h} + \frac{h^2}{12}f'''(x) + \cdots
$$

Sigue siendo "hacia adelante" (solo usa $f$ en $x$ y a su derecha),
pero su error es $O(h^2)$, como el de la diferencia central. El precio
es que necesita **tres** evaluaciones de $f$ en vez de dos. Es útil,
por ejemplo, en el extremo izquierdo de una tabla de datos, donde la
diferencia central no se puede usar.

Nótese que los coeficientes suman cero ($4 - 1 - 3 = 0$): tiene que
ser así, porque la derivada de una función constante debe salir cero.

**Análisis de error.** El error de truncamiento es
$E_{\text{app}} = \frac{h^2}{12}|f'''(x)|$. Para el de redondeo usamos
el mismo modelo que antes: cada evaluación de $f$ trae un error de
aproximadamente $|f(x)|\epsilon_{\text{mach}}$, multiplicado por su
coeficiente; con coeficientes $4$, $1$ y $3$ el error del numerador es
$(4+1+3)|f(x)|\epsilon_{\text{mach}} = 8|f(x)|\epsilon_{\text{mach}}$,
y al dividir entre $h$:

$$
E = \frac{h^2}{12}|f'''(x)| + \frac{8|f(x)|\epsilon_{\text{mach}}}{h}.
$$

Minimizando respecto a $h$:

$$
h_{\text{opt}} = \left(\frac{48\,\epsilon_{\text{mach}}\,f(x)}{f'''(x)}\right)^{1/3},
\qquad
E_{\text{opt}} = \left(36\,\epsilon_{\text{mach}}^2\,[f(x)]^2\,|f'''(x)|\right)^{1/3}.
$$

Con $f$ y $f'''$ de orden 1: $h_{\text{opt}} \approx 2\times10^{-5}$ y
$E_{\text{opt}} \approx 1\times10^{-10}$. Es decir, el mismo orden que
la diferencia central, un poco peor (por el redondeo de una evaluación
más) y con una evaluación extra de $f$.

### Segunda diferencia central

Ahora usamos cuatro puntos, simétricos alrededor de $x$:
$x \pm h/2$ y $x \pm 3h/2$. Para cualquier $a$, restando las series de
$f(x+a)$ y $f(x-a)$ se cancelan todas las derivadas pares:

$$
f(x+a) - f(x-a) = 2a f'(x) + \frac{a^3}{3}f'''(x) + \frac{a^5}{60}f^{(5)}(x) + \cdots
$$

Con $a = h/2$ y $a = 3h/2$:

$$
f\!\left(x+\tfrac{h}{2}\right) - f\!\left(x-\tfrac{h}{2}\right) = h f'(x) + \frac{h^3}{24}f'''(x) + \frac{h^5}{1920}f^{(5)}(x) + \cdots
$$

$$
f\!\left(x+\tfrac{3h}{2}\right) - f\!\left(x-\tfrac{3h}{2}\right) = 3h f'(x) + \frac{27h^3}{24}f'''(x) + \frac{243h^5}{1920}f^{(5)}(x) + \cdots
$$

Multiplicando la primera por 27 y restándole la segunda se cancela el
término de $f'''$ (el que daba el error $O(h^2)$ de la diferencia
central), y queda $24h f'(x) - \frac{216h^5}{1920}f^{(5)}(x)$.
Despejando se obtiene la **segunda diferencia central** (*second
central difference*):

$$
f'(x) = \frac{27f\!\left(x+\tfrac{h}{2}\right) + f\!\left(x-\tfrac{3h}{2}\right) - 27f\!\left(x-\tfrac{h}{2}\right) - f\!\left(x+\tfrac{3h}{2}\right)}{24h} + \frac{3}{640}h^4 f^{(5)}(x) + \cdots
$$

Su error es $O(h^4)$: es la fórmula más precisa de las que hemos
visto, a cambio de **cuatro** evaluaciones de $f$. Otra vez, los
coeficientes suman cero ($27 + 1 - 27 - 1 = 0$).

**Análisis de error.** Con coeficientes $27, 1, 27, 1$ el error de
redondeo del numerador es $56|f(x)|\epsilon_{\text{mach}}$, que
dividido entre $24h$ da $\frac{7|f(x)|\epsilon_{\text{mach}}}{3h}$:

$$
E = \frac{3}{640}h^4|f^{(5)}(x)| + \frac{7|f(x)|\epsilon_{\text{mach}}}{3h}.
$$

Minimizando respecto a $h$:

$$
h_{\text{opt}} = \left(\frac{1120\,\epsilon_{\text{mach}}\,f(x)}{9\,f^{(5)}(x)}\right)^{1/5},
\qquad
E_{\text{opt}} = \frac{3}{128}\,|f^{(5)}(x)|\,h_{\text{opt}}^4.
$$

Con $f$ y $f^{(5)}$ de orden 1: $h_{\text{opt}} \approx 2\times10^{-3}$
y $E_{\text{opt}} \approx 3\times10^{-13}$, unas 100 veces mejor que la
diferencia central. Y otra vez: el mejor método permite usar una $h$
todavía más grande.

### ¿Por qué no se usan siempre?

Se puede seguir así indefinidamente: usar más puntos, combinar más
series de Taylor y cancelar más términos, para obtener fórmulas
(centrales o no) de orden cada vez más alto. El problema es que cada
una necesita más evaluaciones de $f$, y en la práctica evaluar $f$
suele ser lo caro (por ejemplo, si cada evaluación es una simulación
completa). Por eso estas fórmulas de orden alto no son tan comunes en
la práctica como la diferencia central.

**Para pensar:** construyan otra fórmula $O(h^4)$ combinando la
diferencia central con paso $h$,
$\frac{f(x+h/2) - f(x-h/2)}{h}$, con la diferencia central con paso
$2h$, $\frac{f(x+h) - f(x-h)}{2h}$. (Pista: escriban el término de
error $O(h^2)$ de cada una y busquen la combinación que lo cancela.)

### Comparación

| Método    | Evaluaciones de $f$          | Error de truncamiento | $h_{\text{opt}}$ (orden) | $E_{\text{opt}}$ (orden) |
|:----------|:------------------------------|:----------------------:|:-------------------------:|:--------------------------:|
| Adelante  | $f(x)$, $f(x+h)$               | $O(h)$                  | $10^{-8}$                  | $10^{-8}$                    |
| Atrás     | $f(x-h)$, $f(x)$               | $O(h)$                  | $10^{-8}$                  | $10^{-8}$                    |
| Central   | $f(x-h/2)$, $f(x+h/2)$          | $O(h^2)$                | $10^{-5}$                  | $10^{-11}$                   |
| Segunda adelante | $f(x)$, $f(x+h/2)$, $f(x+h)$ | $O(h^2)$           | $10^{-5}$                  | $10^{-10}$                   |
| Segunda central  | $f(x\pm h/2)$, $f(x\pm 3h/2)$ | $O(h^4)$          | $10^{-3}$                  | $10^{-13}$                   |

(órdenes de magnitud para $\epsilon_{\text{mach}}$ de doble precisión y
$f$, sus derivadas, de orden 1)

## Segunda derivada

En física la segunda derivada aparece en todos lados: la aceleración
en la segunda ley de Newton, el laplaciano en la ecuación de Schrödinger
o en la de difusión, la curvatura de un potencial alrededor de un
mínimo... Así que también necesitamos aproximarla con diferencias
finitas.

Una forma de hacerlo es pensar que la segunda derivada es la derivada
de la primera derivada, y aplicar la diferencia central dos veces:

$$
f''(x) \approx \frac{f'\!\left(x+\tfrac{h}{2}\right) - f'\!\left(x-\tfrac{h}{2}\right)}{h}.
$$

Pero, como en las secciones anteriores, es más instructivo partir de
las series de Taylor, porque así obtenemos también el término de
error.

### Diferencia central para la segunda derivada

Volvemos a las series de $f(x+h/2)$ y $f(x-h/2)$ de la diferencia
central. Al **restarlas** se cancelaban las derivadas pares y quedaban
solo las impares (útil para $f'$). Si en cambio las **sumamos**, se
cancelan las derivadas impares y quedan solo las pares:

$$
f\!\left(x+\tfrac{h}{2}\right) + f\!\left(x-\tfrac{h}{2}\right) = 2f(x) + \frac{h^2}{4}f''(x) + \frac{h^4}{192}f^{(4)}(x) + \cdots
$$

Como buscamos $f''$, vamos por buen camino: despejando,

$$
f''(x) = \frac{4\left[f\!\left(x+\tfrac{h}{2}\right) + f\!\left(x-\tfrac{h}{2}\right) - 2f(x)\right]}{h^2} - \frac{h^2}{48}f^{(4)}(x) - \cdots
$$

lo que da la **aproximación de diferencia central para la segunda
derivada**:

$$
f''(x) = \frac{4\left[f\!\left(x+\tfrac{h}{2}\right) + f\!\left(x-\tfrac{h}{2}\right) - 2f(x)\right]}{h^2} + O(h^2).
$$

Usa tres evaluaciones de $f$, y otra vez los coeficientes suman cero
($1 + 1 - 2 = 0$): la segunda derivada de una constante debe salir
cero.

Si en vez de $x \pm h/2$ se usan los puntos $x \pm h$ (con $h$ el
espaciamiento de una rejilla), la misma fórmula toma la forma que
probablemente ya han visto:

$$
f''(x) \approx \frac{f(x+h) - 2f(x) + f(x-h)}{h^2}.
$$

Es exactamente lo que se obtiene al aplicar la diferencia central dos
veces, como al inicio de esta sección.

### Análisis de error

El error de truncamiento es $E_{\text{app}} = \frac{h^2}{48}|f^{(4)}(x)|$.
Para el de redondeo, los coeficientes del numerador son $1, 1, 2$, así
que su error es $(1+1+2)|f(x)|\epsilon_{\text{mach}} =
4|f(x)|\epsilon_{\text{mach}}$, y se multiplica por $4/h^2$:

$$
E = \frac{h^2}{48}|f^{(4)}(x)| + \frac{16|f(x)|\epsilon_{\text{mach}}}{h^2}.
$$

Minimizando respecto a $h$:

$$
h_{\text{opt}} = \left(\frac{768\,\epsilon_{\text{mach}}\,f(x)}{f^{(4)}(x)}\right)^{1/4},
\qquad
E_{\text{opt}} = \sqrt{\frac{4}{3}\epsilon_{\text{mach}}\,f(x)\,f^{(4)}(x)}.
$$

Con $f$ y $f^{(4)}$ de orden 1: $h_{\text{opt}} \approx 6\times10^{-4}$
y $E_{\text{opt}} \approx 2\times10^{-8}$.

Nótese que, aunque el error de truncamiento es $O(h^2)$ como en la
diferencia central de la primera derivada, el error mínimo es mucho
peor ($10^{-8}$ contra $10^{-11}$): ahora el error de redondeo se
divide entre $h^2$, no entre $h$, así que crece mucho más rápido al
achicar $h$. En general, **cada derivada adicional que se aproxima con
diferencias finitas cuesta varios dígitos de precisión**.

### Derivadas de orden más alto

En física casi nunca se necesitan derivadas más allá de la segunda,
pero la receta para obtenerlas ya está clara: las **sumas** de series
de Taylor dan derivadas pares, las **restas** dan derivadas impares, y
combinando suficientes sumas o restas (con distintos pasos) se pueden
cancelar todos los términos que no queremos.

**Para pensar:** aproximen $f^{(4)}(x)$ combinando
$f(x+h/2) + f(x-h/2)$ con $f(x+h) + f(x-h)$ (y $f(x)$). (Pista:
escriban las dos sumas hasta el término de $h^4$ y busquen la
combinación que cancela el término de $f''$.) ¿De qué orden es el
error?

## Contenido
- [`diferencias_finitas.py`](diferencias_finitas.py): punto de
  partida, con `diff_forward` ya implementada, las funciones de
  prueba ($f(x)=5$, $f(x)=x$, $f(x)=x^2$, $f(x)=\sin(x^2)$) y un
  ejemplo de uso. `diff_backward`, `diff_central`, la comparación
  contra la derivada exacta, el barrido de $h$ (potencias de $1/2$)
  que muestra en números la competencia entre el error de
  truncamiento y el error de redondeo descrita arriba, y guardar ese
  barrido en `datos/derivada_sin_x2.dat` son la **práctica** (ver
  [`practica.md`](practica.md)).
- [`graficar_derivada.gp`](graficar_derivada.gp) y
  [`graficar_derivada.py`](graficar_derivada.py) (gnuplot y
  matplotlib, respectivamente — este último con instrucciones de
  instalación en el `.venv` del curso al inicio del archivo): ya
  están completos, listos para usarse una vez terminada la práctica.
  Grafican, en escala log-log, el error relativo de las diferencias
  hacia adelante y central contra $h$ (leyendo
  `datos/derivada_sin_x2.dat`, con el formato de columnas descrito en
  `practica.md`); se ven las pendientes $O(h)$ y $O(h^2)$ del error de
  truncamiento, y cómo cada curva da vuelta y empieza a subir cuando
  el error de redondeo toma el control. Marcan con una línea vertical
  punteada el $h_{\text{opt}}$ "de juguete" de cada método (las
  fórmulas de la sección de análisis de error arriba), justo donde
  cada curva da vuelta. Si además existe
  `datos/derivada_sin_x2_precisas.dat` (los métodos de la sección
  "Diferencias finitas más precisas"), muestran una segunda gráfica
  que compara los cinco métodos, con pendientes $O(h)$, $O(h^2)$ y
  $O(h^4)$. Se corren con `gnuplot graficar_derivada.gp` o
  `python3 graficar_derivada.py`, después de haber completado la
  práctica en `diferencias_finitas.py` (que es quien genera el
  `.dat`).
- [`segunda_derivada.py`](segunda_derivada.py): la diferencia central
  para la segunda derivada (sección "Segunda derivada"), comparada
  contra la segunda derivada exacta de las funciones de prueba, con un
  barrido de $h$ guardado en `datos/segunda_derivada_sin_x2.dat`.
- [`graficar_segunda_derivada.gp`](graficar_segunda_derivada.gp) y
  [`graficar_segunda_derivada.py`](graficar_segunda_derivada.py):
  grafican ese barrido en escala log-log, igual que las gráficas de
  la primera derivada. Se corren después de `segunda_derivada.py`.
