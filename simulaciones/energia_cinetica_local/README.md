# Energía cinética local en mecánica cuántica

Proyecto que usa la segunda derivada con diferencias finitas
([unidad 08](../../unidades/08_diferencias_finitas/notas.md), sección
"Segunda derivada") para calcular la energía cinética de una partícula
cuántica en una dimensión. Si todavía no han llevado mecánica cuántica,
pueden leer la física por encima: la parte de programación (una misma
función que sirve para sistemas físicos distintos) vale por sí sola.

## La física

La ecuación de Schrödinger independiente del tiempo es un problema de
eigenvalores, $\hat{H}\psi = E\psi$, con $\hat{H} = \hat{T} + \hat{V}$
(energía cinética más potencial). En una dimensión, la energía cinética
es una segunda derivada:

$$
\hat{T}\psi = -\frac{\hbar^2}{2m}\frac{d^2\psi}{dx^2}.
$$

La **energía cinética local** se define como

$$
T_L(x) = \frac{\hat{T}\psi(x)}{\psi(x)} = -\frac{\hbar^2}{2m}\frac{\psi''(x)}{\psi(x)}.
$$

No es el valor esperado de la energía cinética: es su valor en cada
punto $x$, así que en general depende de $x$. La normalización de
$\psi$ no importa, porque se cancela entre numerador y denominador.
Aquí la calculamos para dos casos en los que conocemos la respuesta
exacta, para comprobar el método. Pero el método sirve igual para
cualquier $\psi$, aunque no sea eigenfunción o no se conozca
analíticamente; por ejemplo, en cálculos de Monte Carlo cuántico.

En todo el proyecto usamos unidades con $\hbar = m = 1$.

### Oscilador armónico

Con potencial $V(x) = \frac{1}{2}m\omega^2x^2$, las eigenenergías son
$E_n = (n + \tfrac{1}{2})\hbar\omega$, con $n = 0, 1, 2, \dots$ (la
energía del estado base no es cero: es la energía de punto cero). Las
eigenfunciones normalizadas son

$$
\psi_n(x) = \frac{1}{\sqrt{2^n n!}}\left(\frac{m\omega}{\pi\hbar}\right)^{1/4}
H_n\!\left(\sqrt{\frac{m\omega}{\hbar}}\,x\right)e^{-m\omega x^2/(2\hbar)},
$$

donde $H_n$ son los polinomios de Hermite. Se calculan con la relación
de recurrencia

$$
H_{j+1}(x) = 2xH_j(x) - 2jH_{j-1}(x),
\qquad H_0(x) = 1,\quad H_1(x) = 2x.
$$

Como $\hat{T}\psi + \hat{V}\psi = E\psi$, la energía cinética local
exacta es

$$
T_L(x) = E_n - V(x) = \left(n + \tfrac{1}{2}\right)\hbar\omega - \tfrac{1}{2}m\omega^2x^2.
$$

Para $|x|$ suficientemente grande, $T_L < 0$: la partícula cuántica
puede estar en la región **clásicamente prohibida**, algo imposible en
mecánica clásica.

**Principio de correspondencia.** Un oscilador clásico de amplitud
$x_0$ pasa más tiempo cerca de los puntos de retorno $\pm x_0$, donde
es más lento. Su densidad de probabilidad es

$$
P_c(x) = \frac{1}{\pi\sqrt{x_0^2 - x^2}}, \qquad |x| < x_0.
$$

Igualando su energía $\tfrac{1}{2}m\omega^2x_0^2$ con $E_n$ se obtiene
$x_0 = \sqrt{(2n+1)\hbar/(m\omega)}$. Para $n$ grande, la densidad
cuántica $|\psi_n(x)|^2$ oscila muy rápido; si se promedia localmente
(suavizando las oscilaciones), se parece a la clásica.

### Partícula en una caja periódica

Una partícula libre ($V = 0$) en una caja de longitud $L$, de
$-L/2$ a $L/2$, con condiciones periódicas: $\psi(x + L) = \psi(x)$.
Las eigenfunciones son ondas planas,

$$
\psi_k(x) = \frac{1}{\sqrt{L}}e^{ikx}, \qquad k = \frac{2\pi}{L}n,
\quad n = 0, \pm1, \pm2, \dots
$$

y la energía, que aquí es toda cinética y no depende de $x$, es

$$
E = \frac{\hbar^2k^2}{2m} = \frac{\hbar^2}{2m}\left(\frac{2\pi}{L}\right)^2n^2.
$$

## La programación: una interfaz común

Las dos funciones de onda de
[`funciones_de_onda.py`](funciones_de_onda.py) tienen la **misma
interfaz**: reciben la posición `x` y un diccionario `parametros`,
aunque cada una necesita parámetros distintos:

```python
psi_oscilador(x, {"n": 100, "m_w_hbar": 1.0, "alfa": 1.0})
psi_caja(x, {"n": -2, "L": 2 * PI})  # PI de fiscomp.funciones_especiales
```

Gracias a eso, `energia_cinetica_local(psi, x, parametros, h)` (en
[`energia_cinetica_local.py`](energia_cinetica_local.py)) funciona con
cualquiera de las dos, y con cualquier otra función de onda que siga la
misma interfaz, sin saber cuál es ni qué parámetros usa: solo le pasa
el diccionario tal cual a `psi`.

La segunda derivada se calcula con `diff2_central` (de
`fiscomp/derivadas.py`), la diferencia central de la unidad 08:

$$
\psi''(x) \approx \frac{4\left[\psi(x+h/2) + \psi(x-h/2) - 2\psi(x)\right]}{h^2}.
$$

`diff2_central` recibe una función de una sola variable, `f(x)`, pero
nuestras funciones de onda reciben también `parametros`. Por eso
`energia_cinetica_local` define adentro una pequeña función
`psi_de_x(y)` que llama a `psi(y, parametros)`: un **adaptador** entre
las dos interfaces.

## Archivos

- [`funciones_de_onda.py`](funciones_de_onda.py): `hermite`,
  `psi_oscilador` y `psi_caja`.
- [`energia_cinetica_local.py`](energia_cinetica_local.py):
  `energia_cinetica_local`, los resultados exactos y una prueba que
  compara los dos para varios `h`. También guarda los datos de las
  gráficas en `datos/`.
- [`graficar_energia_cinetica.gp`](graficar_energia_cinetica.gp):
  $T_L(x)$ numérica y exacta del oscilador con $n = 3$.
- [`graficar_densidad.gp`](graficar_densidad.gp): densidad cuántica
  contra clásica para $n = 100$ (principio de correspondencia).

## Cómo correrlo

Con el entorno virtual activado (ver
[`recursos/notas_entorno_virtual.md`](../../recursos/notas_entorno_virtual.md)),
desde esta carpeta:

```bash
python3 energia_cinetica_local.py
gnuplot graficar_energia_cinetica.gp
gnuplot graficar_densidad.gp
```

La tabla que imprime compara, en $x = 1$, la $T_L$ numérica contra la
exacta: 100 para el oscilador con $n = 100$, y 2 para la caja con
$n = -2$, $L = 2\pi$. Al reducir $h$ el error baja como $O(h^2)$,
hasta que el redondeo lo alcanza.

## Para explorar

1. **Escala de $h$.** Para el oscilador con $n = 100$ el error con
   $h = 0.1$ es mucho mayor que para la caja. ¿Por qué? (Vean la
   gráfica de la densidad: ¿qué tan separados están los nodos de
   $\psi_{100}$?) ¿Qué pasa con $n = 3$?
2. **La parte imaginaria.** La caja da un número complejo, aunque la
   energía es real. ¿Cómo cambia la parte imaginaria al reducir $h$?
   ¿Por qué crece justo cuando la parte real mejora?
3. **Nodos.** $T_L$ no está definida donde $\psi(x) = 0$. Calculen
   $T_L$ muy cerca de un nodo (por ejemplo, $x = 10^{-6}$ con $n = 3$).
   ¿Sigue coincidiendo con la exacta? ¿Por qué el script evita evaluar
   justo en $x = 0$?
4. **Función de prueba.** Con `"alfa"` distinto de 1, `psi_oscilador`
   ya no es una eigenfunción. Grafiquen $T_L(x) + V(x)$ para
   `alfa = 1` y para `alfa = 0.8`. ¿Qué cambia? (Para una
   eigenfunción, $T_L + V$ es la constante $E_n$ en todo $x$.)
5. **Una función de onda nueva.** Escriban otra función con la misma
   interfaz; por ejemplo, la de una partícula en una caja con paredes
   infinitas, $\psi_n(x) = \sqrt{2/L}\,\sin(n\pi x/L)$ para
   $0 < x < L$. Úsenla con `energia_cinetica_local` sin modificar esa
   función.
