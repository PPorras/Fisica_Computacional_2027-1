# Práctica 1 — Tipos de datos, control de flujo y funciones

Cubre lo visto en las unidades
[00 (tipos de datos)](../unidades/00_tipos_de_datos/),
[01 (control de flujo y punto flotante)](../unidades/01_control_de_flujo_y_punto_flotante/),
[02 (funciones)](../unidades/02_funciones/),
[03 (módulos)](../unidades/03_modulos/)
y [04 (archivos)](../unidades/04_archivos/).

## Contexto

En clase implementamos, en [`fiscomp/funciones_especiales.py`](../fiscomp/funciones_especiales.py):

- `factorial(n)`: de forma iterativa.
- `seno(x)`: con la serie de Taylor alrededor de 0, usando `EPS`
  (de [`fiscomp/precision_numerica.py`](../fiscomp/precision_numerica.py))
  como criterio de convergencia — se suman términos de la serie
  mientras sigan siendo mayores o iguales que el épsilon de la
  máquina, y se corta la suma en cuanto un término deja de aportar
  precisión adicional.

## Ejercicio 1 — Carga del electrón (experimento de Millikan)

En el experimento de la gota de aceite de Millikan (1909) se rocían
gotas de aceite entre dos placas con un campo eléctrico, y se ajusta
el voltaje hasta que cada gota queda suspendida (la fuerza eléctrica
compensa exactamente su peso). De ahí se despeja la carga eléctrica
de cada gota. Lo interesante es que esas cargas medidas siempre
resultan (aproximadamente) múltiplos enteros de una carga mínima
común: la carga del electrón, `e ≈ 1.602176634 × 10⁻¹⁹ C`. El propio
Millikan estimó `e` así: dividiendo cada carga medida entre un número
entero de electrones, y promediando esas estimaciones — el mismo
método que van a programar en la Parte B.

### Parte A — Recolección de datos

Escriban un programa (`recoleccion_datos.py`, en esta misma carpeta)
que use la función `input()` para capturar los datos de una corrida
del experimento (real o inventada), y que use **todos los tipos de
datos vistos en la unidad 00**:

- `str` para el nombre del experimento y del responsable.
- `float`/`int` para las condiciones del experimento y las mediciones
  (recuerden convertir lo que regresa `input()`, que siempre es
  `str`).
- `bool` para alguna condición (por ejemplo, si la gota se considera
  válida — por ejemplo, si su carga medida es positiva y de un orden
  de magnitud razonable).
- una `tupla` para las condiciones del experimento que no deban
  cambiar entre gota y gota (por ejemplo, el voltaje aplicado entre
  las placas, la distancia entre ellas, la viscosidad del aceite).
- una `lista` con la carga medida de cada gota, en Coulombs (pidan al
  menos 3 o 4 gotas, para que el promedio de la Parte B tenga
  sentido).
- un `set` con los valores únicos medidos (por si dos gotas dieron la
  misma carga).
- un `dict` para el resumen final del experimento — lo van a terminar
  de llenar en la Parte B.

### Parte B — Estimar la carga del electrón

Escriban una función `estimar_carga_electron(cargas_medidas)` (en el
mismo script) que, a partir de la lista de cargas de la Parte A,
regrese su mejor estimación de `e` y la desviación estándar de esa
estimación. El método (el mismo que usó Millikan):

1. Tomen como primera aproximación de `e` la carga medida más chica
   (`min(cargas_medidas)`) — suponiendo que corresponde a una sola
   carga elemental.
2. Para cada gota, calculen cuántas cargas elementales tiene:
   `n = round(carga / e_aproximada)` (el entero más cercano).
3. Con ese `n`, calculen una estimación de `e` **por gota**:
   `carga / n`.
4. La estimación final de `e` es el **promedio** de esas estimaciones
   por gota; calculen también su **desviación estándar** (con
   `math.sqrt`, de la unidad 03) para tener una idea de qué tan
   dispersas quedaron.

Comparen su estimación final contra el valor aceptado,
`e = 1.602176634e-19 C`, usando `error_relativo()` (de
`fiscomp.precision_numerica`, que ya usaron en clase) y agréguenlo al
resumen.

### Reporte

Al final, el programa debe generar un archivo `reporte_recoleccion.txt`
(en esta misma carpeta) con un reporte que incluya el resumen del
`dict` (ahora con la carga estimada, la desviación estándar y el
error relativo contra el valor aceptado) y el detalle de cada gota
(carga medida, `n` y estimación individual), usando lo visto en la
[unidad 04](../unidades/04_archivos/).

> **Nota:** con pocas gotas (3-5) y datos inventados a mano, es fácil
> que la estimación quede sorprendentemente cerca del valor real; con
> datos experimentales de verdad entra ruido de medición y la
> desviación estándar sería mucho más grande.

## Ejercicio 2 — pi con el método de Leibniz

La serie de Leibniz aproxima pi con:

```
pi/4 = 1 - 1/3 + 1/5 - 1/7 + 1/9 - ...
```

Escriban un script que sume esta serie y calcule su mejor
aproximación de `pi`.

> **Aviso:** a diferencia de la serie de `seno()`, la serie de
> Leibniz converge muy lento (el error se reduce más o menos como
> `1/numero_de_terminos`). Llegar al épsilon de la máquina con esta
> serie tomaría un número de términos poco práctico de correr. Sumen
> tantos términos como les sea razonable (por ejemplo, unos cuantos
> millones) y, en vez de esperar llegar a `EPS`, **anoten en un
> comentario qué tan cerca quedaron** (comparando contra el `pi` "real"
> de la librería `math`) y por qué creen que no se puede llegar más
> lejos en un tiempo razonable.

Una vez calculado, guarden el resultado en un archivo de constantes
(`fiscomp/constantes.py`, con una línea como `PI = 3.14159...`), para
que el resto del código pueda hacer `from fiscomp.constantes import PI`
sin tener que volver a correr la suma cada vez que se use.

## Ejercicio 3 — El resto de las funciones especiales

Siguiendo el mismo patrón que `seno()` (series, con `EPS` como
criterio de corte, reutilizando `factorial()` donde haga falta),
agreguen a `fiscomp/funciones_especiales.py`:

- `coseno(x)`
- `exponencial(x)`
- `ln(x)`

**Pistas / cosas a investigar:**
- La serie de Taylor de `ln(x)` alrededor de `x = 1` solo converge
  para `0 < x <= 2`. Investiguen una serie alternativa (por ejemplo,
  la de `ln((1+y)/(1-y))` con `y = (x-1)/(x+1)`) si quieren que su
  `ln()` funcione para cualquier `x > 0`.
- No es obligatorio, pero si quieren un reto extra: investiguen qué
  es la "reducción de rango" (llevar `x` a un intervalo más chico,
  usando identidades trigonométricas, antes de sumar la serie) y por
  qué mejora la precisión de `seno()`/`coseno()` para `x` grande.

## Ejercicio 4 — Error de sus funciones especiales

Usando `error_relativo()` (de `fiscomp/precision_numerica.py`),
comparen sus funciones `seno`, `coseno`, `exponencial` y `ln` contra
las funciones "reales" del módulo `math` (`math.sin`, `math.cos`,
`math.exp`, `math.log`), para varios valores de `x`.

Guarden sus resultados en un archivo de texto (por ejemplo
`reporte_ejercicio4.txt`, en esta misma carpeta), usando lo visto en
la [unidad 04](../unidades/04_archivos/) (`open()` en modo `'w'`,
`write()`). Ahí mismo respondan: ¿hay algún valor de `x` donde el
error sea sorprendentemente alto? ¿Por qué creen que pasa (piensen en
qué tan cerca está el valor "real" de cero, y en qué le hace eso al
error *relativo*)?

## Cómo revisar su trabajo

Hay un script de pruebas en
[`pruebas_practica_01.py`](pruebas_practica_01.py) que revisa:

- Ejercicio 1: que `recoleccion_datos.py` corra de principio a fin
  sin crashear y que genere `reporte_recoleccion.txt` (no revisa qué
  tipos de datos usaron ni el contenido del reporte).
- Ejercicio 2: que el `PI` guardado en `fiscomp/constantes.py` esté
  cerca del `pi` real.
- Ejercicio 3: `coseno`, `exponencial` y `ln`, contra sus equivalentes
  de `math`.
- Ejercicio 4: que `reporte_ejercicio4.txt` exista y no esté vacío.

Los ejercicios 1 y 4 son reportes abiertos (no hay una única
respuesta "correcta"), así que esas pruebas solo revisan que el
script/archivo exista y no truene o esté vacío — no su contenido; esa
parte se las revisamos a mano. Córranlo desde la raíz del
repositorio:

```bash
python3 practicas/pruebas_practica_01.py
```

Por cada caso de prueba va a decir si pasó (`OK`) o no (`FALLÓ`), pero
**no** les dice cuál era el valor esperado ni el que obtuvieron — esa
parte les toca investigarla a ustedes. Si falla una de las funciones
(Ejercicio 2 o 3), revisen con `print()` qué está regresando para ese
valor de `x` y compárenlo a mano. Si falla el Ejercicio 1, el mensaje
de la prueba incluye el traceback con el error que tronó su script.
