# Fisica_Computacional_2027-1
Materiales, códigos y actividades del curso de Física Computacional, semestre 2027-1

## Descripción del curso

Este es un curso introductorio de métodos numéricos aplicados a la
física. Se empieza por el lenguaje (Python) y, desde el principio, por entender
cómo representa la computadora los números y de dónde salen los
errores. Con esas bases se construyen, paso a paso, los métodos
numéricos que se usan todos los días en física: derivadas, álgebra
lineal, búsqueda de raíces y mínimos, interpolación y ajuste de datos,
integrales y ecuaciones diferenciales.

La idea del curso es **no usar los métodos como cajas negras**: cada
método se programa desde cero, se compara contra resultados exactos,
se analiza su error y se aplica a un problema físico concreto (un
*proyecto* al final de cada tema). El código que se reutiliza entre
temas se va acumulando en el paquete [`fiscomp/`](fiscomp/).

El temario sigue, de manera general, el libro de Gezerlis (ver
[Bibliografía](#bibliografía)).

El repositorio está organizado así:

- [`unidades/`](unidades/) — notas (`notas.md`) y ejemplos ejecutables
  de cada tema.
- [`practicas/`](practicas/) — enunciados de las prácticas y sus
  archivos de pruebas.
- [`simulaciones/`](simulaciones/) — los proyectos de física de cada
  tema.
- [`fiscomp/`](fiscomp/) — paquete con el código reutilizable del
  curso.
- [`recursos/`](recursos/) — notas de referencia (entorno virtual,
  git, Linux, Vim, gnuplot, etc.).

## Temario

Los temas marcados como *(pendiente)* aún no se han visto en clase.

### 1. Python para física

- Tipos de datos — [unidad 00](unidades/00_tipos_de_datos/)
- Control de flujo — [unidad 01](unidades/01_control_de_flujo/)
- Funciones — [unidad 02](unidades/02_funciones/)
- Módulos — [unidad 03](unidades/03_modulos/)
- Archivos — [unidad 04](unidades/04_archivos/)
- Manejo de excepciones — [unidad 05](unidades/05_manejo_de_excepciones/)
- Programación orientada a objetos — [unidad 07](unidades/07_programacion_orientada_a_objetos/)
- Calidad de código: convenciones de nombres, docstrings, pruebas
  — [convenciones](recursos/convenciones_nombres_python.md)
- Gráficas con gnuplot y matplotlib — [gnuplot](recursos/notas_gnuplot.md),
  [matplotlib](recursos/notas_matplotlib.md)
- **Proyecto:** potencial electrostático de un arreglo de cargas
  — [`simulaciones/simulacion_potencial/`](simulaciones/simulacion_potencial/)
- Prácticas: [práctica 1](practicas/practica_01.md) (tipos de
  datos a archivos), [práctica 2](practicas/practica_02.md)
  (`VectorND` y `Matrix`)

### 2. Números

- Errores absoluto y relativo; error de redondeo contra error de
  truncamiento
- Representación de números reales: punto flotante (IEEE 754),
  precisión de máquina — [unidad 06](unidades/06_aritmetica_punto_flotante/)
- Errores de redondeo en la práctica: sumas largas, suma
  compensada (Kahan)
- Funciones especiales programadas desde cero con series de
  Taylor — [`fiscomp/funciones_especiales.py`](fiscomp/funciones_especiales.py)
- **Proyecto:** expansión multipolar en electromagnetismo *(pendiente)*

### 3. Derivadas

- Diferenciación analítica
- Diferencias finitas: hacia adelante, hacia atrás, central, de
  orden superior y segunda derivada — [unidad 08](unidades/08_diferencias_finitas/)
- Diferenciación automática *(pendiente)*
- **Proyecto:** energía cinética local en mecánica cuántica —
  [`simulaciones/energia_cinetica_local/`](simulaciones/energia_cinetica_local/)
- Práctica: [práctica 3](practicas/practica_03.md)

### 4. Matrices *(pendiente)*

- Análisis de error (número de condición)
- Solución de sistemas de ecuaciones lineales
- Problemas de eigenvalores
- Descomposición en valores singulares (SVD)
- **Proyecto:** la ecuación de Schrödinger como problema de
  eigenvalores

### 5. Ceros y mínimos *(pendiente)*

- Ecuaciones no lineales en una variable
- Ceros de polinomios
- Sistemas de ecuaciones no lineales
- Minimización en una y en varias dimensiones
- **Proyecto:** extremizar la acción en mecánica clásica

### 6. Aproximación *(pendiente)*

- Interpolación polinomial, con *splines* cúbicos y trigonométrica
- Ajuste lineal por mínimos cuadrados e inferencia estadística
- Ajuste no lineal por mínimos cuadrados
- **Proyecto:** comprobar la ley de Stefan–Boltzmann

### 7. Integrales *(pendiente)*

- Métodos de Newton–Cotes (trapecio, Simpson)
- Integración adaptativa y de Romberg
- Cuadratura gaussiana
- Monte Carlo
- **Proyecto:** Monte Carlo cuántico variacional

### 8. Ecuaciones diferenciales *(pendiente)*

- Problemas de valor inicial
- Problemas de valores en la frontera
- Problemas de eigenvalores
- Ecuaciones diferenciales parciales
- **Proyecto:** la ecuación de Poisson en dos dimensiones

## Bibliografía


- A. Gezerlis, *Numerical Methods in Physics with Python*, 2.ª ed.,
  Cambridge University Press, 2023.
- P. O. J. Scherer, *Computational Physics: Simulation of Classical
  and Quantum Systems*, Springer (Graduate Texts in Physics).
- T. Pang, *An Introduction to Computational Physics*, 2.ª ed.,
  Cambridge University Press, 2006.
- M. H. Holmes, *Introduction to Scientific Computing and Data
  Analysis*, 2.ª ed., Springer (Texts in Computational Science and
  Engineering), 2023.

## Instalación

El código reutilizable del curso vive en el paquete `fiscomp/` y se
instala en un entorno virtual (`.venv/`), en modo editable, para que
`import fiscomp` funcione desde cualquier práctica, unidad o
notebook:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -e .
```

Ver [`recursos/notas_entorno_virtual.md`](recursos/notas_entorno_virtual.md)
para más detalle (cómo activar/desactivar el entorno, verificar que
quedó bien instalado, etc.).

## Recursos

Notas de referencia para el curso, en [`recursos/`](recursos/):

- [Entorno virtual y el paquete `fiscomp`](recursos/notas_entorno_virtual.md)
- [Errores comunes](recursos/errores_comunes.md)
- [Comandos básicos de git](recursos/notas_git_basico.md)
- [Clonar un repositorio de GitHub](recursos/notas_clonar_github.md)
- [Comandos de Linux](recursos/notas_comandos_linux.md)
- [Convenciones de nombres en Python](recursos/convenciones_nombres_python.md)
- [Vim](recursos/notas_vim.md)
- [gnuplot](recursos/notas_gnuplot.md)
- [Instalar matplotlib](recursos/notas_matplotlib.md)
