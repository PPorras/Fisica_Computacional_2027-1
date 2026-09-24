# Notas de gnuplot

**gnuplot** es un programa de línea de comandos para hacer gráficas. Es
ligero, rápido y muy común en física computacional: un programa (en
Python, C, Fortran…) escribe resultados a un archivo de texto y gnuplot
los grafica. También sirve para graficar funciones directamente.

Página oficial: <http://www.gnuplot.info/>

## 1. Instalación

### Linux

```bash
sudo dnf install gnuplot      # Fedora / RHEL
sudo apt install gnuplot      # Ubuntu / Debian
sudo pacman -S gnuplot        # Arch
```

En Ubuntu/Debian, si al graficar no se abre ninguna ventana, instala la
versión con soporte gráfico: `sudo apt install gnuplot-qt`.

### macOS

Con [Homebrew](https://brew.sh/) (si no lo tienes, instálalo primero
siguiendo las instrucciones de su página):

```bash
brew install gnuplot
```

### Windows

Hay dos opciones:

1. **Instalador oficial.** Descarga el instalador `.exe` más reciente
   desde <https://sourceforge.net/projects/gnuplot/files/gnuplot/>
   (elige la carpeta de la versión más nueva y el archivo que termina
   en `-win64-mingw.exe`). Durante la instalación marca la opción
   **"Add application directory to your PATH environment variable"**
   para poder usar `gnuplot` desde PowerShell o cmd.exe.
2. **Con winget** (desde PowerShell):

   ```powershell
   winget install gnuplot.gnuplot
   ```

Si usas WSL (Linux dentro de Windows), sigue las instrucciones de Linux.

### Verificar la instalación

Abre una terminal nueva y corre:

```bash
gnuplot --version
```

Debe aparecer algo como `gnuplot 6.0 patchlevel 2`.

## 2. Primeros pasos: modo interactivo

Escribe `gnuplot` en la terminal. El prompt cambia a `gnuplot>` y ahí
escribes comandos. Para salir: `exit` o `quit` (o `Ctrl+D`).

```gnuplot
gnuplot> plot sin(x)
```

Se abre una ventana con la gráfica de sen(x). Algunos ejemplos más:

```gnuplot
plot sin(x), cos(x)              # dos funciones en la misma gráfica
plot x**2 - 3*x + 1              # potencias se escriben con **
plot [0:2*pi] sin(x)             # rango de x de 0 a 2π
plot [-5:5] [-1:30] x**2         # rango de x y de y
plot exp(-x**2) title "gaussiana"  # texto de la leyenda
```

### Funciones y variables propias

```gnuplot
g = 9.81
v0 = 20
altura(t) = v0*t - 0.5*g*t**2
plot [0:4.1] altura(x) title "tiro vertical"
```

> **Cuidado con la división entera:** en gnuplot, `1/2` vale `0`
> (igual que en C). Escribe `1.0/2` o `0.5` para obtener `0.5`.

### Personalizar la gráfica

Los comandos `set` cambian propiedades; afectan a los `plot` que vengan
**después**. Si ya graficaste, usa `replot` para volver a dibujar.

```gnuplot
set title "Oscilador armónico"
set xlabel "t (s)"
set ylabel "x (m)"
set grid                         # cuadrícula
set key top left                 # posición de la leyenda
set xrange [0:10]
set yrange [-1.5:1.5]
set samples 500                  # más puntos → curvas más suaves
plot cos(x) with lines linewidth 2 linecolor "blue" title "x(t)"
```

Muchas palabras tienen abreviaturas: `with lines` → `w l`,
`linewidth` → `lw`, `linecolor` → `lc`, `title` → `t`,
`with points` → `w p`, `pointtype` → `pt`.

Para quitar un ajuste: `unset grid`, `unset key`, etc. Para regresar
todo a los valores iniciales: `reset`.

### Ayuda

`help plot`, `help set xrange`, `help with`… gnuplot trae toda su
documentación integrada.

## 3. Graficar datos de un archivo

gnuplot lee archivos de texto con **columnas separadas por espacios**.
Las líneas que empiezan con `#` se ignoran (sirven para comentarios o
encabezados). Por ejemplo, `caida.dat`:

```text
# t(s)   y(m)    v(m/s)
0.0     100.00   0.00
0.5      98.77  -4.91
1.0      95.10  -9.81
1.5      88.96 -14.72
2.0      80.38 -19.62
2.5      69.34 -24.53
3.0      55.86 -29.43
```

Para graficarlo:

```gnuplot
plot "caida.dat"                          # usa columnas 1 (x) y 2 (y)
plot "caida.dat" using 1:2 with points    # lo mismo, explícito
plot "caida.dat" using 1:3 with lines     # t contra v
plot "caida.dat" u 1:2 w lp title "y(t)"  # lp = líneas y puntos
```

`using 1:2` significa "x = columna 1, y = columna 2". El archivo debe
estar en el directorio donde abriste gnuplot (o da la ruta completa).
Puedes ver y cambiar el directorio con `pwd` y `cd "carpeta"`.

### Operaciones con columnas

Dentro de `using`, `$2` significa "el valor de la columna 2", y se
puede operar con él (hay que poner la expresión entre paréntesis):

```gnuplot
plot "caida.dat" using 1:(-$3)            # rapidez (cambia el signo)
plot "caida.dat" using 1:(0.5*$3**2)      # energía cinética por unidad de masa
plot "caida.dat" using ($1*1000):2        # tiempo en milisegundos
```

### Varias series y comparación con la teoría

```gnuplot
plot "caida.dat" u 1:2 w p pt 7 t "simulación", \
     100 - 0.5*9.81*x**2 w l t "teoría"
```

La `\` al final de una línea indica que el comando continúa en la
siguiente.

### Barras de error

Si el archivo tiene una columna con la incertidumbre (por ejemplo, la
columna 3):

```gnuplot
plot "mediciones.dat" using 1:2:3 with yerrorbars
```

### Escalas logarítmicas

Útiles, por ejemplo, para ver cómo decrece un error numérico:

```gnuplot
set logscale y        # solo eje y
set logscale xy       # ambos ejes
unset logscale
```

### Ajuste de curvas (fit)

gnuplot puede ajustar parámetros por mínimos cuadrados:

```gnuplot
f(t) = y0 - 0.5*a*t**2
y0 = 100; a = 10                 # valores iniciales
fit f(x) "caida.dat" using 1:2 via y0, a
plot "caida.dat" u 1:2 w p t "datos", f(x) t "ajuste"
```

Al terminar, `fit` imprime los valores de `y0` y `a` con sus errores.

### Gráficas en 3D

Para archivos con tres columnas (x, y, z) se usa `splot`:

```gnuplot
splot "trayectoria.dat" using 1:2:3 with lines
splot sin(x)*cos(y)
```

## 4. Guardar la gráfica en un archivo

Por defecto gnuplot dibuja en una ventana. Para guardar una imagen hay
que cambiar la **terminal** (el formato de salida) y el archivo de
**salida**:

```gnuplot
set terminal pngcairo size 800,600
set output "caida.png"
replot                        # o el comando plot completo
set output                    # cierra el archivo (¡importante!)
set terminal qt               # regresa a la ventana (wxt en Windows,
                              # qt o aqua en macOS)
```

Otras terminales útiles: `pdfcairo` (PDF, ideal para reportes),
`svg` y `epslatex` (para LaTeX). Para ver cuáles tienes: `set terminal`.

## 5. Scripts de gnuplot

Escribir los comandos en un archivo es mucho mejor que teclearlos cada
vez: puedes corregirlos, reutilizarlos y guardarlos junto con tus datos.
Por convención los scripts terminan en `.gp` (también se ven `.gnu` o
`.plt`).

Ejemplo, `grafica_caida.gp`:

```gnuplot
# grafica_caida.gp — posición de un objeto en caída libre
reset

set terminal pngcairo size 800,600 font "Sans,12"
set output "caida.png"

set title "Caída libre desde 100 m"
set xlabel "t (s)"
set ylabel "y (m)"
set grid
set key top right

g = 9.81
y_teo(t) = 100 - 0.5*g*t**2

plot "caida.dat" using 1:2 with points pointtype 7 title "simulación", \
     y_teo(x) with lines linewidth 2 title "teoría"

set output
```

### Correr un script

Desde la terminal (en la carpeta donde están el script y los datos):

```bash
gnuplot grafica_caida.gp
```

Esto genera `caida.png` sin abrir ninguna ventana. Desde el modo
interactivo también se puede con:

```gnuplot
gnuplot> load "grafica_caida.gp"
```

### Script que muestra una ventana

Si el script grafica en pantalla (sin `set terminal`/`set output`), la
ventana se cierra en cuanto termina el script. Para mantenerla abierta:

```bash
gnuplot -persist grafica.gp
```

o agrega al final del script `pause -1 "Presiona Enter para salir"`.

### Pasar argumentos a un script

```bash
gnuplot -e "archivo='caida.dat'" grafica.gp
```

y dentro del script usa la variable: `plot archivo using 1:2`.

### Llamar gnuplot desde Python

Un flujo típico del curso: Python calcula y escribe los datos, y luego
llama al script de gnuplot:

```python
import subprocess

# ... escribir los resultados en caida.dat ...
subprocess.run(["gnuplot", "grafica_caida.gp"], check=True)
```

## 6. Resumen de comandos

| Comando                              | Qué hace                                  |
|--------------------------------------|-------------------------------------------|
| `plot f(x)`                          | Grafica una función                       |
| `plot "datos.dat" u 1:2 w l`         | Grafica columnas 1 y 2 de un archivo      |
| `splot`                              | Gráfica en 3D                             |
| `replot`                             | Vuelve a dibujar la última gráfica        |
| `set xrange [a:b]`                   | Rango del eje x                           |
| `set xlabel "..."` / `set title`     | Etiquetas y título                        |
| `set grid` / `set key`               | Cuadrícula / leyenda                      |
| `set logscale y`                     | Escala logarítmica                        |
| `fit f(x) "datos" u 1:2 via a,b`     | Ajuste por mínimos cuadrados              |
| `set terminal pngcairo`              | Salida a imagen PNG                       |
| `set output "archivo.png"`           | Nombre del archivo de salida              |
| `load "script.gp"`                   | Ejecuta un script desde gnuplot           |
| `reset`                              | Regresa los ajustes a sus valores iniciales |
| `help <tema>`                        | Ayuda integrada                           |
| `exit`                               | Salir                                     |
