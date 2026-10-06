#!/usr/bin/env gnuplot
# Grafica, en escala log-log, el tiempo de la sustitucion (O(n^2)) y
# del producto Matrix * Matrix (O(n^3)) contra n, junto con rectas de
# referencia proporcionales a n^2 y n^3. En log-log, t = C n^p es una
# recta de pendiente p.
#
# Uso (despues de completar complejidad.py -- ver la Practica 4 -- y
# correrlo, para que genere el .dat que se lee aqui):
#
#   gnuplot graficar_complejidad.gp

archivo = "datos/tiempos.dat"

# Primer renglon de datos (n0, t0) de cada columna, para que las
# rectas de referencia salgan del primer punto medido.
stats archivo using 1 every ::0::0 nooutput
n0 = STATS_min
stats archivo using 2 every ::0::0 nooutput
t0_sust = STATS_min
stats archivo using 3 every ::0::0 nooutput
t0_mm = STATS_min

set title "Tiempo de ejecucion vs. tamano del problema n"
set xlabel "n"
set ylabel "tiempo (s)"
set logscale xy
set format y "10^{%L}"
set grid
set key top left

# Columnas del .dat: n  t_sustitucion  t_mat_mat
plot archivo using 1:2 with linespoints lc rgb "blue" pt 7 title "sustitucion", \
     t0_sust * (x / n0)**2 with lines lc rgb "blue" dt 2 title "proporcional a n^2", \
     archivo using 1:3 with linespoints lc rgb "red" pt 7 title "Matrix * Matrix", \
     t0_mm * (x / n0)**3 with lines lc rgb "red" dt 2 title "proporcional a n^3"

pause -1 "Presiona Enter para cerrar..."
