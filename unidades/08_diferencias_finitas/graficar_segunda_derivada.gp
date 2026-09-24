#!/usr/bin/env gnuplot
# Grafica, en escala log-log, el error relativo de diff2_central (la
# diferencia central para la segunda derivada) contra h, para
# f(x) = sin(x^2) en x0 = 1. Se ve la misma competencia que en la
# primera derivada: el error de truncamiento baja como O(h^2) hasta
# que el error de redondeo (que ahora crece como 1/h^2) lo detiene.
#
# Uso (despues de correr segunda_derivada.py, que genera el .dat que
# se lee aqui):
#
#   gnuplot graficar_segunda_derivada.gp
#
# o, como el archivo es ejecutable y empieza con #!/usr/bin/env gnuplot:
#
#   ./graficar_segunda_derivada.gp

archivo = "datos/segunda_derivada_sin_x2.dat"

# Epsilon de la maquina en doble precision IEEE 754 (la misma
# constante que fiscomp.precision_numerica.EPS).
EPS = 2.220446049250313e-16

# h_opt "de juguete" de notas.md (f y f^(4) de orden 1).
h_opt = (768.0 * EPS)**(1.0 / 4.0)

set title "Segunda derivada: error vs. tamano de paso h"
set xlabel "h"
set ylabel "error relativo"

set logscale xy
set format x "10^{%L}"
set format y "10^{%L}"
set grid

set key bottom left

set arrow from h_opt, graph 0 to h_opt, graph 1 \
    nohead lc rgb "dark-green" dt 2 lw 1.5
set label "h_{opt}" at h_opt, graph 0.05 \
    rotate by 90 offset 1,0 tc rgb "dark-green"

# Columnas del .dat: h  diff2_central  error_relativo
plot archivo using 1:3 with linespoints lc rgb "dark-green" pt 7 ps 0.5 \
     title "diferencia central para f'', O(h^2)"

pause -1 "Presiona Enter para cerrar..."
