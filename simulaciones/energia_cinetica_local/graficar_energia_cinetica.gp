#!/usr/bin/env gnuplot
# Grafica la energía cinética local T_L(x) del oscilador armónico con
# n = 3: la calculada con diferencias finitas (puntos) y la exacta
# (línea).
#
# Uso (después de correr energia_cinetica_local.py, que genera el
# .dat que se lee aquí):
#
#   gnuplot graficar_energia_cinetica.gp

archivo = "datos/energia_cinetica_oscilador.dat"

set title "Energía cinética local, oscilador armónico (n = 3)"
set xlabel "x"
set ylabel "T_L(x)"
set grid
set key bottom center

# Puntos de retorno clásicos: donde T_L = 0, es decir x0 = sqrt(2n + 1).
x0 = sqrt(7.0)
set arrow from -x0, graph 0 to -x0, graph 1 nohead lc rgb "gray" dt 2
set arrow from  x0, graph 0 to  x0, graph 1 nohead lc rgb "gray" dt 2
set label "prohibida" at graph 0.99, graph 0.95 right tc rgb "gray"
set label "prohibida" at graph 0.01, graph 0.95 left tc rgb "gray"
set label "permitida" at 0, graph 0.6 center tc rgb "gray"

# T_L < 0 fuera de [-x0, x0]: algo imposible en mecánica clásica.
set arrow from graph 0, first 0 to graph 1, first 0 nohead lc rgb "black"

plot archivo using 1:3 with lines lw 2 lc rgb "blue" title "exacta: (n + 1/2) - x^2/2", \
     archivo using 1:2 every 8 with points pt 7 ps 0.6 lc rgb "red" title "diferencias finitas"

pause -1 "Presiona Enter para cerrar..."
