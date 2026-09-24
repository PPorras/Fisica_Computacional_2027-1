#!/usr/bin/env gnuplot
# Principio de correspondencia: densidad de probabilidad cuántica
# |psi_100(x)|^2 del oscilador armónico contra la densidad clásica
# P_c(x) = 1 / (pi sqrt(x0^2 - x^2)) de un oscilador con la misma
# energía.
#
# Uso (después de correr energia_cinetica_local.py, que genera el
# .dat que se lee aquí):
#
#   gnuplot graficar_densidad.gp

archivo = "datos/densidad_oscilador.dat"

set title "Oscilador armónico, n = 100: densidad cuántica vs. clásica"
set xlabel "x"
set ylabel "densidad de probabilidad"
set grid
set key top center
set yrange [0:0.2]

# La columna 3 (clásica) tiene "?" fuera de [-x0, x0], donde no está
# definida; gnuplot lo toma como dato faltante y no dibuja ahí.
set datafile missing "?"

plot archivo using 1:2 with lines lw 1 lc rgb "blue" title "cuántica, |psi_{100}(x)|^2", \
     archivo using 1:3 with lines lw 2 lc rgb "red" dt 2 title "clásica, P_c(x)"

pause -1 "Presiona Enter para cerrar..."
