#!/usr/bin/env gnuplot
# Niveles de energia de tres espines 1/2 contra omega_I (Fig. 4.1 de
# [1], panel derecho), con omega_II = 2 omega_I, omega_III = 3 omega_I
# y gamma = hbar = 1, calculados con el metodo QR.
#
# Se grafican puntos, no lineas: las energias de cada renglon estan
# ordenadas de menor a mayor, asi que una columna no sigue siempre al
# mismo estado (en un cruce, dos estados intercambian columna).
#
# Uso (despues de correr niveles_de_energia.py, que genera el .dat):
#
#   gnuplot graficar_tres_espines.gp
#
# [1] A. Gezerlis, Numerical Methods in Physics with Python, vol. 1.
#     Cambridge University Press, Cambridge, UK, 2023.

archivo = "datos/tres_espines.dat"

set title "Tres espines: niveles de energia"
set xlabel "omega_I / gamma"
set ylabel "E / gamma"
set grid
unset key

plot for [k=2:9] archivo using 1:k with points pt 7 ps 0.6 lc rgb "blue"

pause -1 "Presiona Enter para cerrar..."
