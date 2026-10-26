#!/usr/bin/env gnuplot
# Niveles de energia de cuatro espines 1/2 contra omega_I (problema 4.50
# de [1]), con omega_k = k omega_I (k = 1, 2, 3, 4)
# y gamma = hbar = 1, calculados con el metodo QR.
#
# Se grafican puntos, no lineas: las energias de cada renglon estan
# ordenadas de menor a mayor, asi que una columna no sigue siempre al
# mismo estado (en un cruce, dos estados intercambian columna).
#
# Uso (despues de correr niveles_de_energia.py, que genera el .dat):
#
#   gnuplot graficar_cuatro_espines.gp
#
# [1] A. Gezerlis, Numerical Methods in Physics with Python, vol. 1.
#     Cambridge University Press, Cambridge, UK, 2023.

archivo = "datos/cuatro_espines.dat"

set title "Cuatro espines: niveles de energia"
set xlabel "omega_I / gamma"
set ylabel "E / gamma"
set grid
unset key

plot for [k=2:17] archivo using 1:k with points pt 7 ps 0.6 lc rgb "blue"

pause -1 "Presiona Enter para cerrar..."
