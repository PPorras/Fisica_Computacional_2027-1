#!/usr/bin/env gnuplot
# Niveles de energia de dos espines 1/2 contra omega_I (Fig. 4.1 de [1],
# panel izquierdo), con omega_II = 2 omega_I y gamma = hbar = 1.
#
# Puntos: los eigenvalores del metodo QR. Lineas: los exactos, del
# bloque 2x2; cada linea sigue a un mismo estado, asi que se ve donde
# se cruzan dos niveles.
#
# Uso (despues de correr niveles_de_energia.py, que genera el .dat):
#
#   gnuplot graficar_dos_espines.gp
#
# [1] A. Gezerlis, Numerical Methods in Physics with Python, vol. 1.
#     Cambridge University Press, Cambridge, UK, 2023.

archivo = "datos/dos_espines.dat"

set title "Dos espines: niveles de energia"
set xlabel "omega_I / gamma"
set ylabel "E / gamma"
set grid
set key outside right

# En omega = 0 tres niveles coinciden (el triplete, espin total 1).
set label "triplete (E = 1/4)" at 0.1, 2.0 tc rgb "gray40"
set arrow from 0.1, 1.9 to 0.01, 0.35 lc rgb "gray40"

# Columnas: omega_I  E_0..E_3 (QR)  E_arriba_arriba  E_bloque_menor  E_bloque_mayor  E_abajo_abajo
plot archivo using 1:6 with lines lw 2 lc rgb "red"         title "arriba-arriba", \
     archivo using 1:7 with lines lw 2 lc rgb "blue"        title "bloque 2x2 (menor)", \
     archivo using 1:8 with lines lw 2 lc rgb "dark-green"  title "bloque 2x2 (mayor)", \
     archivo using 1:9 with lines lw 2 lc rgb "dark-orange" title "abajo-abajo", \
     for [k=2:5] archivo using 1:k every 4 with points pt 6 ps 0.8 lc rgb "black" \
         title (k == 2 ? "metodo QR" : "")

pause -1 "Presiona Enter para cerrar..."
