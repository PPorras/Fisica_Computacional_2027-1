#!/usr/bin/env gnuplot
# Niveles de energia de dos espines 1/2 contra omega_I (Fig. 4.1 del
# libro, panel izquierdo), con omega_II = 2 omega_I y gamma = hbar = 1.
#
# Uso (despues de correr niveles_de_energia.py, que genera el .dat):
#
#   gnuplot graficar_dos_espines.gp

archivo = "datos/dos_espines.dat"

set title "Dos espines: niveles de energia"
set xlabel "omega_I / gamma"
set ylabel "E / gamma"
set grid
set key bottom left

# En omega = 0 tres niveles coinciden (el triplete, espin total 1).
set label "triplete (E = 1/4)" at 0.05, 0.6 tc rgb "gray"

# Columnas: omega_I  E_arriba_arriba  E_bloque_menor  E_bloque_mayor  E_abajo_abajo
plot archivo using 1:2 with lines lw 2 lc rgb "red"         title "arriba-arriba", \
     archivo using 1:3 with lines lw 2 lc rgb "blue"        title "bloque 2x2 (menor)", \
     archivo using 1:4 with lines lw 2 lc rgb "dark-green"  title "bloque 2x2 (mayor)", \
     archivo using 1:5 with lines lw 2 lc rgb "dark-orange" title "abajo-abajo"

pause -1 "Presiona Enter para cerrar..."
