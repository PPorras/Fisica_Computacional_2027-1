# Grafica el potencial electrostático del arreglo de 36 cargas
# (potencial_arreglo_cargas.py) como un mapa de calor, con curvas de
# nivel (equipotenciales) encima.
#
# Uso (después de correr potencial_arreglo_cargas.py, que genera el
# .dat que se lee aquí):
#
#   gnuplot graficar_potencial.gp

archivo = "datos/potencial_arreglo_cargas.dat"

set title "Potencial electrostático, arreglo de 36 cargas"
set xlabel "x"
set ylabel "y"
set size ratio -1  # misma escala en x y y, para no deformar el arreglo

set pm3d map
set palette defined (-1 "blue", 0 "white", 1 "red")  # rojo=positivo, azul=negativo
# El potencial diverge junto a cada carga; sin acotar el rango de
# color, esos picos dominan la escala y no se alcanza a ver el resto.
set cbrange [-2e2:2e2]
set cblabel "phi (V)"

set contour base
set cntrparam levels 15
unset clabel

splot archivo using 1:2:3 with pm3d nocontours notitle, \
      archivo using 1:2:3 with lines lc rgb "black" lw 0.5 nosurface notitle

pause -1 "Presiona Enter para cerrar..."
