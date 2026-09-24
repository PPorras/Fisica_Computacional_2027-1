#!/usr/bin/env gnuplot
# Grafica, en escala log-log, el error relativo de diff_forward y
# diff_central contra h, para ver las dos fuentes de error de
# notas.md compitiendo: el error de truncamiento baja con h
# (pendiente O(h) o O(h^2)) hasta que el error de redondeo lo detiene
# y lo hace volver a subir.
#
# Despues, si existe datos/derivada_sin_x2_precisas.dat (la parte de
# "diferencias finitas mas precisas" de notas.md), abre una segunda
# grafica que compara los cinco metodos: adelante O(h), central y
# segunda adelante O(h^2), segunda central y Richardson O(h^4).
#
# Uso (despues de completar la practica en diferencias_finitas.py --
# ver practica.md -- y correrlo, para que genere el .dat que se lee
# aqui):
#
#   gnuplot graficar_derivada.gp
#
# o, como el archivo es ejecutable y empieza con #!/usr/bin/env gnuplot:
#
#   ./graficar_derivada.gp

archivo = "datos/derivada_sin_x2.dat"

# Epsilon de la maquina en doble precision IEEE 754, escrito a mano
# (es un valor fijo del estandar, no depende de la corrida): la misma
# constante que fiscomp.precision_numerica.EPS.
EPS = 2.220446049250313e-16

# h_opt "de juguete" de notas.md (asumiendo f y sus derivadas de orden
# 1): el minimo teorico de cada curva de error, donde el error de
# truncamiento y el de redondeo quedan balanceados.
h_opt_adelante = sqrt(4.0 * EPS)
h_opt_central  = (24.0 * EPS)**(1.0 / 3.0)

set title "Error de diferencias finitas vs. tamano de paso h"
set xlabel "h"
set ylabel "error relativo"

# h y el error recorren varios ordenes de magnitud: log-log es la
# escala natural para ver las pendientes O(h) y O(h^2) como lineas
# rectas.
set logscale xy
set format x "10^{%L}"
set format y "10^{%L}"
set grid

set key top left

# Lineas verticales punteadas en cada h_opt teorico, con etiqueta,
# para comparar contra donde cada curva realmente da vuelta.
set arrow from h_opt_adelante, graph 0 to h_opt_adelante, graph 1 \
    nohead lc rgb "red" dt 2 lw 1.5
set arrow from h_opt_central, graph 0 to h_opt_central, graph 1 \
    nohead lc rgb "blue" dt 2 lw 1.5
set label "h_{opt} adelante" at h_opt_adelante, graph 0.05 \
    rotate by 90 offset 1,0 tc rgb "red"
set label "h_{opt} central" at h_opt_central, graph 0.05 \
    rotate by 90 offset 1,0 tc rgb "blue"

# Columnas del .dat: h  diff_forward  error_forward  diff_central  error_central
plot archivo using 1:3 with linespoints lc rgb "red"  pt 7 ps 0.5 title "adelante, O(h)", \
     archivo using 1:5 with linespoints lc rgb "blue" pt 7 ps 0.5 title "central, O(h^2)"

###############################################
# Segunda grafica: los cinco metodos
###############################################

# Solo si existe el archivo: `stats` intenta leerlo, y si lo logra
# define la variable STATS_records (el numero de lineas de datos).
archivo_precisas = "datos/derivada_sin_x2_precisas.dat"
stats archivo_precisas using 1 nooutput
if (!exists("STATS_records")) {
    pause -1 "Presiona Enter para cerrar..."
    exit
}

pause -1 "Presiona Enter para ver la comparacion de los cinco metodos..."

# h_opt "de juguete" de notas.md para los dos metodos nuevos (f y sus
# derivadas de orden 1).
h_opt_adelante2 = (48.0 * EPS)**(1.0 / 3.0)
h_opt_central2  = (1120.0 * EPS / 9.0)**(1.0 / 5.0)

set title "Diferencias finitas mas precisas: error vs. h"
unset arrow
unset label
set key top center

set arrow from h_opt_adelante2, graph 0 to h_opt_adelante2, graph 1 \
    nohead lc rgb "dark-orange" dt 2 lw 1.5
set arrow from h_opt_central2, graph 0 to h_opt_central2, graph 1 \
    nohead lc rgb "dark-green" dt 2 lw 1.5
set label "h_{opt} segunda adelante" at h_opt_adelante2, graph 0.05 \
    rotate by 90 offset 1,0 tc rgb "dark-orange"
set label "h_{opt} segunda central" at h_opt_central2, graph 0.05 \
    rotate by 90 offset 1,0 tc rgb "dark-green"

# Columnas: h  error_forward  error_central  error_forward2  error_central2  error_richardson
# Richardson va con linea punteada porque casi se encima con la
# segunda central (ambas son O(h^4)).
plot archivo_precisas using 1:2 with linespoints lc rgb "red"         pt 7 ps 0.5 title "adelante, O(h)", \
     archivo_precisas using 1:3 with linespoints lc rgb "blue"        pt 7 ps 0.5 title "central, O(h^2)", \
     archivo_precisas using 1:4 with linespoints lc rgb "dark-orange" pt 7 ps 0.5 title "segunda adelante, O(h^2)", \
     archivo_precisas using 1:5 with linespoints lc rgb "dark-green"  pt 7 ps 0.5 title "segunda central, O(h^4)", \
     archivo_precisas using 1:6 with lines       lc rgb "purple" dt 2 lw 2       title "Richardson, O(h^4)"

pause -1 "Presiona Enter para cerrar..."
