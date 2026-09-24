#!/usr/bin/env python3
"""Energía cinética local en mecánica cuántica (ver README.md).

La energía cinética local es

    T_L(x) = (T psi)(x) / psi(x) = -(hbar^2 / 2m) psi''(x) / psi(x),

y la segunda derivada se aproxima con diff2_central de
fiscomp/derivadas.py, la diferencia central de la unidad 08
(unidades/08_diferencias_finitas/notas.md, sección "Segunda
derivada"). energia_cinetica_local() funciona con *cualquier* función
de onda que siga la interfaz psi(x, parametros) de
funciones_de_onda.py: aquí se usa con el oscilador armónico y con la
caja periódica, y se compara contra el resultado exacto.

Al final guarda dos archivos en datos/, para graficarlos:

- energia_cinetica_oscilador.dat: T_L(x) numérica y exacta del
  oscilador con n = 3 (graficar_energia_cinetica.gp).
- densidad_oscilador.dat: |psi_100(x)|^2 y la densidad de
  probabilidad clásica, para ver el principio de correspondencia
  (graficar_densidad.gp).

Unidades: hbar = m = 1.
"""

from pathlib import Path

from fiscomp.derivadas import diff2_central, h_optimo_segunda
from fiscomp.funciones_especiales import PI, raiz_cuadrada
from fiscomp.precision_numerica import error_relativo
from funciones_de_onda import psi_caja, psi_oscilador

CARPETA_DATOS = Path(__file__).resolve().parent / "datos"

# hbar^2 / m, en unidades con hbar = m = 1.
HBAR2_SOBRE_M = 1.0

# h óptimo "de juguete" de la diferencia central para la segunda
# derivada (unidad 08): (768 EPS)^(1/4) ~ 6e-4.
H_OPT = h_optimo_segunda()

###############################################
# Energía cinética local
###############################################


def energia_cinetica_local(psi, x: float, parametros: dict, h: float = H_OPT):
    """Energía cinética local T_L(x) = -(hbar^2 / 2m) psi''(x) / psi(x).

    psi: función de onda con la interfaz psi(x, parametros).
    x: posición (donde psi(x) != 0: en un nodo de psi, T_L no está
       definida).
    parametros: diccionario que se le pasa tal cual a psi.
    h: tamaño de paso de la segunda derivada; por defecto H_OPT.

    La normalización de psi no importa: aparece arriba y abajo, y se
    cancela.
    """

    # diff2_central espera una función de una sola variable, f(x), pero
    # psi necesita también los parámetros. psi_de_x "fija" parametros
    # y deja solo x como variable: es un adaptador entre las dos
    # interfaces.
    def psi_de_x(y):
        return psi(y, parametros)

    return -0.5 * HBAR2_SOBRE_M * diff2_central(psi_de_x, x, h) / psi(x, parametros)


###############################################
# Resultados exactos, para comparar
###############################################


def energia_cinetica_local_oscilador_exacta(x: float, parametros: dict) -> float:
    """T_L exacta del oscilador armónico (solo para alfa = 1, cuando
    psi es una eigenfunción). Como T + V = E:

        T_L(x) = E_n - V(x) = (n + 1/2) hbar w - (1/2) m w^2 x^2

    Con hbar = m = 1, w = m_w_hbar.
    """
    n = parametros["n"]
    w = parametros["m_w_hbar"]
    return (n + 0.5) * w - 0.5 * w**2 * x**2


def energia_caja_exacta(parametros: dict) -> float:
    """Energía de la onda plana en la caja periódica (no hay potencial,
    así que es toda cinética, y no depende de x):

        E = (hbar^2 / 2m) (2 pi n / L)^2
    """
    k = 2.0 * PI * parametros["n"] / parametros["L"]
    return 0.5 * HBAR2_SOBRE_M * k**2


###############################################
# Prueba: T_L en un punto, para varios h
###############################################


def probar_energia_cinetica():
    """Compara T_L numérica contra la exacta, en x = 1, para el
    oscilador (n = 100) y la caja (n = -2, L = 2 pi), con varios h."""
    x = 1.0
    parametros_oscilador = {"n": 100, "m_w_hbar": 1.0, "alfa": 1.0}
    parametros_caja = {"n": -2, "L": 2 * PI}

    exacta_oscilador = energia_cinetica_local_oscilador_exacta(x, parametros_oscilador)
    exacta_caja = energia_caja_exacta(parametros_caja)

    print(f"Energía cinética local en x = {x}")
    print(f"  oscilador (n=100): exacta = {exacta_oscilador}")
    print(f"  caja (n=-2, L=2pi): exacta = {exacta_caja}\n")
    print(f"{'h':>9}{'oscilador':>22}{'error':>10}{'caja (parte real)':>22}{'error':>10}{'parte imaginaria':>18}")

    hs = [10.0**-i for i in range(1, 6)] + [H_OPT]
    for h in hs:
        oscilador = energia_cinetica_local(psi_oscilador, x, parametros_oscilador, h)
        caja = energia_cinetica_local(psi_caja, x, parametros_caja, h)
        print(
            f"{h:9.1e}{oscilador:22.16f}{error_relativo(oscilador, exacta_oscilador):10.1e}"
            f"{caja.real:22.16f}{error_relativo(caja.real, exacta_caja):10.1e}{caja.imag:18.1e}"
        )
    print(f"(el último renglón es H_OPT = (768 EPS)^(1/4) = {H_OPT:.1e})")

    # Lo que se ve:
    # - Al reducir h, las dos aproximaciones mejoran (error de
    #   truncamiento O(h^2): cada factor de 10 en h baja el error ~100
    #   veces), hasta que el error de redondeo las alcanza (con
    #   h = 1e-5 ya empeoran).
    # - Para el oscilador con n = 100, psi tiene 100 nodos y oscila
    #   rápido (ver graficar_densidad.gp): con h = 0.1 el error es de
    #   ~4%, mucho peor que en la caja. h tiene que ser mucho más chica
    #   que la escala en la que cambia psi.
    # - Por lo mismo, H_OPT no es el mejor h para el oscilador: la
    #   fórmula de juguete supone f^(4)/f ~ 1, pero para n = 100 es
    #   mucho mayor, y el h óptimo real es más chico (~1e-4).
    # - La caja da un número complejo, aunque la energía exacta es
    #   real: la parte imaginaria es puro error numérico. Es muy chica
    #   para h grande (la resta cancela casi exacto) y crece al reducir
    #   h, cuando domina el redondeo; en la práctica se descarta.


###############################################
# T_L(x) del oscilador, en toda una región
###############################################


def guardar_energia_cinetica_oscilador(ruta, n: int = 3, limite: float = 5.0, puntos: int = 400):
    """Guarda "x  T_L numérica  T_L exacta" para el oscilador con
    número cuántico n, en [-limite, limite].

    Los puntos se toman en el centro de cada intervalo
    (x_i = -limite + (i + 1/2) dx), para no caer exactamente en un
    nodo de psi (por ejemplo, x = 0 para n impar), donde T_L no está
    definida.
    """
    parametros = {"n": n, "m_w_hbar": 1.0, "alfa": 1.0}
    dx = 2.0 * limite / puntos
    with open(ruta, "w") as archivo:
        archivo.write(f"# oscilador armonico, n = {n}, hbar = m = w = 1\n")
        archivo.write("# x  T_L_numerica  T_L_exacta\n")
        for i in range(puntos):
            x = -limite + (i + 0.5) * dx
            numerica = energia_cinetica_local(psi_oscilador, x, parametros)
            exacta = energia_cinetica_local_oscilador_exacta(x, parametros)
            archivo.write(f"{x} {numerica} {exacta}\n")


###############################################
# Principio de correspondencia
###############################################


def densidad_clasica(x: float, amplitud: float):
    """Densidad de probabilidad clásica del oscilador con amplitud x0:

        P_c(x) = 1 / (pi sqrt(x0^2 - x^2)),   para |x| < x0.

    Fuera de [-x0, x0] (región clásicamente prohibida) regresa None.
    """
    if abs(x) >= amplitud:
        return None
    return 1.0 / (PI * raiz_cuadrada(amplitud**2 - x**2))


def guardar_densidad_oscilador(ruta, n: int = 100, puntos: int = 2000):
    """Guarda "x  |psi_n(x)|^2  P_c(x)" para comparar la densidad de
    probabilidad cuántica con la clásica de la misma energía.

    La amplitud clásica sale de igualar (1/2) m w^2 x0^2 con
    (n + 1/2) hbar w: x0 = sqrt((2n + 1) hbar / (m w)). Donde P_c no
    está definida se escribe "?", que gnuplot interpreta como dato
    faltante.
    """
    parametros = {"n": n, "m_w_hbar": 1.0, "alfa": 1.0}
    amplitud = raiz_cuadrada((2 * n + 1) / parametros["m_w_hbar"])
    limite = 1.2 * amplitud
    dx = 2.0 * limite / puntos
    with open(ruta, "w") as archivo:
        archivo.write(f"# oscilador armonico, n = {n}, amplitud clasica x0 = {amplitud}\n")
        archivo.write("# x  densidad_cuantica  densidad_clasica\n")
        for i in range(puntos + 1):
            x = -limite + i * dx
            cuantica = psi_oscilador(x, parametros) ** 2
            clasica = densidad_clasica(x, amplitud)
            archivo.write(f"{x} {cuantica} {'?' if clasica is None else clasica}\n")


if __name__ == "__main__":
    probar_energia_cinetica()

    CARPETA_DATOS.mkdir(exist_ok=True)
    ruta = CARPETA_DATOS / "energia_cinetica_oscilador.dat"
    guardar_energia_cinetica_oscilador(ruta)
    print(f"\nT_L(x) del oscilador (n=3) guardada en {ruta.name}")

    ruta = CARPETA_DATOS / "densidad_oscilador.dat"
    guardar_densidad_oscilador(ruta)
    print(f"Densidades cuántica y clásica (n=100) guardadas en {ruta.name}")
