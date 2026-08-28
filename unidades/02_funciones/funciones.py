#!/usr/bin/env python3

###############################################
# Definiendo funciones 
###############################################

def energia_cinetica(masa, velocidad):
    """Calcula la energía cinética (J) de un objeto de cierta masa y velocidad."""
    return 0.5 * masa * velocidad**2

print(f"E_cinetica = {energia_cinetica(2.0, 3.0)} J")

# Una función puede llamar a otra. Aquí reutilizamos el ejemplo del 
# rebote de la pelota de la unidad anterior, pero como una función 
# recursiva que regresa la lista de alturas alcanzadas en cada bote.
def alturas_de_rebote(altura, coef_restitucion, altura_minima):
    """Regresa la lista de alturas de una pelota que rebota hasta que
    la altura cae por debajo de altura_minima."""
    nueva_altura = altura * coef_restitucion
    if nueva_altura < altura_minima:
        return [nueva_altura]
    return [nueva_altura] + alturas_de_rebote(nueva_altura, coef_restitucion, altura_minima)

print(f"Alturas de rebote: {alturas_de_rebote(10.0, 0.6, 0.5)}")

###############################################
# Más sobre definición de funciones
###############################################

#################################################
# Valores por default
#################################################

def energia_potencial(masa, altura, gravedad=9.8):
    """gravedad tiene un valor por default: el de la Tierra al nivel del mar."""
    return masa * gravedad * altura

print(f"En la Tierra: {energia_potencial(2.0, 10.0)} J")        # usa gravedad=9.8
print(f"En la Luna:   {energia_potencial(2.0, 10.0, 1.62)} J")  # gravedad de la Luna

# Los valores por default se evalúan una sola vez, al definir
# la función. Por eso no conviene usar una lista (u otro objeto mutable)
# como valor por default: todas las llamadas compartirían la misma lista.
def registrar_medicion(valor, historial=None):
    """historial=None y luego se crea una lista nueva en cada llamada."""
    if historial is None:
        historial = []
    historial.append(valor)
    return historial

print(registrar_medicion(9.81))
print(registrar_medicion(9.79))   # no arrastra la medición anterior

#################################################
# Argumentos por palabra clave (Keyword Arguments)
#################################################

def describir_particula(nombre, masa, carga=0.0, espin=0.5):
    """Describe una partícula: nombre, masa (kg), carga (C) y espín."""
    return f"{nombre}: masa={masa} kg, carga={carga} C, espín={espin}"

# Por posición
print(describir_particula("electrón", 9.11e-31, -1.602e-19, 0.5))
# Por palabra clave, en cualquier orden
print(describir_particula(masa=1.67e-27, nombre="protón", carga=1.602e-19))
# Mezclando ambas (las posicionales siempre van primero)
print(describir_particula("neutrón", 1.675e-27, espin=0.5))

# Llamadas que fallarían (se dejan comentadas):
# describir_particula()                          # faltan argumentos requeridos
# describir_particula(masa=1.0, "fotón")          # posicional después de keyword: SyntaxError
# describir_particula("neutrino", 0, 0, 0.5, 0)   # sobran argumentos: TypeError

#################################################
# Parámetros especiales (Special parameters)
#################################################

### Positional-or-Keyword Arguments
# Es el caso por default: describir_particula (arriba) usa este tipo de
# parámetros, por eso se pudo llamar tanto por posición como por nombre.

### Positional-Only Parameters
def redondear_medicion(valor, cifras, /):
    """valor y cifras son positional-only (nota el '/'): sus nombres no
    son parte de la interfaz, solo su orden."""
    return round(valor, cifras)

print(f"Medición redondeada: {redondear_medicion(9.80665, 2)}")
# redondear_medicion(valor=9.80665, cifras=2)   # TypeError: son positional-only

### Keyword-Only Arguments
def simular_caida(altura_inicial, *, gravedad=9.8, dt=0.1):
    """gravedad y dt son keyword-only (nota el '*'): hay que nombrarlos
    al llamar la función, para que la llamada sea inequívoca.

    Integra la caída con el método de Euler (velocidad y altura se
    actualizan paso a paso), hasta que la altura llega a cero."""
    tiempo = 0.0
    altura = altura_inicial
    velocidad = 0.0
    while altura > 0.0:
        velocidad += gravedad * dt
        altura -= velocidad * dt
        tiempo += dt
    return tiempo

print(f"Tiempo de caída en la Luna: {simular_caida(20.0, gravedad=1.62, dt=0.05):.2f} s")
# simular_caida(20.0, 1.62, 0.05)   # TypeError: gravedad y dt deben darse por nombre

### Function Examples
def medir_particula(id_particula, /, masa, *, unidades="kg"):
    """id_particula es positional-only, masa es positional-or-keyword,
    y unidades es keyword-only."""
    return f"Partícula {id_particula}: masa = {masa} {unidades}"

print(medir_particula(1, 9.11e-31))
print(medir_particula(2, masa=1.67e-27, unidades="kg"))

# - Usa positional-only ('/') cuando el nombre del parámetro no aporta
#   nada (p. ej. valores genéricos, como en redondear_medicion).
# - Usa keyword-only ('*') en parámetros opcionales que conviene nombrar
#   explícitamente para que la llamada sea clara (p. ej. gravedad, dt).
# - Deja positional-or-keyword (el caso por default) para el resto.

#################################################
# Listas de argumentos arbitrarias (Arbitrary Argument Lists)
#################################################

def fuerza_neta(*fuerzas):
    """Suma (en 1D) un número arbitrario de fuerzas, en N."""
    return sum(fuerzas)

print(f"Fuerza neta: {fuerza_neta(10.0, -3.0, 5.5)} N")
print(f"Fuerza neta: {fuerza_neta(1.0, 2.0, 3.0, 4.0, 5.0)} N")

#################################################
# Desempaquetado de listas de argumentos (Unpacking Argument Lists)
#################################################

componentes_de_fuerza = [3.0, -4.0, 0.0]
print(f"Fuerza neta (desempaquetando lista): {fuerza_neta(*componentes_de_fuerza)} N")

parametros_caida_luna = {"gravedad": 1.62, "dt": 0.05}
print(f"Caída en la Luna (desempaquetando dict): {simular_caida(20.0, **parametros_caida_luna):.2f} s")

#################################################
# Expresiones lambda (Lambda Expressions)
#################################################

energia_cinetica_lambda = lambda masa, velocidad: 0.5 * masa * velocidad**2
print(f"E_cinetica (lambda) = {energia_cinetica_lambda(2.0, 3.0)} J")

# Uso típico de lambda: como 'key' al ordenar mediciones
particulas = [("electrón", 9.11e-31), ("protón", 1.67e-27), ("neutrón", 1.675e-27)]
particulas_por_masa = sorted(particulas, key=lambda p: p[1])
print(f"Partículas ordenadas por masa: {particulas_por_masa}")

#################################################
# Cadenas de documentación (Documentation Strings)
#################################################

def energia_en_reposo(masa, velocidad_luz=3e8):
    """Calcula la energía en reposo de una partícula (E = m c^2).

    Sigue la convención de docstrings de PEP 257: una línea de resumen,
    una línea en blanco, y luego una descripción más detallada.

    Parámetros
    ----------
    masa : float
        Masa en reposo de la partícula, en kg.
    velocidad_luz : float, opcional
        Velocidad de la luz en m/s (por default, ~3e8 m/s).

    Regresa
    -------
    float
        La energía en reposo, en Joules.
    """
    return masa * velocidad_luz**2

print(f"E_reposo del electrón = {energia_en_reposo(9.11e-31):.3e} J")
print(energia_en_reposo.__doc__)
