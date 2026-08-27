# Convenciones de nombres en Python (PEP 8)

Python no obliga a seguir un estilo de nombres, pero la comunidad sigue casi
universalmente la guía de estilo oficial **PEP 8**
(<https://peps.python.org/pep-0008/>). Usar estas convenciones hace que el
código sea legible para cualquier persona que lo revise, incluyendo a tu
"yo" del futuro.

## Resumen rápido

| Elemento                     | Convención                          | Ejemplo                          |
|-------------------------------|--------------------------------------|-----------------------------------|
| Variables                    | `snake_case`                        | `velocidad_luz`, `masa_electron` |
| Constantes                   | `UPPER_SNAKE_CASE`                  | `PI`, `VELOCIDAD_LUZ`            |
| Funciones                    | `snake_case`                        | `calcular_energia()`             |
| Métodos                      | `snake_case`                        | `self.calcular_area()`           |
| Clases                       | `PascalCase` (o `CapWords`)         | `class Particula:`               |
| Módulos / archivos `.py`     | `snake_case`, minúsculas, cortos    | `tipos_de_datos_python.py`       |
| Paquetes                     | minúsculas, sin guion bajo          | `numpy`, `scipy`                 |
| Atributo/método "protegido"  | un guion bajo al inicio             | `_valor_interno`                 |
| Atributo/método "privado"    | dos guiones bajos al inicio         | `__valor_oculto`                 |
| Métodos/atributos especiales | dos guiones bajos antes y después   | `__init__`, `__str__`            |
| Booleanos (estilo, no PEP 8) | prefijo `es_` / `tiene_`            | `es_estable`, `tiene_carga`      |

## Detalle por elemento

### Variables
Minúsculas, palabras separadas por guion bajo (`snake_case`). El nombre
debe describir qué contiene la variable, no cómo se usa.

```python
tiempo_inicial = 0.0
posicion_x = 3.2
numero_particulas = 10
```

### Constantes
También en mayúsculas cuando el valor no debe cambiar durante la
ejecución del programa. Es una convención (Python no las hace realmente
inmutables), pero comunica intención.

```python
VELOCIDAD_LUZ = 3e8       # m/s
CONSTANTE_PLANCK = 6.626e-34  # J·s
```

### Funciones y métodos
`snake_case`, usualmente un verbo que indica la acción que realizan.
Los métodos siguen la misma regla; el primer parámetro de un método de
instancia se llama `self` por convención.

```python
def calcular_energia(masa, velocidad):
    return masa * velocidad**2

class Particula:
    def calcular_momento(self):
        return self.masa * self.velocidad
```

### Clases
`PascalCase` (también llamado `CapWords`): cada palabra empieza con
mayúscula y no se usan guiones bajos.

```python
class ParticulaCargada:
    ...

class OsciladorArmonico:
    ...
```

### Archivos / módulos
Igual que las variables: minúsculas y `snake_case`, nombres cortos y
descriptivos. Evita mayúsculas y espacios en los nombres de archivo
(`.py`).

```
tipos_de_datos_python.py   # correcto
IntroPython.py              # evitar (PascalCase en un archivo)
Intro Python.py             # evitar (espacios)
```

### Atributos y métodos "privados"
Python no tiene privacidad real, pero por convención:

- `_nombre`: uso interno, "no lo toques desde fuera" (protegido).
- `__nombre`: activa *name mangling*, dificulta el acceso accidental
  desde fuera de la clase (privado).
- `__nombre__`: reservado para métodos/atributos especiales del
  lenguaje ("dunder", de *double underscore*), por ejemplo
  `__init__`, `__str__`, `__len__`. No inventes nombres con este
  formato para tu propio código.

### Otras recomendaciones de PEP 8
- Evita nombres de una sola letra ambigua como `l`, `O`, `I` (se
  confunden con `1` y `0`).
- Prefiere nombres descriptivos sobre abreviaciones poco claras:
  `temperatura` es mejor que `temp` si hay riesgo de confusión.
- No mezcles estilos: no uses `camelCase` (propio de otros lenguajes
  como Java o JavaScript) para variables o funciones en Python.
