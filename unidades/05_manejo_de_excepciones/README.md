# Unidad 05 — Manejo de excepciones

## Temas
- **Python (Tutorial oficial, capítulo 8 — Errores y excepciones):**
  - Errores de sintaxis vs. excepciones
  - Manejo de excepciones: `try` / `except`, capturar varios tipos de
    excepción, la cláusula `else`
  - Lanzar excepciones con `raise`
  - Encadenamiento de excepciones (`raise ... from ...`)
  - Excepciones definidas por el usuario (heredando de `Exception`)
  - Acciones de limpieza definidas: `finally`
  - Acciones de limpieza predefinidas: el bloque `with` (ya visto en
    la [unidad 04](../04_archivos/), ahora desde el punto de vista de
    manejo de errores)
- **Física computacional:** distinguir errores de programación (que
  hay que corregir) de condiciones esperables en tiempo de ejecución
  (datos faltantes o mal capturados, archivos que no existen,
  mediciones físicamente imposibles); validar entradas de un
  "experimento" y señalar condiciones inválidas con excepciones
  propias, en vez de dejar que el programa truene con un traceback
  críptico o, peor, que siga corriendo con datos corruptos.

## Contenido
Código de ejemplo en [`excepciones.py`](excepciones.py).
