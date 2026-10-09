# Solución de la Práctica 2

Una solución posible de la [Práctica 2](../practica_02.md).

- **Parte 1, `VectorND`:** la versión de clase está en
  [`fiscomp/vectores.py`](../../../fiscomp/vectores.py).
- **Parte 2, `Matrix`:** [`matrices.py`](matrices.py), con `__add__`,
  `__sub__` y `__mul__` completos. Es una copia de
  `fiscomp/matrices.py` tal como quedó al terminar la práctica. Desde
  entonces, el `fiscomp/matrices.py` del repositorio creció: tiene la
  clase hija `Identity`, y en la Práctica 4 se le agregan las
  matrices triangulares.

Puntos a notar en la solución:

- `__sub__` no repite el recorrido elemento por elemento: es
  `self + other * -1`, reutilizando `__add__` y `__mul__`.
- `__mul__` distingue los casos con `isinstance`: primero escalar
  (`int` o `float`), luego `Matrix`, y si no es ninguno levanta
  `TypeError`. El orden importa: la última línea (`raise`) solo se
  alcanza si los dos `if` anteriores no regresaron nada.
- Con `__rmul__ = __mul__`, `2 * A` funciona igual que `A * 2`. Ojo:
  para dos matrices eso no importa, porque Python llama primero al
  `__mul__` de la de la izquierda, y `A * B` nunca llega a
  `B.__rmul__(A)`.

## Cómo correrla

Desde la raíz del repositorio, con el entorno virtual activado:

```bash
python3 practicas/practica2/solucion_practica_02/matrices.py
```

No copien este archivo encima de `fiscomp/matrices.py`: el del
repositorio ya tiene más cosas (`Identity`, y sus matrices
triangulares de la Práctica 4) que se perderían. Para comparar,
basta con ver lado a lado sus `__add__`, `__sub__` y `__mul__` y los
de aquí.
