"""Utilidades compartidas entre las unidades del curso.

Uso típico desde un script de una unidad, por ejemplo
unidades/02_funciones_y_metodo_euler/metodo_euler.py:

    import sys
    from pathlib import Path
    sys.path.append(str(Path(__file__).resolve().parents[2]))  # raíz del curso

    from comun.simulacion import guardar_csv

    guardar_csv("datos/euler_dt_0.1.csv", ["t", "y"], filas)
"""

import csv
from pathlib import Path


def guardar_csv(ruta_relativa, encabezados, filas):
    """Guarda filas de datos en un CSV, relativo al script que llama.

    ruta_relativa: por ejemplo "datos/euler_dt_0.1.csv"
    encabezados: lista de nombres de columna, p.ej. ["t", "posicion"]
    filas: lista de tuplas/listas con los valores de cada renglón
    """
    ruta = Path(ruta_relativa)
    ruta.parent.mkdir(parents=True, exist_ok=True)

    with open(ruta, "w", newline="") as archivo_csv:
        escritor = csv.writer(archivo_csv)
        escritor.writerow(encabezados)
        escritor.writerows(filas)

    print(f"Datos guardados en {ruta.resolve()}")
