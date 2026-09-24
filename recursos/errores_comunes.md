# Errores comunes

Recopilación de errores frecuentes en el curso, con su causa y cómo
resolverlos. Si te topas con alguno de estos, revisa primero aquí
antes de preguntar.

## `ModuleNotFoundError: No module named 'fiscomp'`

### El error

Al correr un script de una práctica o unidad, por ejemplo:

```
unidades/06_aritmetica_punto_flotante via 🐍 v3.13.7
🕒 10:32 sh❯ ./aritmetica_punto_flotante.py
Traceback (most recent call last):
  File "/home/.../unidades/06_aritmetica_punto_flotante/./aritmetica_punto_flotante.py", line 11, in <module>
    from fiscomp.precision_numerica import EPS, error_relativo
ModuleNotFoundError: No module named 'fiscomp'
```

### Por qué pasa

El paquete `fiscomp` solo está instalado dentro del entorno virtual
del curso (`.venv/`), en modo editable. Si el entorno **no está
activado** en esa terminal, el comando `python3` (o el `#!/usr/bin/env
python3` del script) usa el Python del sistema en vez del de `.venv`,
y ese Python del sistema no conoce `fiscomp`.

El ícono 🐍 que muestra la versión de Python en el prompt (starship u
otro) **no** indica que el entorno esté activado — solo indica qué
versión de Python hay disponible. Lo que sí indica que está activado
es el prefijo `(.venv)` al inicio del prompt.

### Cómo se arregla

Activa el entorno virtual antes de correr el script:

```bash
source .venv/bin/activate   # Linux/macOS
./aritmetica_punto_flotante.py
```

En Windows: `.venv\Scripts\Activate.ps1` (PowerShell) o
`.venv\Scripts\activate.bat` (cmd.exe).

Después de activarlo, el prompt debe mostrar `(.venv)` al inicio, y
`python3 -c "import fiscomp"` ya no debe dar error.

Para más detalle sobre crear/activar el entorno e instalar `fiscomp`
en modo editable, ver
[`recursos/notas_entorno_virtual.md`](notas_entorno_virtual.md).
