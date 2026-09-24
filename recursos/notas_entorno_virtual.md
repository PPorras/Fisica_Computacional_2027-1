# Entorno virtual y el paquete `fiscomp`

El código que se reutiliza en todo el curso (por ejemplo
`error_relativo`, `EPS`, `seno`, `coseno`…) vive en el paquete
`fiscomp/`, en la raíz del repositorio. Para que cualquier script de
una unidad o práctica pueda hacer

```python
from fiscomp.precision_numerica import EPS, error_relativo
```

sin importar en qué carpeta esté, `fiscomp` se instala dentro de un
**entorno virtual**. Estas notas explican qué es eso, cómo crearlo y
cómo usarlo en el día a día.

## 1. ¿Qué es un entorno virtual?

Es una carpeta (en el curso, `.venv/`) que contiene **su propia copia
de Python y sus propios paquetes**, separados del Python del sistema.
Sirve para:

- Instalar paquetes (`fiscomp`, `matplotlib`, `numpy`…) sin permisos
  de administrador y sin tocar el Python del sistema operativo.
- Que cada proyecto tenga sus propias versiones de paquetes, sin
  chocar con otros proyectos.
- Poder borrarlo y volver a crearlo desde cero si algo se descompone.

La carpeta `.venv/` **no se sube a git** (está en `.gitignore`): cada
quien crea la suya en su computadora.

## 2. Requisitos

Python 3.9 o más reciente. Para verificarlo:

```bash
python3 --version        # Linux / macOS
python --version         # Windows
```

Si no lo tienes:

- **Linux:** normalmente ya viene instalado. En Ubuntu/Debian quizá
  falte el módulo de entornos virtuales: `sudo apt install python3-venv`.
- **macOS:** `brew install python`, o el instalador de
  <https://www.python.org/downloads/>.
- **Windows:** el instalador de <https://www.python.org/downloads/>.
  Durante la instalación marca **"Add python.exe to PATH"**.

## 3. Crear el entorno (una sola vez)

Desde la **raíz del repositorio** (la carpeta que contiene
`README.md` y `pyproject.toml`):

```bash
python3 -m venv .venv          # Linux / macOS
python -m venv .venv           # Windows
```

Esto crea la carpeta `.venv/`. No hace falta repetirlo: el entorno se
queda ahí hasta que lo borres.

## 4. Activar el entorno (cada vez que abras una terminal)

| Sistema                  | Comando                          |
|--------------------------|----------------------------------|
| Linux / macOS            | `source .venv/bin/activate`      |
| Windows (PowerShell)     | `.venv\Scripts\Activate.ps1`     |
| Windows (cmd.exe)        | `.venv\Scripts\activate.bat`     |

Cuando está activado, el prompt empieza con `(.venv)`:

```
(.venv) usuario@maquina:~/Fisica_Computacional_2027-1$
```

Mientras está activado, los comandos `python3`, `python` y `pip`
usan el Python del entorno, no el del sistema.

> **Importante:** la activación solo dura en **esa** terminal. Si
> abres una terminal nueva (o una pestaña nueva), tienes que volver a
> activarlo. Es la causa más común de errores en el curso (ver
> [`errores_comunes.md`](errores_comunes.md)).

> **Windows:** si PowerShell dice que "la ejecución de scripts está
> deshabilitada en este sistema", corre esto una vez y vuelve a
> intentarlo:
>
> ```powershell
> Set-ExecutionPolicy -Scope CurrentUser RemoteSigned
> ```

Se puede activar desde cualquier carpeta dando la ruta hacia `.venv`;
por ejemplo, si estás en `unidades/08_diferencias_finitas/`:

```bash
source ../../.venv/bin/activate
```

### Desactivar

```bash
deactivate
```

(o simplemente cierra la terminal).

## 5. Instalar `fiscomp` en modo editable (una sola vez)

Con el entorno activado y desde la raíz del repositorio:

```bash
python3 -m pip install -e .
```

- El `.` significa "el proyecto de esta carpeta": pip lee
  `pyproject.toml` y encuentra el paquete `fiscomp/`.
- El `-e` (**editable**) hace que pip no copie el código, sino que
  apunte a la carpeta `fiscomp/` del repositorio. Así, cuando agregues
  o cambies una función en `fiscomp/` (o hagas `git pull` y lleguen
  funciones nuevas), el cambio se ve de inmediato, **sin reinstalar**.

## 6. Verificar que todo quedó bien

Con el entorno activado:

```bash
which python3        # Linux/macOS: debe terminar en .venv/bin/python3
where python         # Windows: la primera línea debe estar en .venv\Scripts\
python3 -c "import fiscomp; print(fiscomp.__file__)"
```

El último comando debe imprimir una ruta que termine en
`fiscomp/__init__.py` **dentro de tu repositorio**. Si da
`ModuleNotFoundError`, el entorno no está activado o falta el paso 5.

Para ver todo lo instalado en el entorno:

```bash
python3 -m pip list
```

## 7. Instalar otros paquetes

Algunos scripts del curso necesitan paquetes extra. Se instalan igual,
con el entorno activado:

```bash
python3 -m pip install matplotlib
```

Se recomienda `python3 -m pip` en lugar de solo `pip`, porque
garantiza que el paquete se instala en el mismo Python con el que vas
a correr los scripts. Ver [`notas_matplotlib.md`](notas_matplotlib.md)
para el caso de matplotlib.

## 8. Resumen: el día a día

La primera vez:

```bash
cd ruta/al/repositorio
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -e .
```

Cada vez que abras una terminal para trabajar en el curso:

```bash
cd ruta/al/repositorio
source .venv/bin/activate
```

(En Windows, cambia `python3` por `python` y
`source .venv/bin/activate` por `.venv\Scripts\Activate.ps1`.)

## 9. Problemas comunes

### `ModuleNotFoundError: No module named 'fiscomp'`

El entorno no está activado en esa terminal (revisa que el prompt
muestre `(.venv)`), o no se hizo el paso 5. Detalle completo en
[`errores_comunes.md`](errores_comunes.md).

### `error: externally-managed-environment`

Intentaste instalar con `pip` **fuera** del entorno virtual, en el
Python del sistema. No uses `--break-system-packages`: activa el
entorno y vuelve a instalar.

### `The virtual environment was not created successfully because ensurepip is not available`

En Ubuntu/Debian falta el módulo de entornos virtuales:

```bash
sudo apt install python3-venv
```

Borra la carpeta `.venv/` que quedó a medias y vuelve a crearla.

### Moví o renombré la carpeta del repositorio y el entorno dejó de funcionar

Los entornos virtuales guardan rutas absolutas, así que no sobreviven
a que se mueva la carpeta. Solución: bórralo y créalo de nuevo (pasos
3 a 5). No se pierde nada de tu código; solo hay que reinstalar los
paquetes.

```bash
rm -rf .venv                               # Linux / macOS
Remove-Item -Recurse -Force .venv          # Windows (PowerShell)
```

### El editor (VS Code, etc.) no encuentra `fiscomp`

El editor también tiene que usar el Python del entorno. En VS Code:
`Ctrl+Shift+P` → **Python: Select Interpreter** → elige el que dice
`.venv`.
