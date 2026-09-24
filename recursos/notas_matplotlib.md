# Instalar matplotlib en el entorno virtual del curso

**matplotlib** es la biblioteca más usada en Python para hacer
gráficas. No viene incluida en el entorno base del curso: el
`README.md` solo instala `fiscomp`. Algunos scripts de las unidades la
usan (por ejemplo,
[`unidades/08_diferencias_finitas/graficar_derivada.py`](../unidades/08_diferencias_finitas/graficar_derivada.py)).
Si la corres sin tenerla instalada, verás un mensaje como este:

```
Este script necesita matplotlib, que no está instalado en este entorno virtual. ...
```

o, en un script sin ese aviso:

```
ModuleNotFoundError: No module named 'matplotlib'
```

matplotlib se instala **una sola vez** dentro del entorno virtual
`.venv/`. Después queda disponible para todos los scripts del curso.

## 1. Ir a la raíz del repositorio

Abre una terminal y entra a la carpeta del repositorio del curso (la
que contiene `README.md`, `pyproject.toml` y `.venv/`):

```bash
cd ruta/al/repositorio
ls        # debes ver README.md, pyproject.toml, fiscomp/, unidades/, ...
```

Si no existe la carpeta `.venv/`, primero crea el entorno siguiendo el
`README.md`.

## 2. Activar el entorno virtual

| Sistema                  | Comando                          |
|--------------------------|----------------------------------|
| Linux / macOS            | `source .venv/bin/activate`      |
| Windows (PowerShell)     | `.venv\Scripts\Activate.ps1`     |
| Windows (cmd.exe)        | `.venv\Scripts\activate.bat`     |

Cuando está activado, el prompt empieza con `(.venv)`.

> **Windows:** si PowerShell dice que "la ejecución de scripts está
> deshabilitada", corre esto una vez y vuelve a intentarlo:
>
> ```powershell
> Set-ExecutionPolicy -Scope CurrentUser RemoteSigned
> ```

## 3. Instalar matplotlib

Con el entorno activado (el comando es igual en los tres sistemas):

```bash
python3 -m pip install matplotlib
```

En Windows, si `python3` no existe, usa `python`:
`python -m pip install matplotlib`.

Se recomienda escribir `python3 -m pip` en lugar de solo `pip`, porque
así pip instala el paquete en el **mismo** Python que vas a usar para
correr los scripts. Si solo escribes `pip`, a veces se usa el de otro
Python del sistema.

## 4. Verificar la instalación

```bash
python3 -c "import matplotlib; print(matplotlib.__version__)"
```

Si imprime una versión (por ejemplo `3.10.6`), ya quedó. Ahora sí:

```bash
cd unidades/08_diferencias_finitas
python3 graficar_derivada.py
```

Debe abrirse una ventana con la gráfica.

## 5. Problemas comunes

### Instalé matplotlib, pero sigue diciendo que no está

Casi siempre es porque el entorno **no está activado en esa terminal**.
Cada terminal nueva empieza sin entorno: hay que activarlo otra vez
(paso 2). Revisa que el prompt muestre `(.venv)` y que Python sea el
del entorno:

```bash
which python3          # Linux/macOS: debe terminar en .venv/bin/python3
where python           # Windows: la primera línea debe estar en .venv\Scripts\
```

Ver también [`errores_comunes.md`](errores_comunes.md), en el error
`ModuleNotFoundError: No module named 'fiscomp'`: la causa es la misma.

### El script termina pero no aparece ninguna ventana (Linux)

matplotlib necesita un "backend" gráfico para abrir ventanas. El más
común usa `tkinter`, que en algunas distribuciones viene en un paquete
aparte del sistema (no se instala con pip):

```bash
sudo dnf install python3-tkinter     # Fedora / RHEL
sudo apt install python3-tk          # Ubuntu / Debian
sudo pacman -S tk                    # Arch
```

Para comprobarlo: `python3 -c "import tkinter"` no debe dar error.

Si trabajas en un servidor o en WSL sin interfaz gráfica, en vez de
`plt.show()` guarda la figura en un archivo y ábrela aparte:

```python
plt.savefig("grafica.png", dpi=150)
```

### `error: externally-managed-environment`

Aparece si intentas instalar con `pip` **fuera** del entorno virtual,
en el Python del sistema. No lo fuerces: activa el entorno (paso 2) y
vuelve a instalar.

## 6. Un primer ejemplo

Para comprobar que todo funciona, guarda esto en un archivo `prueba.py`
y córrelo con el entorno activado:

```python
import math

import matplotlib.pyplot as plt

x = [i * 0.1 for i in range(100)]
y = [math.sin(v) for v in x]

plt.plot(x, y, label="sin(x)")
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.grid(True)
plt.show()
```

Documentación y galería de ejemplos: <https://matplotlib.org/stable/>.
