# Entorno virtual y el paquete `fiscomp`

Estas notas explican cómo crear y usar el entorno virtual del curso,
y cómo instalar en él el paquete `fiscomp` (el código que vamos
construyendo juntos en `fiscomp/`) para poder hacer `import fiscomp`
desde cualquier práctica, unidad o notebook sin trucos de rutas.

## 1. ¿Qué es un entorno virtual y por qué lo usamos?

Un entorno virtual es una copia aislada de Python con sus propios
paquetes instalados, separada del Python del sistema. Sirve para que:

- Lo que instalemos para el curso no choque con otras cosas que
  tengas instaladas (o con las que necesite el propio sistema
  operativo).
- Todos en el curso tengamos exactamente las mismas versiones,
  evitando el clásico "en mi máquina sí funciona".
- Puedas borrar todo el entorno (`rm -rf .venv`) y volver a crearlo
  desde cero sin miedo a romper nada más.

Vive en la carpeta `.venv/` en la raíz del repositorio, y **no se
sube a git** (ya está listado en [`.gitignore`](../.gitignore)) —
cada quien lo crea una vez en su propia máquina.

## 2. Crear el entorno virtual (una sola vez por máquina)

Desde la raíz del repositorio (la carpeta que contiene `fiscomp/`,
`practicas/`, `pyproject.toml`, etc.):

```bash
python3 -m venv .venv
```

Esto crea la carpeta `.venv/` con una copia de Python adentro. No
hace falta volver a correr este comando salvo que borres `.venv/` o
quieras empezar de cero.

## 3. Activar el entorno virtual

| Sistema | Comando |
|---|---|
| Linux / macOS (bash o zsh) | `source .venv/bin/activate` |
| Windows (PowerShell) | `.venv\Scripts\Activate.ps1` |
| Windows (cmd.exe) | `.venv\Scripts\activate.bat` |

Cuando está activado, el prompt de la terminal cambia y muestra
`(.venv)` al inicio:

```
$ source .venv/bin/activate
(.venv) $
```

Eso es lo que indica que los comandos `python3` y `pip` de ahora en
adelante son los del entorno virtual, no los del sistema.

> **Hay que activarlo en cada terminal nueva.** Activar el entorno no
> es permanente: si cierras la terminal o abres una nueva, hay que
> volver a correr `source .venv/bin/activate` (no hay que volver a
> *crearlo*, solo a activarlo).

## 4. Instalar `fiscomp` en modo editable

Con el entorno activado, y desde la raíz del repositorio:

```bash
pip install -e .
```

El `-e` es por *editable*: instala el paquete `fiscomp` "apuntando"
directamente a la carpeta `fiscomp/` del repositorio, en vez de
copiarlo. Esto quiere decir que:

- Cualquier archivo `.py` que agreguemos o modifiquemos dentro de
  `fiscomp/` se refleja de inmediato, sin volver a instalar nada.
- Después de instalarlo, `import fiscomp` funciona desde cualquier
  carpeta (una práctica, una unidad, un notebook), no solo desde la
  raíz del repositorio.

Solo hay que correr `pip install -e .` una vez por máquina (a menos
que borres y recrees `.venv/`).

## 5. Comprobar que quedó bien

```bash
python3 -c "import fiscomp; print(fiscomp.__file__)"
```

Debe imprimir una ruta dentro de tu copia del repositorio, algo como
`/ruta/a/Cursos/fiscomp/__init__.py`.

## 6. Desactivar el entorno virtual

```bash
deactivate
```

Con esto el prompt vuelve a la normalidad (desaparece el `(.venv)`) y
`python3`/`pip` vuelven a apuntar al Python del sistema. No borra
nada — el entorno sigue ahí, listo para activarse de nuevo la
próxima vez con `source .venv/bin/activate`.

## 7. Resumen rápido

```bash
# Una sola vez por máquina
python3 -m venv .venv
source .venv/bin/activate
pip install -e .

# Cada vez que abras una terminal nueva para trabajar en el curso
source .venv/bin/activate
...
deactivate    # al terminar, opcional
```

> **Consejo:** si algo se ve raro (`ModuleNotFoundError: No module
> named 'fiscomp'`, o corre con una versión de Python que no
> esperabas), lo primero que hay que revisar es si el entorno está
> activado — busca el `(.venv)` al inicio del prompt.
