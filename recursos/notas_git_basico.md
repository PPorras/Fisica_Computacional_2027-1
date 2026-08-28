# Comandos básicos de git

Estas notas resumen el flujo de trabajo mínimo de git que se necesita
para llevar el control de versiones del código durante el curso:
guardar el progreso, ver qué cambió, y (más adelante) trabajar con
ramas y con GitHub.

## 1. Configuración inicial (una sola vez por máquina)

```bash
git config --global user.name "Tu Nombre"
git config --global user.email "tu_correo@ejemplo.com"
```

Esto es lo que aparece como autor en cada commit que hagas.

## 2. Crear o entrar a un repositorio

| Comando | Qué hace | Ejemplo |
|---|---|---|
| `git init` | Convierte la carpeta actual en un repositorio de git (nuevo, vacío) | `git init` |
| `git clone <url>` | Descarga una copia de un repositorio existente | ver [notas de cómo clonar de GitHub](notas_clonar_github.md) |

## 3. El flujo básico: modificar → agregar → guardar

Git trabaja en tres "zonas": tu carpeta de trabajo, el *staging area*
(o índice) y el historial de commits.

```
carpeta de trabajo  --git add-->  staging area  --git commit-->  historial
```

| Comando | Qué hace | Ejemplo |
|---|---|---|
| `git status` | Muestra qué archivos cambiaron, cuáles están en staging y cuáles no | `git status` |
| `git add` | Manda un archivo (o varios) al staging area | `git add simulacion.py`, `git add .` (todo lo que cambió) |
| `git commit -m "mensaje"` | Guarda lo que está en staging como un nuevo punto en el historial | `git commit -m "agrega integrador de Euler"` |
| `git diff` | Muestra los cambios que aún no están en staging | `git diff` |
| `git diff --staged` | Muestra los cambios que ya están en staging pero no se han commiteado | `git diff --staged` |

> **Consejo:** usa `git status` seguido, sobre todo antes de un
> `git add .` — así ves exactamente qué vas a agregar y evitas subir
> archivos que no querías (como `.venv/` o archivos de datos grandes).

## 4. Ver el historial

| Comando | Qué hace | Ejemplo |
|---|---|---|
| `git log` | Muestra el historial de commits | `git log` |
| `git log --oneline` | Una línea por commit, más fácil de leer | `git log --oneline -n 10` |
| `git show <commit>` | Muestra los cambios de un commit específico | `git show a1b2c3d` |

## 5. Deshacer cosas

| Comando | Qué hace | Ejemplo |
|---|---|---|
| `git restore <archivo>` | Descarta cambios sin guardar en un archivo (vuelve a la última versión del *commit*) | `git restore simulacion.py` |
| `git restore --staged <archivo>` | Saca un archivo del staging area (sin perder los cambios) | `git restore --staged simulacion.py` |
| `git commit --amend` | Corrige el mensaje (o el contenido) del último commit | `git commit --amend -m "mensaje corregido"` |

> **Cuidado:** `git restore <archivo>` borra cambios sin guardar y no
> se puede deshacer. Si tienes dudas, primero haz `git status` para
> confirmar qué se va a perder.

## 6. Ramas (branches)

Las ramas sirven para probar cosas sin afectar tu código principal.

| Comando | Qué hace | Ejemplo |
|---|---|---|
| `git branch` | Lista las ramas que existen (`*` marca en la que estás) | `git branch` |
| `git branch <nombre>` | Crea una rama nueva | `git branch integrador-rk4` |
| `git switch <nombre>` | Cambia a otra rama | `git switch integrador-rk4` |
| `git switch -c <nombre>` | Crea una rama y cambia a ella en un solo paso | `git switch -c integrador-rk4` |
| `git merge <nombre>` | Junta los cambios de otra rama a la rama actual | `git merge integrador-rk4` (estando en `main`) |

## 7. Convenciones de nombres: ramas y commits

Seguir una convención hace que el historial sea fácil de leer, tanto
para ti como para quien revise tu trabajo (o para ti mismo dentro de
unos meses).

### Ramas

Formato: `<tipo>/<descripcion-corta-en-kebab-case>`

| Tipo | Cuándo usarlo | Ejemplo |
|---|---|---|
| `feat/` | Agregas algo nuevo (un ejemplo, una práctica, una unidad) | `feat/integrador-rk4` |
| `fix/` | Corriges un error | `fix/signo-gravedad` |
| `docs/` | Cambios de documentación o notas | `docs/notas-git-basico` |
| `refactor/` | Reorganizas código sin cambiar su comportamiento | `refactor/separa-funciones-fisica` |

> La descripción va en minúsculas, sin espacios ni acentos, separando
> palabras con guiones (`-`), como en `feat/integrador-rk4`.

### Mensajes de commit (Conventional Commits)

Formato: `<tipo>(<alcance>): <descripción en presente, minúsculas>`

| Tipo | Cuándo usarlo |
|---|---|
| `feat` | Agregas una funcionalidad, ejemplo o práctica nueva |
| `fix` | Corriges un error |
| `docs` | Cambios solo en documentación o notas |
| `refactor` | Reorganizas código existente sin cambiar su comportamiento |
| `style` | Cambios de formato (espacios, comas, etc.) que no afectan el código |
| `test` | Agregas o corriges pruebas |

El `<alcance>` es opcional: suele ser la unidad o el tema al que
pertenece el cambio (por ejemplo `class`, `tipos-de-datos`,
`control-flujo`).

Ejemplos:

```
feat(class): add boolean type and comparison operators example
feat(class): add if/elif/else, while, for-range, break and continue examples
docs(recursos): add branch and commit naming conventions
fix(unidad-01): corrige signo de la gravedad en caída libre
```

> **Consejo:** el mensaje debe completar la frase "Si se aplica, este
> commit ...". Por ejemplo: "Si se aplica, este commit **agrega
> ejemplo de tipo booleano**", no "Si se aplica, este commit
> **arreglos varios**".

## 8. Trabajar con un repositorio remoto (GitHub)

| Comando | Qué hace | Ejemplo |
|---|---|---|
| `git remote -v` | Muestra a qué URL remota está conectado el repo | `git remote -v` |
| `git push` | Sube tus commits al remoto | `git push origin main` |
| `git pull` | Trae y combina los cambios del remoto | `git pull` |
| `git fetch` | Trae los cambios del remoto sin combinarlos todavía | `git fetch` |

Si `git pull` te dice que hay un **conflicto** (dos cambios distintos
sobre la misma línea), git marca el archivo con algo así:

```
<<<<<<< HEAD
tu versión del código
=======
la versión que venía del remoto
>>>>>>> origin/main
```

Edita el archivo a mano para dejar la versión correcta (borrando esas
marcas `<<<<<<<`, `=======`, `>>>>>>>`), y luego:

```bash
git add archivo_con_conflicto.py
git commit
```

## 9. Ejemplo integrador

Un flujo típico para guardar el avance de una práctica y subirlo a
GitHub:

```bash
git status                              # ver qué cambió
git add integrador_euler.py             # agregar el archivo nuevo
git commit -m "implementa método de Euler para la práctica 2"
git push origin main                    # subirlo al repositorio remoto
```

## 10. Archivos que no quieres versionar: `.gitignore`

Crea un archivo `.gitignore` en la raíz del repositorio y lista ahí
lo que git debe ignorar siempre (entornos virtuales, datos generados,
archivos temporales). El repositorio de este curso ya trae uno
([`.gitignore`](../.gitignore)) que ignora el entorno virtual
(`.venv/`), la caché de Python (`__pycache__/`, `*.pyc`) y las
salidas de simulaciones y gráficas de cada unidad (`datos/`,
`figuras/`), entre otras cosas.

Si un archivo ya estaba versionado *antes* de agregarlo al
`.gitignore`, agregarlo ahí no basta — hay que quitarlo del
seguimiento explícitamente (sin borrarlo del disco):

```bash
git rm --cached archivo_que_ya_no_quieres_versionar.txt
```
