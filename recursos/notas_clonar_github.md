# Cómo clonar el repositorio del curso desde GitHub

Estas notas explican cómo obtener una copia local del repositorio del
curso (o de cualquier repositorio de GitHub) usando `git clone`, y qué
hacer si necesitas actualizarla más adelante.

## 1. Requisitos previos

Verifica que tengas `git` instalado:

```bash
git --version
```

Si no aparece nada, instálalo según tu distribución:

```bash
sudo dnf install git      # Fedora / RHEL
sudo apt install git      # Ubuntu / Debian
sudo pacman -S git        # Arch
```

## 2. Clonar por HTTPS (la forma más sencilla)

Copia la URL del repositorio desde el botón verde **Code** en GitHub
(elige la pestaña **HTTPS**) y corre:

```bash
git clone https://github.com/usuario/nombre-del-repo.git
```

Esto crea una carpeta nueva llamada `nombre-del-repo/` con todo el
historial del proyecto adentro. Para entrar:

```bash
cd nombre-del-repo
```

Si quieres que la carpeta tenga otro nombre:

```bash
git clone https://github.com/usuario/nombre-del-repo.git otro-nombre
```

> **Repositorios privados por HTTPS:** GitHub te pedirá usuario y
> contraseña, pero ya no acepta tu contraseña normal — necesitas un
> *personal access token* (se genera en GitHub, en
> Settings → Developer settings → Personal access tokens) y lo usas
> en vez de la contraseña.

## 3. Clonar por SSH (recomendado si vas a hacer `push` seguido)

Con SSH no te pide usuario/contraseña cada vez, pero primero hay que
configurar una llave:

```bash
# 1. Generar una llave SSH (si no tienes una ya)
ssh-keygen -t ed25519 -C "tu_correo@ejemplo.com"

# 2. Copiar la llave pública
cat ~/.ssh/id_ed25519.pub
```

Pega el contenido en GitHub: Settings → SSH and GPG keys → New SSH key.

Después, clona usando la URL SSH (pestaña **SSH** del botón **Code**):

```bash
git clone git@github.com:usuario/nombre-del-repo.git
```

## 4. Verificar que quedó bien

```bash
cd nombre-del-repo
git status      # debería decir "nothing to commit, working tree clean"
git log -n 5    # los últimos 5 commits, para confirmar que sí bajó el historial
```

## 5. Mantener tu copia actualizada

Si el repositorio original recibe cambios después de que lo clonaste
(por ejemplo, subo material nuevo al curso), tráelos con:

```bash
git pull
```

Si tienes cambios propios sin guardar y `git pull` se queja, revisa
las [notas de comandos básicos de git](notas_git_basico.md), sección
de ramas y conflictos.

## 6. Crear tu propio repositorio y subirlo a GitHub

Esto es lo contrario de clonar: en vez de bajar un repositorio que ya
existe, partes de código que ya tienes en tu computadora y lo subes a
GitHub por primera vez.

### Paso 1: crear el repositorio en GitHub

En [github.com](https://github.com) da clic en el botón **New**
(o el `+` de arriba a la derecha → **New repository**). Ponle nombre,
elige si es público o privado, y **no** marques la opción de agregar
un `README` si ya tienes archivos localmente (así evitas conflictos
al hacer el primer `push`). Da clic en **Create repository**.

GitHub te muestra una URL como `https://github.com/usuario/mi-repo.git`
(o la versión SSH) — la vas a necesitar en el paso 3.

### Paso 2: convertir tu carpeta local en un repositorio de git

Párate en la carpeta con tu código y corre:

```bash
cd mi-proyecto
git init
git add .
git commit -m "primer commit"
```

(Ver las [notas de comandos básicos de git](notas_git_basico.md) si
`git add`/`git commit` no te suenan.)

### Paso 3: conectarla con GitHub y subir el código

```bash
git remote add origin https://github.com/usuario/mi-repo.git
git branch -M main
git push -u origin main
```

- `git remote add origin <url>` — le dice a tu repositorio local
  dónde está la copia remota (`origin` es solo un nombre, por
  convención).
- `git branch -M main` — se asegura de que tu rama principal se
  llame `main` (algunas versiones viejas de git usan `master` por
  default).
- `git push -u origin main` — sube tus commits y, gracias a `-u`,
  deja tu rama local conectada con `origin/main` para que después
  baste con `git push` a secas.

Si ya tenías un repositorio local (por ejemplo, uno que empezaste con
`git init` hace tiempo) y solo te faltaba conectarlo a GitHub, basta
con el paso 3.

## 7. Errores comunes

| Mensaje | Causa probable | Qué hacer |
|---|---|---|
| `fatal: repository not found` | La URL está mal escrita, o el repo es privado y no tienes acceso | Revisa la URL, o pide acceso al dueño del repositorio |
| `Permission denied (publickey)` | Estás usando la URL SSH pero no configuraste tu llave en GitHub | Sigue la sección 3, o usa la URL HTTPS en su lugar |
| `fatal: destination path '...' already exists and is not an empty directory` | Ya existe una carpeta con ese nombre | Bórrala, muévela, o clona con otro nombre (`git clone <url> otro-nombre`) |
| `error: failed to push some refs` (al hacer `git push`) | El repositorio remoto ya tiene commits que no tienes localmente (por ejemplo, creaste el repo en GitHub marcando "agregar README") | `git pull origin main --allow-unrelated-histories`, resuelve conflictos si aparecen, y vuelve a hacer `git push` |
| `remote origin already exists` (al hacer `git remote add`) | Ya habías conectado un remoto antes | Revisa cuál con `git remote -v`; si es el correcto no hace falta agregarlo de nuevo, o cámbialo con `git remote set-url origin <url>` |
