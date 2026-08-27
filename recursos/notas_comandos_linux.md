# Notas de Linux para Física Computacional

Estas notas resumen los comandos básicos de la terminal que ya vimos en
clase, y agregan algunos más que son muy útiles cuando se trabaja con
simulaciones, scripts de Python y datos (por ejemplo, en un clúster o
servidor remoto).

## 1. Los comandos que ya vimos

### Navegación

| Comando | Qué hace | Ejemplo |
|---|---|---|
| `pwd` | Muestra el directorio en el que estás (*print working directory*) | `pwd` |
| `cd` | Cambia de directorio | `cd datos/`, `cd ..` (subir un nivel), `cd ~` (ir al home), `cd -` (volver al directorio anterior) |
| `ls` | Lista el contenido de un directorio | `ls`, `ls -l` (formato largo), `ls -a` (incluye ocultos), `ls -lh` (tamaños legibles: KB/MB), `ls -la` |
| `whoami` | Muestra el usuario con el que estás conectado | `whoami` |

### Archivos y directorios

| Comando | Qué hace | Ejemplo |
|---|---|---|
| `touch` | Crea un archivo vacío (o actualiza su fecha de modificación) | `touch resultados.txt` |
| `cp` | Copia archivos o directorios | `cp datos.txt respaldo.txt`, `cp -r carpeta/ copia/` (recursivo, para directorios) |
| `mv` | Mueve o renombra archivos/directorios | `mv datos.txt datos_finales.txt` |
| `rm` | Elimina archivos | `rm resultados.txt`, `rm -r carpeta/` (recursivo), `rm -i` (pide confirmación) |
| `cat` | Muestra el contenido completo de un archivo | `cat script.py` |

> **Cuidado con `rm`:** no manda los archivos a una papelera, los borra
> directamente. No existe un `Ctrl+Z` para `rm`. Si tienes dudas usa
> `rm -i` para que te pregunte antes de borrar cada archivo.

### Permisos

| Comando | Qué hace | Ejemplo |
|---|---|---|
| `chmod` | Cambia permisos de un archivo | `chmod +x simulacion.py` (lo hace ejecutable), `chmod 755 script.sh` |

Recordatorio rápido de los tres bloques de permisos (`usuario`,
`grupo`, `otros`) y sus valores:

```
r (leer)    = 4
w (escribir) = 2
x (ejecutar) = 1

chmod 754 archivo
      │││
      ││└─ otros:  4       (r--)
      │└── grupo:  5 = 4+1 (r-x)
      └─── usuario: 7 = 4+2+1 (rwx)
```

## 2. Comandos nuevos recomendados

### Crear y explorar directorios

| Comando | Qué hace | Ejemplo |
|---|---|---|
| `mkdir` | Crea un directorio | `mkdir simulaciones`, `mkdir -p proyecto/datos/2026` (crea toda la ruta de una vez) |
| `tree` | Muestra el árbol de directorios (puede no venir instalado) | `tree -L 2` (solo 2 niveles) |

### Ver contenido de archivos sin abrir un editor

`cat` funciona bien para archivos cortos, pero para archivos largos
(por ejemplo, la salida de una simulación con miles de líneas) conviene
usar:

| Comando | Qué hace | Ejemplo |
|---|---|---|
| `less` | Muestra el archivo página por página (`q` para salir, `/palabra` para buscar) | `less resultados.dat` |
| `head` | Muestra las primeras líneas de un archivo | `head -n 20 datos.csv` |
| `tail` | Muestra las últimas líneas de un archivo | `tail -n 20 log.txt`, `tail -f log.txt` (sigue el archivo en vivo, útil mientras corre una simulación) |
| `wc` | Cuenta líneas, palabras o caracteres | `wc -l datos.csv` (número de líneas, útil para saber cuántos puntos de datos hay) |

### Buscar

| Comando | Qué hace | Ejemplo |
|---|---|---|
| `grep` | Busca texto dentro de archivos | `grep "error" log.txt`, `grep -r "def " .` (busca recursivamente en todos los archivos) |
| `find` | Busca archivos por nombre, tipo, etc. | `find . -name "*.py"` (todos los `.py` desde aquí), `find . -mtime -1` (modificados en el último día) |

### Comodines (wildcards)

Muy útiles para trabajar con varios archivos a la vez:

```bash
ls *.py          # todos los archivos que terminan en .py
rm datos_*.csv   # todos los que empiezan con datos_
cp *.dat respaldo/
```

### Redirección y tuberías (pipes)

Fundamentales para guardar resultados y encadenar comandos:

| Símbolo | Qué hace | Ejemplo |
|---|---|---|
| `>` | Redirige la salida a un archivo (lo sobrescribe) | `python3 simulacion.py > resultados.txt` |
| `>>` | Redirige agregando al final del archivo | `echo "corrida 2" >> bitacora.txt` |
| `<` | Toma la entrada de un archivo | `python3 script.py < entrada.txt` |
| `\|` | Manda la salida de un comando como entrada de otro | `cat log.txt \| grep "error"`, `ls -l \| wc -l` |

### Procesos (útil para simulaciones que tardan)

| Comando | Qué hace | Ejemplo |
|---|---|---|
| `&` | Corre un comando en segundo plano | `python3 simulacion_larga.py &` |
| `jobs` | Lista los procesos en segundo plano de la sesión actual | `jobs` |
| `fg` / `bg` | Trae un proceso al primer plano / lo manda a segundo plano | `fg %1` |
| `ps` | Lista procesos en ejecución | `ps aux` |
| `top` / `htop` | Monitor de procesos en tiempo real (uso de CPU/memoria) | `top` (`q` para salir) |
| `kill` | Termina un proceso por su PID | `kill 12345`, `kill -9 12345` (forzado) |
| `nohup` | Deja un proceso corriendo aunque cierres la terminal | `nohup python3 simulacion.py &` |

### Otros útiles

| Comando | Qué hace | Ejemplo |
|---|---|---|
| `man` | Manual de un comando | `man ls` (`q` para salir) |
| `comando --help` | Ayuda rápida de un comando | `ls --help` |
| `history` | Muestra el historial de comandos usados | `history`, o con flecha ↑ / `Ctrl+R` para buscar |
| `echo` | Imprime texto o el valor de una variable | `echo "Hola"`, `echo $PATH` |
| `which` | Dice en dónde está instalado un programa | `which python3` |
| `df -h` | Espacio disponible en disco (legible) | `df -h` |
| `du -sh carpeta/` | Tamaño que ocupa una carpeta | `du -sh datos/` |
| `tar` | Empaqueta/comprime archivos | `tar -czvf datos.tar.gz datos/` (comprimir), `tar -xzvf datos.tar.gz` (descomprimir) |

### Para trabajar en un servidor o clúster remoto

Muy relevante si en algún momento del curso corren simulaciones en un
servidor remoto en vez de su computadora:

| Comando | Qué hace | Ejemplo |
|---|---|---|
| `ssh` | Conectarse a otra máquina | `ssh usuario@servidor.unam.mx` |
| `scp` | Copiar archivos entre tu máquina y el servidor | `scp resultados.txt usuario@servidor:~/datos/` |

### Editor de texto en terminal

Para editar archivos sin salir de la terminal (por ejemplo, al estar
conectado por `ssh`):

| Comando | Qué hace |
|---|---|
| `nano archivo.py` | Editor sencillo: `Ctrl+O` guarda, `Ctrl+X` sale |
| `vim archivo.py` | Editor más potente pero con curva de aprendizaje mayor |

## 3. Ejemplo integrador

Un flujo típico para correr una simulación larga y revisar sus
resultados sin quedarte pegado a la terminal:

```bash
mkdir -p simulaciones/corrida_01
cd simulaciones/corrida_01
chmod +x simulacion.py
nohup python3 simulacion.py > salida.log 2>&1 &
tail -f salida.log        # Ctrl+C para dejar de seguir el archivo (no mata el proceso)
grep "error" salida.log   # revisar si algo falló
```
