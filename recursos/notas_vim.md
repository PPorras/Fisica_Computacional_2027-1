# Notas de vim

`vim` es un editor de texto que se maneja casi todo con el teclado y
que está disponible en prácticamente cualquier sistema Linux (incluidos
servidores remotos a los que te conectas por `ssh`, donde no tienes un
editor gráfico). Vale la pena saber lo básico aunque uses otro editor
para tu trabajo diario.

## 1. La idea central: modos

A diferencia de un editor normal, en vim el teclado se comporta
distinto según el **modo** en el que estás. Los tres que vas a usar
todo el tiempo:

| Modo | Para qué sirve | Cómo entrar |
|---|---|---|
| **Normal** | moverte, borrar, copiar, pegar (es el modo por defecto) | `Esc` (desde cualquier otro modo) |
| **Inserción** | escribir texto, como en cualquier editor | `i`, `a`, `o`, `O` (ver abajo) |
| **Visual** | seleccionar texto para copiar/borrar/indentar | `v`, `V`, `Ctrl+v` |
| **Línea de comando** | comandos como guardar, salir, buscar y reemplazar | `:` |

> La confusión más común de quien empieza con vim es tratar de escribir
> estando en modo normal (y que "coma" el teclado como comandos). Si
> algo raro pasa, presiona `Esc` para regresar a modo normal y respira.

## 2. Abrir, guardar y salir

```bash
vim archivo.py        # abre (o crea) el archivo
```

Estos comandos se escriben en modo normal, empezando con `:`:

| Comando | Qué hace |
|---|---|
| `:w` | Guarda (*write*) |
| `:q` | Sale (*quit*), solo si no hay cambios sin guardar |
| `:wq` o `:x` | Guarda y sale |
| `:q!` | Sale sin guardar, descarta cambios |
| `:w nombre.py` | Guarda una copia con otro nombre |

## 3. Moverte en modo normal

| Tecla | Movimiento |
|---|---|
| `h` `j` `k` `l` | izquierda, abajo, arriba, derecha (las flechas también funcionan) |
| `w` / `b` | siguiente palabra / palabra anterior |
| `0` / `$` | inicio / fin de la línea |
| `gg` / `G` | inicio / fin del archivo |
| `:N` | ir a la línea `N`, p. ej. `:42` |
| `Ctrl+f` / `Ctrl+b` | avanzar / retroceder una pantalla |

## 4. Editar

Para escribir necesitas pasar a modo inserción; `Esc` te regresa a
normal:

| Tecla | Qué hace |
|---|---|
| `i` | insertar antes del cursor |
| `a` | insertar después del cursor |
| `o` / `O` | abrir una línea nueva abajo / arriba y entrar en inserción |
| `x` | borrar el caracter bajo el cursor |
| `dd` | borrar (cortar) la línea completa |
| `dw` | borrar la palabra desde el cursor |
| `cw` | cambiar la palabra (borra y entra en inserción) |
| `yy` | copiar (*yank*) la línea |
| `p` / `P` | pegar después / antes del cursor |
| `u` | deshacer |
| `Ctrl+r` | rehacer |
| `.` | repetir el último cambio |

Muchos de estos comandos aceptan un número antes, por ejemplo `3dd`
borra 3 líneas y `5j` baja 5 líneas.

## 5. Modo visual (seleccionar)

| Tecla | Qué selecciona |
|---|---|
| `v` | selección normal, caracter por caracter |
| `V` | selección de líneas completas |
| `Ctrl+v` | selección en bloque (columna), útil para editar varias líneas a la vez |

Una vez seleccionado, usa `d` para borrar, `y` para copiar, o `>`/`<`
para indentar/desindentar.

## 6. Buscar y reemplazar

| Comando | Qué hace |
|---|---|
| `/patron` | buscar hacia adelante |
| `?patron` | buscar hacia atrás |
| `n` / `N` | repetir la búsqueda en la misma dirección / en la contraria |
| `:%s/viejo/nuevo/g` | reemplazar todas las ocurrencias de "viejo" por "nuevo" en todo el archivo |
| `:%s/viejo/nuevo/gc` | igual, pero pide confirmación en cada ocurrencia |

## 7. Personalizar vim (`.vimrc`)

La configuración de vim vive en `~/.vimrc`. En [`vimconfig.txt`](vimconfig.txt)
tienes un ejemplo comentado con opciones típicas para trabajar con
código: números de línea, resaltado de sintaxis, usar espacios en vez
de tabs, corrector ortográfico, etc. Para usarlo:

```bash
cp recursos/vimconfig.txt ~/.vimrc
```

## 8. Por qué te sirve en este curso

- Cuando trabajes en un servidor remoto por `ssh` (ver
  [`notas_comandos_linux.md`](notas_comandos_linux.md)), normalmente
  no vas a tener VS Code ni un editor gráfico: `vim` (o `nano`, más
  simple) va a estar ahí.
- Para editar rápido un script `.py` o un `.gp` de gnuplot sin salir
  de la terminal.
- `:%s/.../.../g` es muy útil para cambiar, por ejemplo, el nombre de
  una variable en todo un script.
