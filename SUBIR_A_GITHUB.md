# Subir este proyecto a GitHub desde Windows

## 1. Preparar la carpeta

1. Descargá `nsf-energy_github.zip` y elegí **Extraer todo**.
2. Abrí la carpeta `nsf-energy`. Adentro deben aparecer `README.md`,
   `graficos_nsf_local.py`, `analysis`, `docs`, etc.
3. El ZIP `nsf_reference_data.zip` contiene las muestras numéricas.
   Guardalo aparte; va como archivo de una Release.

## 2. Subir el código desde el navegador

1. Entrá en <https://github.com/new> con tu cuenta.
2. Poné, por ejemplo, `nsf-energy` como nombre del repositorio.
3. Elegí **Private** para trabajar primero con acceso restringido, o
   **Public** si ya querés que cualquiera pueda verlo.
4. Crealo vacío: dejá desactivadas las opciones de agregar README,
   `.gitignore` y licencia; este paquete ya trae los dos primeros.
5. En el repositorio vacío, elegí **uploading an existing file**.
   Si ya tiene archivos, usá **Add file > Upload files**.
6. Arrastrá **el contenido descomprimido** de `nsf-energy`, conservando
   sus carpetas. Así `README.md` queda en la raíz del repositorio.
7. Escribí `Initial NSF numerical code` como mensaje y pulsá
   **Commit changes**.

Subí el código desde la copia limpia recién extraída. El navegador no
aplica `.gitignore` al elegir archivos; los resultados de ejecuciones
locales no forman parte del paquete de código.

## 3. Adjuntar los datos para reproducibilidad

1. Abrí **Releases** y elegí **Draft a new release**.
2. Creá la etiqueta `v0.1.0` apuntando a la rama `main`.
3. Poné un título, por ejemplo `Reference numerical data`.
4. En el área de archivos adjuntos arrastrá `nsf_reference_data.zip`.
5. Publicá la Release cuando quieras que esté disponible a quienes
   tengan acceso al repositorio.

El código funciona sin estos datos para generar las tablas desde el
resumen incluido. Para reanalizar las muestras originales, descargá el
asset y ejecutá los tres comandos de la sección de reproducción del README.

## 4. Actualizaciones con GitHub Desktop

Para cambios frecuentes podés usar <https://desktop.github.com/>:

1. Iniciá sesión y cloná el repositorio que acabás de crear.
2. Hacé los cambios en esa copia local.
3. En **Changes**, revisá los archivos y escribí un resumen del cambio.
4. Pulsá **Commit to main** y luego **Push origin**.

GitHub Desktop respeta `.gitignore`, que mantiene fuera de los commits
las muestras grandes, los entornos Python y los resultados generados.

## 5. Comprobar que quedó bien

- El README debe verse en la página principal.
- `analysis`, `algebra`, `docs` y `tests` deben ser carpetas navegables.
- En **Actions** aparecerá `Python checks`: se ejecutarán pruebas en
  Linux y Windows. El paquete se comprobó localmente en Linux; la
  ejecución de GitHub verificará también el entorno Windows.
- La Release debe ofrecer el archivo `nsf_reference_data.zip`.

La licencia de reutilización del código queda a elección de los autores;
este paquete no adopta una automáticamente.

Referencias oficiales consultadas:

- <https://docs.github.com/en/repositories/creating-and-managing-repositories/creating-a-new-repository>
- <https://docs.github.com/en/repositories/working-with-files/managing-files/adding-a-file-to-a-repository>
- <https://docs.github.com/en/repositories/releasing-projects-on-github/managing-releases-in-a-repository>
- <https://docs.github.com/en/desktop/overview/creating-your-first-repository-using-github-desktop>
