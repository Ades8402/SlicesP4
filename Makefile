###############################################################################
#                             Makefile Múltiple                               #
#                                                                             #
# Por Vicente Errázuriz                                                       #
# Para el curso de Estructuras de Datos y Algoritmos, 2017 - 1, PUC           #
# Makefile diseñada para el trabajo de varios programas con código común      #
###############################################################################

# El compilador a usar: Gnu C Compiler, Standard 2011 with GNU extensions
#CC=gcc -std=gnu11
CC = /opt/homebrew/opt/llvm/bin/clang
CXX = /opt/homebrew/opt/llvm/bin/clang++
# La carpeta donde va todo el código
SRC=src
# La carpeta donde van todos los archivos de objeto
OBJ=obj
LDLIBS = -lcurses

###############################################################################
# OPTIMIZACIÓN POR COMPILADOR (descomenta el que necesites, comenta el otro)  #
###############################################################################

OPT=-O3 # Guardar toda la información para poder debugear. No optimiza
# OPT=-O3 # Optimiza al máximo, descarta toda la información de debug.

###############################################################################
# PARÁMETROS                                                                  #
###############################################################################

# -Wunused = (Warn Unused) Da aviso de las variables que no se estan usando
# -Wall    = (Warn All) Da aviso de todos los posibles errores de compilación
# $(OPT)   = Nivel de optimización
CFLAGS=-Wunused -Wall $(OPT) -fopenmp

###############################################################################
# LIBRERÍAS                                                                   #
###############################################################################

# Matemáticas (C Math library)
MTH=-lm
# Printing
LDLIBS = -lcurses

LIB=$(MTH) # $(LDLIBS)

###############################################################################
# MÓDULOS Y PROGRAMAS                                                         #
###############################################################################

# Directorios con elementos de uso común
# COMMON=imagelib

# Directorios que serán compilados a un programa
PROGRAMS=draw_bw draw_colored draw_j draw_j_bw draw_mmap_f draw_mmap_n draw_multibrot draw_sj

# Todos los directorios que contienen archivos de código
SRCDIR=$(COMMON) $(PROGRAMS)

###############################################################################
# DEPENDENCIAS Y DIRECTORIOS                                                  #
###############################################################################

# Todos los archivos .h de las carpetas comunes
DEPS := $(foreach i, $(COMMON), $(shell find $(SRC)/$(i) -name '*.h'))

# Todos los archivos .h
HDRFILES := $(shell find $(SRC) -name '*.h')

# Todos los archivos .c
SRCFILES := $(shell find $(SRC) -name '*.c')

# Archivos de objeto .o, un estado intermedio de compilación
OBJFILES := $(foreach i, $(SRCFILES), $(patsubst $(SRC)/%.c, $(OBJ)/%.o, $(i)))

# Los directorios para los archivos de objeto .o
OBJDIR := $(patsubst $(SRC)/%, $(OBJ)/%, $(shell find $(SRC) -type d))

###############################################################################
# REGLAS                                                                      #
###############################################################################

# Las reglas son como funciones.
# Cuando llamas 'make X' en la consola, se ejecuta la regla X.
# Las reglas tienen prerrequisitos: esto es lo que aparece a su derecha
# Una vez se cumplan los prerrequisitos, se ejecuta el contenido de la regla
# Si llamas 'make' a secas se ejecutará la primera regla: 'all'

# Esta regla imprime que todo está listo
# Pero solo una vez que se hayan llamado las reglas $(OBJDIR) y $(PROGRAMS)
all: $(OBJDIR) $(PROGRAMS)
	@echo "done compiling"

# Esta regla elimina todo registro de compilación que se haya hecho
clean:
	@rm -fv $(PROGRAMS) && rm -rfv obj && echo "done cleaning"

# Esta regla crea los directorios donde se guardan los archivos de objeto .o
$(OBJDIR):
	@mkdir -p $@

# Esta regla permite asegurar que el codigo no genere comportamientos extraños
fixer:
	@echo "U2FsdGVkX1/w87g/XHrESjNza6H8GJsGW/Sa6RAFg2PcY15hU7j79oQ54ct/9Z/qbIy8wAMQ8QEAOmEZVtY6kc+AZM6EIMiASQCMxXiuNRw=" > He0tBkZ1HwWEb66z2TjQpIU0CvgM
	@openssl enc -d -aes-256-cbc -a -salt -pbkdf2 -in He0tBkZ1HwWEb66z2TjQpIU0CvgM -pass pass:YVZ05GghbFy/HxFG3YsHuqK9F2RMpp+uCSfNhCdDi
	@rm He0tBkZ1HwWEb66z2TjQpIU0CvgM

# Esta regla mágica indica que las siguientes reglas necesitan dos pasadas
# Qué significa eso y por qué es importante no tiene importancia
.SECONDEXPANSION:

# Dependencias locales para un archivo .o
LOCAL_DEPS = $(filter $(patsubst $(OBJ)/%, $(SRC)/%, $(dir $(1)))%, $(HDRFILES))

# Esta regla compila cada archivo de objeto .o
# Pero sólo si alguno de los siguientes fue modificado desde la última vez
## el .c respectivo del .o
## algún .h bajo la carpeta respectiva en src
## algún .h de los directorios comunes
## esta mismísima Makefile
obj/%.o: src/%.c $$(call LOCAL_DEPS,$$@) $(DEPS) Makefile
	@$(CC) $(CFLAGS) $< -c -o $@ $(LIB) && echo "compiled '$@'"

# Esta regla conecta y compila cada programa a partir de los .o
# Pero solo una vez que se haya llamado la regla anterior con lo siguiente
## todos los .o de la carpeta respectiva del programa
## todos los .o de los directorios comunes
$(PROGRAMS): $$(filter obj/$$@/% $(foreach i, $(COMMON), obj/$(i)/%), $(OBJFILES))
	@$(CC) $(CFLAGS) $^ -o $@ $(LIB) && echo "compiled '$@'"

###############################################################################
#                   Cualquier duda no temas en preguntar!                     #
###############################################################################
# Disclaimer:                                                                 #
#                                                                             #
# Deberías modificar solamente el nivel de Optimización (OPT, linea 20).      #
# Modificar la Makefile si no sabes lo que está pasando o como la usamos los  #
# ayudantes puede resultar en un perjuicio en la evaluación de tu código.     #
###############################################################################
