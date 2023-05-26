# -fPIC process independent code (for shared libraries)
# -Wall to warn about all my bugs :-)
# -std=c99 because I define ints in my loops and use complex numbers
# -march=native gives me an extra 10%
CC = gcc -Wall -fPIC -shared -std=c99 -march=native -O2
SRC = src/*.c
HDR = src/fof.h

# Find the cpython suffix for the library
CPYTHON_SUFFIX := $(shell python -c 'import sysconfig; print(sysconfig.get_config_var("EXT_SUFFIX"))')
LIB = build/libhfof$(CPYTHON_SUFFIX)

all: build ${LIB}
build:
	mkdir -p build

${LIB}: ${SRC} ${HDR} Makefile
	${CC} -o ${LIB} ${SRC}
clean:
	rm -rf build

