ARCH ?= $(shell arch)

CC = g++
CDEBUG = -ggdb -g# -pg
OPTM = -O3
WARN = -W -Wall
CFLAGS = $(WARN) -I.

BIN = $(HOME)/bin/$(ARCH)
RM = rm -rf

D = DATA

ALLSRC = makeBlocks

all: $D $(ALLSRC)

$D:
	-mkdir $D
	ln -s ../config.file $D/config.file

%: %.cpp
	$(CC) $(OPTM) $(addsuffix .cpp, $(basename $@)) -o $@

.PHONY: clean
clean:
	$(RM) *.dSYM $(ALLSRC)
	
