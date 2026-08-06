CC = gcc

# Flags
CFLAGS = -D_GNU_SOURCE -Wall -Wextra -std=c11 -O2 -Iinclude -Icsparse
LDFLAGS = -lgsl -lopenblas -lm

SRCDIR = src
INCDIR = include
OBJDIR = build

EXEC = project

SRC = $(wildcard $(SRCDIR)/*.c)
OBJ = $(SRC:$(SRCDIR)/%.c=$(OBJDIR)/%.o)

CSPARSE_SRC = csparse/csparse.c
CSPARSE_OBJ = $(OBJDIR)/csparse.o

all: $(EXEC)

$(EXEC): $(OBJ) $(CSPARSE_OBJ)
	$(CC) $(OBJ) $(CSPARSE_OBJ) -o $@ $(LDFLAGS)

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJDIR)/csparse.o: $(CSPARSE_SRC)
	mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(OBJDIR) $(EXEC)
