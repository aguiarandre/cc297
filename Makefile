#
# Makefile
#
# CC 297 - CFD
# Projeto 1
#


# compiler to use
CC = clang

# Flags para passar ao compilar (usando -Wall, -Werror : Nível hard!)
CFLAGS = -ggdb3 -O0 -Qunused-arguments -std=c11 -Wall -Werror

# nome do executável.
EXE = projeto1

# lista de arquivos 'Header', separados por 'espaços'. 
HDRS = mesh.h helpers.h definitions.h solvers.h

# lista de 'Libraries', separadas por espaço.
# cada qual com acompanhada pelo prefixo -l 
LIBS = -lm 

# lista de arquivos fonte, separados por espaço.
SRCS = mesh.c helpers.c solvers.c projeto1.c

# lista automática de arquivos .o (object files) 
OBJS = $(SRCS:.c=.o)

# target principal (que vai ser executado em default pelo comando 'make')
$(EXE): $(OBJS) $(HDRS) Makefile
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LIBS)

# dependências 
$(OBJS): $(HDRS) Makefile

# apagar tudo!
clean:
	rm -f core $(EXE) *.o
