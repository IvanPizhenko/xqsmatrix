# Makefile for the matrix test program

TARGET:=mtest
SRC:=mtest.cpp
OBJ:=$(SRC:.cpp=.o)
DEP:=$(OBJ:.o=.d)

CXX:=g++
LD:=g++
CXXFLAGS:=-std=gnu++11 -pedantic -Wall -Wextra -Werror -fmax-errors=5
LDFLAGS:=
LIBS:=-lm

ifeq ($(DEBUG),1)
CXXFLAGS+=-O0 -g3 -DDEBUG -D_DEBUG
LDFLAGS+=-g3
else
CXXFLAGS+=-O2
endif

all: $(TARGET)

clean:
	-rm -f $(TARGET)
	-rm -f *.d
	-rm -f *.o

debug:
	$(MAKE) DEBUG=1

-include $(DEP)

$(TARGET): $(OBJ)
	$(LD) -o $@ $(LDFLAGS) $^ $(LIBS)
