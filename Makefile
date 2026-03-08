###############################################################################
# Copyright (c) 2026 Ivan Pizhenko. All rights reserved.
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT ANY WARRANTIES OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
# OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH
# THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
###############################################################################

# Makefile for the matrix test program

TARGET:=mtest
SRC:=mtest.cpp
OBJ:=$(SRC:.cpp=.o)
DEP:=$(OBJ:.o=.d)

CXX:=g++
LD:=g++
CXXFLAGS:=-std=gnu++23 -pedantic -Wall -Wextra -Werror -fmax-errors=5 -MMD -MP
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
