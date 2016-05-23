# Copyright (C) 2015 Tomas Flouri
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Contact: Tomas Flouri <Tomas.Flouri@h-its.org>,
# Heidelberg Institute for Theoretical Studies,
# Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany

# Profiling
#PROFILING=-g -pg
PROFILING=-g

# Compiler warnings
WARN=-Wall -Wsign-compare

CC = gcc
CFLAGS = -g -O3 $(WARN) 
LINKFLAGS=$(PROFILING)
LIBS=-lm -lpthread

BISON = bison
FLEX = flex

PROG=fastdate

all: $(PROG)

OBJS=util.o maps.o fastdate.o tree.o dp.o gamma.o bd.o parse_newick.o \
     lex_newick.o arch.o exp.o ln.o uni.o norm.o nodeprior.o lca.o optimize.o \
     lbfgsb/lbfgsb.o lbfgsb/linesearch.o lbfgsb/linpack.o lbfgsb/miniCBLAS.o \
     lbfgsb/subalgorithms.o sample.o parse_prior.o lex_prior.o

$(PROG): $(OBJS)
	$(CC) -Wall $(LINKFLAGS) $+ -o $@ $(LIBS)

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

parse_%.c: parse_%.y
	$(BISON) -p $*_ -d -o $@ $<

lex_%.c: lex_%.l
	$(FLEX) -P $*_ -o $@ $<

clean:
	rm -f *~ $(OBJS) gmon.out $(PROG) parse_newick.c lex_newick.c parse_newick.h parse_prior.h
