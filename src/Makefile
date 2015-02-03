GCC = gcc
GCCFLAGS = -g -O2 
all: fastdate

fastdate: util.c fastdate.c maps.c
	$(GCC) $(GCCFLAGS) $+ -o $@

clean:
	rm qualign
