objects = recon.o excom.o max.o
CFLAGS = -c -Wall -O3

all : recon

recon : $(objects)
	g++ -o recon $(objects) -Wall -O3
recon.o : src/ccrecon.cc src/structs.hh
	g++ -c src/ccrecon.cc -o recon.o $(CFLAGS)
excom.o : src/excom.cc src/structs.hh
	g++ -c src/excom.cc -o excom.o -lm $(CFLAGS)
max.o : src/max.cc src/structs.hh
	g++ -c src/max.cc -o max.o -lm $(CFLAGS)

clean :
	rm -f $(objects)
